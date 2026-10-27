"""Descriptive annotated-gene backgrounds and prespecified sensitivity.

Coverage and sampling algorithms follow the MIT-licensed
cacao-marker-harmonization 3.0.0 coverage/{analysis,numerics}.py by Ezekiel Ahn.
This module independently rebuilds inputs and validates the frozen outputs.
"""
from collections import Counter, defaultdict
import hashlib
import csv
from contextlib import ExitStack
import json
from pathlib import Path
import numpy as np

BG_WORKBOOK = 'alignment_coverage/06_matched_gene_background.xlsx'
LENGTH_MODEL = 'chromosome_biotype_length'
DENSITY_MODEL = 'chromosome_biotype_length_gene_density'
TARGETS = ('unique_candidate_genes', 'protein_coding_candidate_genes')
SETTINGS = {
    'seed': 1729, 'replicates': 1000, 'length_bins': [3, 5, 10],
    'gene_density_bins': [2, 3, 5], 'gene_density_radius_bp': [100000, 250000, 500000],
    'primary_contrast_flank_bp': 50000,
    'percentile_definition': '100 * (n_draws_below + 0.5 * n_draws_equal) / n_draws',
    'quantile_method': 'linear',
    'sampling': 'Without replacement within chromosome/biotype/length[/gene-density] strata; one draw shared across all panels, states, and flanks.',
    'insufficient_strata': 'Block entire target/model; no relaxation or candidate omission.',
}


def rng_for(seed, label):
    h = hashlib.sha256((str(seed) + '\x1f' + label).encode('utf-8')).digest()
    return np.random.Generator(np.random.PCG64(int.from_bytes(h[:16], 'big')))


def gene_density(genes, radius):
    """Other annotated-gene midpoint counts at inclusive +/- radius endpoints."""
    if radius < 0:
        raise ValueError('Negative gene-density radius')
    result = np.zeros(len(genes), dtype=np.int64)
    groups = defaultdict(list)
    for i, g in enumerate(genes):
        groups[g['chromosome']].append(i)
    for indices in groups.values():
        mids = np.array([(genes[i]['start_bp'] + genes[i]['end_bp']) // 2 for i in indices], dtype=np.int64)
        ordered = np.sort(mids)
        result[indices] = (np.searchsorted(ordered, mids + radius, side='right')
                           - np.searchsorted(ordered, mids - radius, side='left') - 1)
    return result


def quantile_cuts(values, n_bins):
    if n_bins < 1 or not len(values):
        raise ValueError('Empty pool or invalid bin count')
    if n_bins == 1:
        return np.array([], dtype=float)
    return np.unique(np.quantile(np.asarray(values), np.arange(1, n_bins) / n_bins, method='linear'))


def nearest_distances(genes, points):
    result = np.full(len(genes), np.inf)
    groups = defaultdict(list)
    for i, g in enumerate(genes):
        groups[g['chromosome']].append(i)
    for chrom, ix in groups.items():
        p = points.get(chrom)
        if p is None or not len(p):
            continue
        if np.any(np.diff(p) <= 0):
            raise ValueError('Marker positions must be unique and sorted')
        starts = np.array([genes[i]['start_bp'] for i in ix], dtype=np.int64)
        ends = np.array([genes[i]['end_bp'] for i in ix], dtype=np.int64)
        k = np.searchsorted(p, starts, side='left')
        left, right = np.full(len(ix), np.inf), np.full(len(ix), np.inf)
        valid = k > 0
        left[valid] = starts[valid] - p[k[valid] - 1]
        valid = k < len(p)
        right[valid] = np.maximum(p[k[valid]] - ends[valid], 0)
        result[ix] = np.minimum(left, right)
    return result


def sample_strata(pools, required, rng):
    picks = []
    for key in sorted(required):
        n = required[key]
        pool = pools.get(key, np.array([], dtype=np.int64))
        if len(pool) < n:
            raise ValueError('Insufficient stratum; no relaxation or replacement')
        if n:
            picks.extend(rng.choice(pool, size=n, replace=False).tolist())
    if len(picks) != len(set(picks)):
        raise ValueError('Duplicated sampled gene')
    return np.asarray(picks, dtype=np.int64)


def empirical_percentile(draw_counts, observed_count):
    a = np.asarray(draw_counts)
    if not len(a):
        raise ValueError('Empty background distribution')
    below, equal = int(np.sum(a < observed_count)), int(np.sum(a == observed_count))
    above = len(a) - below - equal
    return dict(n_draws_below=below, n_draws_equal=equal, n_draws_above=above,
                percentile_strict=100.0 * below / len(a),
                percentile_midrank=100.0 * (below + 0.5 * equal) / len(a),
                percentile_inclusive=100.0 * (below + equal) / len(a))


def model_specs():
    out = []
    for n_length in SETTINGS['length_bins']:
        out.append(dict(model_id='L%d' % n_length, model=LENGTH_MODEL,
                        length_bins_requested=n_length, gene_density_bins_requested=None,
                        gene_density_radius_bp=None,
                        canonical=n_length == 5))
        for n_density in SETTINGS['gene_density_bins']:
            for radius in SETTINGS['gene_density_radius_bp']:
                out.append(dict(model_id='L%d_D%d_R%d' % (n_length, n_density, radius),
                                model=DENSITY_MODEL, length_bins_requested=n_length,
                                gene_density_bins_requested=n_density, gene_density_radius_bp=radius,
                                canonical=n_length == 5 and n_density == 3 and radius == 250000))
    for r in out:
        r['rng_label'] = ('annotation-background:' + r['model'] if r['canonical']
                          else 'annotation-background-sensitivity:' + r['model_id'])
    return out


def build_catalog(raw, lengths, excluded_ids):
    genes, exclusions, seen = [], [], set()
    for row in raw:
        gid = str(row.get('gene_id') or '').strip()
        if not gid or gid in seen:
            raise ValueError('Empty or duplicate annotation gene ID')
        seen.add(gid)
        chrom = str(row.get('chromosome') or '').strip()
        start, end = row.get('start_bp'), row.get('end_bp')
        biotype = str(row.get('biotype') or '').strip()
        valid = (chrom in lengths and isinstance(start, int) and isinstance(end, int)
                 and 1 <= start <= end <= lengths[chrom] and bool(biotype))
        if not valid:
            exclusions.append(dict(gene_id=gid, reason='NONPRIMARY_INVALID_INTERVAL_OR_MISSING_BIOTYPE'))
            continue
        genes.append(dict(query_id=gid, gene_id=gid, chromosome=chrom,
                          start_bp=start, end_bp=end, biotype=biotype,
                          gene_length_bp=end - start + 1, excluded_candidate=gid in excluded_ids))
    genes.sort(key=lambda g: g['gene_id'])
    return genes, exclusions


def assignments(genes, target_ix, spec, densities):
    """Cutpoints use all eligible noncandidates of each target biotype."""
    types = sorted(set(genes[i]['biotype'] for i in target_ix))
    lengths = np.asarray([g['gene_length_bp'] for g in genes])
    radius = spec['gene_density_radius_bp'] or 250000
    d_values = densities[radius]
    length_assign = np.full(len(genes), -1, dtype=np.int64)
    density_assign = np.full(len(genes), -1, dtype=np.int64)
    cuts, bin_diagnostics = [], []
    for biotype in types:
        all_ix = [i for i, g in enumerate(genes) if g['biotype'] == biotype]
        pool_ix = [i for i in all_ix if not genes[i]['excluded_candidate']]
        if not pool_ix:
            raise ValueError('No background biotype: ' + biotype)
        features = [('gene_length_bp', lengths, spec['length_bins_requested'], length_assign)]
        if spec['model'] == DENSITY_MODEL:
            features.append(('local_annotated_gene_count', d_values, spec['gene_density_bins_requested'], density_assign))
        for feature, values, requested, dest in features:
            boundary = quantile_cuts(values[pool_ix], requested)
            dest[all_ix] = np.searchsorted(boundary, values[all_ix], side='right')
            for j, v in enumerate(boundary, 1):
                cuts.append(dict(biotype=biotype, feature=feature, cut_number=j, cut_value=float(v)))
            bin_diagnostics.append(dict(biotype=biotype, feature=feature, requested_bins=requested,
                                        distinct_cutpoints=len(boundary), theoretical_bin_labels=len(boundary) + 1,
                                        occupied_noncandidate_bins=len(set(dest[pool_ix].tolist())),
                                        n_noncandidate_genes=len(pool_ix)))
    def key(i):
        g = genes[i]
        k = (g['chromosome'], g['biotype'], int(length_assign[i]))
        return k + (int(density_assign[i]),) if spec['model'] == DENSITY_MODEL else k
    required = Counter(key(i) for i in target_ix)
    pools = defaultdict(list)
    for i, g in enumerate(genes):
        if not g['excluded_candidate'] and length_assign[i] >= 0:
            pools[key(i)].append(i)
    pools = {k: np.array(v, dtype=np.int64) for k, v in pools.items()}
    strata = []
    for k, n in sorted(required.items()):
        available = len(pools.get(k, []))
        strata.append(dict(chromosome=k[0], biotype=k[1], length_bin=k[2],
                           gene_density_bin=k[3] if len(k) > 3 else None,
                           candidate_genes=n, available_noncandidate_genes=available, adequate=available >= n))
    return required, pools, strata, cuts, bin_diagnostics, length_assign, density_assign


def _numeric_close(a, b):
    if isinstance(a, (int, float, np.integer, np.floating)) and not isinstance(a, (bool, np.bool_)):
        return b is not None and bool(np.isclose(a, b, rtol=1e-12, atol=1e-12))
    return a == b


def _compare_rows(ctx, label, actual, frozen, key_fields, fields):
    a = {tuple(r[k] for k in key_fields): r for r in actual}
    b = {tuple(r[k] for k in key_fields): r for r in frozen}
    ctx.check(label + ':actual_unique', len(a), len(actual))
    ctx.check(label + ':frozen_unique', len(b), len(frozen))
    ctx.check(label + ':keys_match', set(a) == set(b), True)
    mismatches = []
    mismatch_count = 0
    max_error = 0.0
    for key in set(a).intersection(b):
        for f in fields:
            av, bv = a[key].get(f), b[key].get(f)
            if not _numeric_close(av, bv):
                mismatch_count += 1
                if len(mismatches) < 10:
                    mismatches.append(dict(key=list(key), field=f, calculated=av, frozen=bv))
            if isinstance(av, (float, int)) and isinstance(bv, (float, int)):
                max_error = max(max_error, abs(av - bv))
    ctx.check(label + ':cell_agreement', mismatch_count, 0)
    return dict(comparison=label, rows_recomputed=len(actual), rows_frozen=len(frozen),
                compared_fields=len(fields), matched_keys=len(set(a).intersection(b)),
                max_absolute_numeric_difference=max_error, mismatched_cells=mismatch_count, mismatch_examples=mismatches,
                status='PASS' if set(a) == set(b) and not mismatches and len(a) == len(actual) and len(b) == len(frozen) else 'FAIL')


def _run(ctx, stack):
    recorded_settings = json.loads((Path(__file__).resolve().parent / 'SENSITIVITY_SETTINGS.json').read_text(encoding='utf-8'))
    if recorded_settings != dict(SETTINGS, models=model_specs()):
        raise ValueError('Recorded sensitivity settings differ from implementation')
    cfg = ctx.baseline_config['coverage']
    ctx.check('background_seed_preserved', cfg['seed'], SETTINGS['seed'])
    ctx.check('background_replicates_preserved', cfg['background']['replicates'], SETTINGS['replicates'])
    windows = cfg['windows_bp']
    keys = [(r, arm) for r in cfg['comparison_resources'] for arm in ('A', 'B', 'C')]
    lengths = {str(r['chromosome']): int(r['length_bp']) for r in ctx.read('reference_mapping/04_reference_annotation.xlsx', 'sequence_inventory') if r.get('chromosome') in ctx.baseline_config['reference']['chromosomes']}
    excluded_ids = {r['gene_id'] for r in ctx.read('source_data/03_candidate_units.xlsx', 'unique_genes')}
    genes, exclusions = build_catalog(ctx.read('reference_mapping/04_reference_annotation.xlsx', 'gene_catalog'), lengths, excluded_ids)
    id_index = {g['gene_id']: i for i, g in enumerate(genes)}
    ctx.check('background_catalog_unique', len(id_index), len(genes))
    ctx.check('all_source_candidates_excluded', sum(g['excluded_candidate'] for g in genes), len(excluded_ids))
    targets = {}
    for target, sheet, expected in [(TARGETS[0], 'unique_candidate_genes', 457), (TARGETS[1], 'protein_coding_candidates', 453)]:
        rows = ctx.read('alignment_coverage/03_fixed_query_sets.xlsx', sheet)
        target_ids = {r['query_id'] for r in rows}
        ctx.check('fixed_target_size:' + target, len(target_ids), expected)
        ctx.check('fixed_target_rows_unique:' + target, len(target_ids), len(rows))
        ctx.check('target_is_candidate:' + target, target_ids <= excluded_ids, True)
        for r in rows:
            if r['query_id'] not in id_index or any(genes[id_index[r['query_id']]][f] != r[f] for f in ('chromosome', 'start_bp', 'end_bp', 'biotype')):
                raise ValueError('Target differs from annotation catalog: ' + r['query_id'])
        targets[target] = np.asarray([i for i, g in enumerate(genes) if g['gene_id'] in target_ids], dtype=np.int64)
    coordinate_lists = {k: defaultdict(list) for k in keys}
    for r in ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx', 'marker_record_ledger'):
        if r['resource'] not in cfg['comparison_resources']:
            continue
        for arm in ('A', 'B', 'C'):
            if not r['in_' + arm]:
                continue
            chrom = str(r['reference_chromosome'] if arm == 'C' else r['original_chromosome'])
            pos = int(r['reference_snp_bp'] if arm == 'C' else r['original_position_bp'])
            if chrom not in lengths or not 1 <= pos <= lengths[chrom]:
                raise ValueError('Marker outside primary chromosome bounds')
            coordinate_lists[(r['resource'], arm)][chrom].append(pos)
    points = {k: {c: np.unique(np.asarray(v, dtype=np.int64)) for c, v in p.items()} for k, p in coordinate_lists.items()}
    distances = np.stack([nearest_distances(genes, points[k]) for k in keys], axis=1)
    density_values = {r: gene_density(genes, r) for r in SETTINGS['gene_density_radius_bp']}
    covered = distances[:, :, None] <= np.asarray(windows)[None, None, :]
    frozen = {sheet: ctx.read(BG_WORKBOOK, sheet) for sheet in ('status', 'summary', 'stratum_requirements', 'bin_cutpoints', 'replicate_counts', 'draw_balance', 'background_gene_distances', 'catalog_exclusions')}
    replay = []
    canonical_catalog = [dict(g, local_annotated_gene_count=int(density_values[250000][i])) for i, g in enumerate(genes)]
    _, _, _, _, _, lbin, dbin = assignments(genes, targets[TARGETS[0]], next(s for s in model_specs() if s['model_id'] == 'L5_D3_R250000'), density_values)
    for i, row in enumerate(canonical_catalog):
        row['length_bin'] = int(lbin[i]) if lbin[i] >= 0 else None
        row['gene_density_bin'] = int(dbin[i]) if dbin[i] >= 0 else None
        for j, k in enumerate(keys):
            row['%s_%s_distance_bp' % k] = int(distances[i, j]) if np.isfinite(distances[i, j]) else None
    replay.append(_compare_rows(ctx, 'baseline_catalog_distances', canonical_catalog, frozen['background_gene_distances'], ('gene_id',), tuple(canonical_catalog[0])))
    replay.append(_compare_rows(ctx, 'baseline_catalog_exclusions', exclusions, frozen['catalog_exclusions'], ('gene_id',), ('reason',)))
    catalog_out = []
    for i, g in enumerate(genes):
        row = dict(g)
        for r in SETTINGS['gene_density_radius_bp']:
            row['other_gene_midpoints_within_%d_bp' % r] = int(density_values[r][i])
        for j, k in enumerate(keys):
            row['%s_%s_distance_bp' % k] = int(distances[i, j]) if np.isfinite(distances[i, j]) else None
        catalog_out.append(row)
    ctx.write_csv('gene_catalog_and_distances.csv', catalog_out)
    ctx.write_csv('catalog_exclusions.csv', exclusions, ['gene_id', 'reason'])
    specs = model_specs()
    ctx.write_csv('background_model_settings.csv', specs)
    status_rows, strata_rows, cut_rows, bin_rows = [], [], [], []
    summaries, observed_rows = [], []
    count_n, balance_n = 0, 0
    count_stream = None
    balance_stream = None
    canonical = {k: [] for k in ('status', 'summary', 'stratum_requirements', 'bin_cutpoints', 'replicate_counts', 'draw_balance')}
    lengths_array = np.asarray([g['gene_length_bp'] for g in genes], dtype=float)
    candidate_flags = np.asarray([g['excluded_candidate'] for g in genes])
    ids = np.asarray([g['gene_id'] for g in genes], dtype=object)
    replay_fields = {
        'status': (('target_set', 'model'), ('status', 'n_target_genes', 'deficient_strata', 'replicates')),
        'summary': (('target_set', 'model', 'resource', 'arm', 'flank_bp'), ('n_candidates', 'replicates', 'observed_fraction', 'background_mean_fraction', 'background_q025', 'background_q975', 'observed_minus_background_mean')),
        'stratum_requirements': (('target_set', 'model', 'chromosome', 'biotype', 'length_bin', 'gene_density_bin'), ('candidate_genes', 'available_noncandidate_genes', 'adequate')),
        'bin_cutpoints': (('target_set', 'biotype', 'feature', 'cut_number'), ('cut_value',)),
        'replicate_counts': (('target_set', 'model', 'resource', 'arm', 'replicate'), tuple(['n_genes'] + ['n_recovered_%d_bp' % w for w in windows])),
        'draw_balance': (('target_set', 'model', 'replicate'), ('n_genes', 'sample_gene_id_sha256', 'candidate_mean_length_bp', 'background_mean_length_bp', 'candidate_mean_local_gene_count', 'background_mean_local_gene_count')),
    }
    for canonical_phase in (True, False):
        for target, target_ix in targets.items():
            n = len(target_ix)
            observed = covered[target_ix].sum(axis=0, dtype=np.int64)
            for j, k in enumerate(keys):
                for w, count in zip(windows, observed[j]):
                    if canonical_phase:
                        observed_rows.append(dict(target_set=target, resource=k[0], arm=k[1], flank_bp=w,
                                              n_candidates=n, observed_count=int(count), observed_fraction=int(count) / n))
            for spec in (s for s in specs if s['canonical'] == canonical_phase):
                print('[BACKGROUND]', target, spec['model_id'], flush=True)
                prefix = dict(target_set=target, model_id=spec['model_id'], model=spec['model'],
                              length_bins_requested=spec['length_bins_requested'],
                              gene_density_bins_requested=spec['gene_density_bins_requested'],
                              gene_density_radius_bp=spec['gene_density_radius_bp'])
                req, pools, strata, cuts, diagnostics, _, _ = assignments(genes, target_ix, spec, density_values)
                for row in strata:
                    strata_rows.append(dict(prefix, **row))
                for row in cuts:
                    cut_rows.append(dict(prefix, **row))
                for row in diagnostics:
                    bin_rows.append(dict(prefix, **row))
                deficits = sum(not r['adequate'] for r in strata)
                status = dict(prefix, status='BLOCKED_INSUFFICIENT_STRATA' if deficits else 'COMPLETED_DESCRIPTIVE_BACKGROUND',
                              n_target_genes=n, deficient_strata=deficits, replicates=0 if deficits else SETTINGS['replicates'])
                status_rows.append(status)
                if spec['canonical']:
                    canonical['status'].append(dict(target_set=target, model=spec['model'], **{k: status[k] for k in ('status', 'n_target_genes', 'deficient_strata', 'replicates')}))
                    canonical['stratum_requirements'].extend(dict(target_set=target, model=spec['model'], **r) for r in strata)
                    if spec['model'] == DENSITY_MODEL:
                        canonical['bin_cutpoints'].extend(dict(target_set=target, **r) for r in cuts)
                if deficits:
                    continue
                rng = rng_for(SETTINGS['seed'], spec['rng_label'])
                arr = np.zeros((SETTINGS['replicates'], len(keys), len(windows)), dtype=np.int64)
                r_density = spec['gene_density_radius_bp'] or 250000
                draw_length, draw_density = np.zeros(SETTINGS['replicates']), np.zeros(SETTINGS['replicates'])
                length_observed = float(lengths_array[target_ix].mean())
                density_observed = float(density_values[r_density][target_ix].mean())
                min_draw_n, max_draw_n, excluded_draw_cells = n, 0, 0
                for rep in range(SETTINGS['replicates']):
                    ix = sample_strata(pools, req, rng)
                    if len(ix) != n or np.any(candidate_flags[ix]):
                        raise ValueError('Background draw violates candidate exclusion or fixed size')
                    min_draw_n, max_draw_n = min(min_draw_n, len(ix)), max(max_draw_n, len(ix))
                    excluded_draw_cells += int(candidate_flags[ix].sum())
                    arr[rep] = covered[ix].sum(axis=0, dtype=np.int64)
                    draw_length[rep] = float(lengths_array[ix].mean())
                    draw_density[rep] = float(density_values[r_density][ix].mean())
                    digest = hashlib.sha256('\n'.join(sorted(ids[ix])).encode('utf-8')).hexdigest()
                    balance = dict(prefix, replicate=rep + 1, n_genes=n, sample_gene_id_sha256=digest,
                                   balance_density_radius_bp=r_density,
                                   candidate_mean_length_bp=length_observed, background_mean_length_bp=float(draw_length[rep]),
                                   candidate_mean_local_gene_count=density_observed,
                                   background_mean_local_gene_count=float(draw_density[rep]))
                    if balance_stream is None:
                        balance_path = ctx.out / 'tables' / 'draw_balance.csv'
                        balance_file = stack.enter_context(balance_path.open('w', newline='', encoding='utf-8-sig'))
                        balance_stream = csv.DictWriter(balance_file, fieldnames=list(balance))
                        balance_stream.writeheader()
                    balance_stream.writerow(balance)
                    balance_n += 1
                    if spec['canonical']:
                        canonical['draw_balance'].append(balance)
                    for j, k in enumerate(keys):
                        count = dict(prefix, resource=k[0], arm=k[1], replicate=rep + 1, n_genes=n,
                                     **{'n_recovered_%d_bp' % w: int(arr[rep, j, h]) for h, w in enumerate(windows)})
                        if count_stream is None:
                            count_path = ctx.out / 'tables' / 'background_replicate_counts.csv'
                            count_file = stack.enter_context(count_path.open('w', newline='', encoding='utf-8-sig'))
                            count_stream = csv.DictWriter(count_file, fieldnames=list(count))
                            count_stream.writeheader()
                        count_stream.writerow(count)
                        count_n += 1
                        if spec['canonical']:
                            canonical['replicate_counts'].append(count)
                ctx.check('valid_draw_size_exclusion:' + target + ':' + spec['model_id'], [min_draw_n, max_draw_n, excluded_draw_cells], [n, n, 0])
                ctx.check('monotone_window_coverage:' + target + ':' + spec['model_id'], bool(np.all(np.diff(arr, axis=2) >= 0)), True)
                ctx.check('B_is_subset_of_A:' + target + ':' + spec['model_id'], bool(all(np.all(arr[:, 3*j+1, :] <= arr[:, 3*j, :]) for j in range(len(cfg['comparison_resources'])))), True)
                f = arr / n
                lo, hi = np.quantile(f, [0.025, 0.975], axis=0, method='linear')
                for j, k in enumerate(keys):
                    for h, w in enumerate(windows):
                        obs = int(observed[j, h])
                        percentile = empirical_percentile(arr[:, j, h], obs)
                        observed_fraction = obs / n
                        location = 'ABOVE_Q975' if observed_fraction > hi[j, h] + 1e-12 else ('BELOW_Q025' if observed_fraction < lo[j, h] - 1e-12 else 'WITHIN_Q025_Q975')
                        summary = dict(prefix, resource=k[0], arm=k[1], query_set=target, flank_bp=w,
                                       n_candidates=n, replicates=SETTINGS['replicates'], observed_count=obs,
                                       observed_fraction=observed_fraction, background_mean_fraction=float(f[:, j, h].mean()),
                                       background_q025=float(lo[j, h]), background_q975=float(hi[j, h]),
                                       observed_minus_background_mean=observed_fraction - float(f[:, j, h].mean()),
                                       location_relative_to_draw_interval=location, **percentile)
                        summaries.append(summary)
                        if spec['canonical']:
                            canonical['summary'].append(summary)
        if canonical_phase:
            for table, (key_fields, fields) in replay_fields.items():
                replay.append(_compare_rows(ctx, 'baseline_' + table, canonical[table], frozen[table], key_fields, fields))
            ctx.write_csv('baseline_replay_comparisons.csv', replay)
            if any(r['status'] != 'PASS' for r in replay):
                raise ValueError('Canonical background replay failed; sensitivity outputs must not be accepted')
    if count_stream is not None:
        count_file.close()
    if balance_stream is not None:
        balance_file.close()
    ctx.write_csv('baseline_replay_comparisons.csv', replay)
    ctx.write_csv('model_evaluation_status.csv', status_rows)
    ctx.write_csv('stratum_requirements.csv', strata_rows)
    ctx.write_csv('bin_cutpoints.csv', cut_rows)
    ctx.write_csv('bin_occupancy.csv', bin_rows)
    ctx.write_csv('candidate_observed_coverage.csv', observed_rows)
    ctx.write_csv('background_sensitivity_summary.csv', summaries)
    primary = [r for r in summaries if r['flank_bp'] == SETTINGS['primary_contrast_flank_bp'] and r['arm'] == 'C']
    ctx.write_csv('background_50kb_corrected_summary.csv', primary)
    ranges = []
    for target in targets:
        for resource in cfg['comparison_resources']:
            for model in (LENGTH_MODEL, DENSITY_MODEL):
                rr = [r for r in primary if r['target_set'] == target and r['resource'] == resource and r['model'] == model]
                if rr:
                    ranges.append(dict(target_set=target, resource=resource, model=model, arm='C', flank_bp=50000,
                                       evaluated_settings=len(rr), min_percentile_midrank=min(r['percentile_midrank'] for r in rr),
                                       max_percentile_midrank=max(r['percentile_midrank'] for r in rr),
                                       min_observed_minus_background_mean=min(r['observed_minus_background_mean'] for r in rr),
                                       max_observed_minus_background_mean=max(r['observed_minus_background_mean'] for r in rr),
                                       settings_above_q975=sum(r['location_relative_to_draw_interval'] == 'ABOVE_Q975' for r in rr),
                                       settings_within_interval=sum(r['location_relative_to_draw_interval'] == 'WITHIN_Q025_Q975' for r in rr),
                                       settings_below_q025=sum(r['location_relative_to_draw_interval'] == 'BELOW_Q025' for r in rr)))
    ctx.write_csv('background_50kb_model_ranges.csv', ranges)
    ctx.summary['background_sensitivity'] = dict(settings=SETTINGS, n_model_settings=len(specs), n_target_sets=len(targets),
        annotation_genes=len(genes), all_original_candidate_ids_excluded=len(excluded_ids),
        completed_target_models=sum(r['status'] == 'COMPLETED_DESCRIPTIVE_BACKGROUND' for r in status_rows),
        blocked_target_models=sum(r['status'] == 'BLOCKED_INSUFFICIENT_STRATA' for r in status_rows),
        total_draws=balance_n, replicate_count_rows=count_n, summary_rows=len(summaries),
        canonical_replay_status='PASS', canonical_replay_tables=len(replay),
        primary_50kb_corrected_canonical=[r for r in primary if r['model_id'] in ('L5', 'L5_D3_R250000')],
        primary_50kb_corrected_ranges=ranges)
    ctx.note('Density means other annotated-gene midpoints on the same chromosome within the specified inclusive radius; it is not marker density.')
    ctx.note('All original candidate IDs remain excluded from every pool, including coding-only analyses. Sparse strata block the full model without omitting candidates or relaxing matching.')
    ctx.note('Empirical percentiles are descriptive with explicit ties; the 2.5th-97.5th draw interval is not a confidence interval or an enrichment significance test.')
    ctx.note('The candidate list was assembled around published associations, including G7 and MAF5 discoveries; G7 and MAF5 contrasts retain discovery ascertainment. Diversity shares the sequence/marker resource and is not an independent effect-validation cohort. Matching annotation features does not recreate the tested-GWAS universe or remove spatial dependence.')
    ctx.note('The 3-bin, 250-kb original density choice is a fixed descriptive setting, not an empirically optimized or biologically privileged model. The full prespecified grid is reported, including blocked models.')


def run(ctx):
    with ExitStack() as stack:
        return _run(ctx, stack)
