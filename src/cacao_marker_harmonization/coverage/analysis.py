"""Coordinate-impact estimates and transparent, non-inferential reference contrasts."""
from __future__ import annotations
from collections import Counter, defaultdict
import hashlib
import numpy as np
from .common import DataError, text, integer
from .inputs import valid_interval
from .numerics import (Queries, positions, nearest, counts, paired_counts, rng_for,
                       quotas_by_chromosome, sample_points, gene_neighborhood_counts,
                       bin_cuts, bin_of, sample_strata)

ARM_DEFINITIONS = {
    'A': 'All original eligible marker records on supplied coordinates.',
    'B': 'Paired eligible alignment-retained marker records, original coordinates.',
    'C': 'Exactly the B records, reference-derived SNP coordinates.'
}


def position_sets(cohorts):
    return {r: {arm: positions(rr) for arm, rr in arms.items()} for r, arms in cohorts.items()}


def representation_analysis(cohorts, queries, windows, log):
    pointsets = position_sets(cohorts)
    summary, transition_summary, distances_rows, changes = [], [], [], []
    distances = {}
    core = {'unique_candidate_genes', 'master_unique_peaks', 'master_fruity_subset_peaks', 'master_source_intervals'}
    for resource, arms in pointsets.items():
        log.info('Paired coverage: %s | %s', resource, ','.join(arms))
        for name, q in queries.items():
            result = {}
            for arm, p in arms.items():
                d = nearest(q, p)
                result[arm] = d
                distances[resource, arm, name] = d
                for w, n in zip(windows, counts(d, windows)):
                    summary.append(dict(resource=resource, arm=arm, query_set=name,
                        flank_bp=w, n_queries=len(q), n_recovered=int(n), fraction_recovered=int(n)/len(q),
                        n_marker_records=len(cohorts[resource][arm]), n_unique_positions=sum(len(v) for v in p.values()),
                        interpretation='Fixed-query positional coverage; not genetic-effect portability.'))
            if 'C' in arms:
                for w in windows:
                    s = paired_counts(result['A'], result['B'], result['C'], w)
                    s.update(resource=resource, query_set=name, flank_bp=w,
                        selection_delta_fraction=s['selection_delta_count']/len(q),
                        coordinate_delta_fraction=s['coordinate_delta_count']/len(q),
                        total_delta_fraction=s['total_delta_count']/len(q))
                    transition_summary.append(s)
            if name in core:
                for i, row in enumerate(q.rows):
                    out = dict(resource=resource, query_set=name, query_id=row['query_id'],
                        chromosome=row['chromosome'], start_bp=row['start_bp'], end_bp=row['end_bp'])
                    for arm in arms:
                        d = result[arm][i]
                        out['distance_' + arm + '_bp'] = int(d) if np.isfinite(d) else None
                        out['same_chromosome_marker_in_' + arm] = bool(np.isfinite(d))
                    distances_rows.append(out)
            if name == 'unique_candidate_genes' and 'C' in arms:
                for i, row in enumerate(q.rows):
                    for w in windows:
                        hb, hc = bool(result['B'][i] <= w), bool(result['C'][i] <= w)
                        if hb != hc:
                            changes.append(dict(resource=resource, gene_id=row['query_id'], chromosome=row['chromosome'],
                                start_bp=row['start_bp'], end_bp=row['end_bp'], pathways=row.get('pathways',''), flank_bp=w,
                                B_recovered=hb, C_recovered=hc, transition='GAIN' if hc else 'LOSS',
                                B_distance_bp=int(result['B'][i]) if np.isfinite(result['B'][i]) else None,
                                C_distance_bp=int(result['C'][i]) if np.isfinite(result['C'][i]) else None))
    return dict(summary=summary, paired_differences=transition_summary,
                query_distances=distances_rows, changed_gene_windows=changes), pointsets, distances


def cohort_summaries(cohorts, lengths):
    rows = []
    for resource, arms in cohorts.items():
        for arm, rr in arms.items():
            p = positions(rr)
            for chrom in sorted(lengths):
                nrec = sum(x['chromosome'] == chrom for x in rr)
                rows.append(dict(resource=resource, arm=arm, chromosome=chrom,
                                 marker_records=nrec, unique_positions=len(p.get(chrom, [])),
                                 definition=ARM_DEFINITIONS[arm]))
    return rows


def count_matching(pointsets, queries, cfg, log):
    windows, nrep = cfg['windows_bp'], cfg['rarefaction_replicates']
    resources = cfg['comparison_resources']
    selected = cfg['rarefaction_query_sets']
    for name in selected:
        if name not in queries:
            raise DataError('Count-matching query set missing: ' + name)
    summaries, replicates, quota_rows = [], [], []
    for arm in ('A', 'C'):
        ps = {r: pointsets[r][arm] for r in resources}
        chromosomes = sorted(set().union(*(set(p) for p in ps.values())))
        quotas = quotas_by_chromosome(ps, chromosomes)
        for c in chromosomes:
            for r in resources:
                quota_rows.append(dict(arm=arm, chromosome=c, resource=r,
                    available_unique_positions=len(ps[r].get(c, [])), quota=quotas[c]))
        if sum(quotas.values()) == 0:
            raise DataError('No positions in the chromosome-matched quota.')
        for resource in resources:
            log.info('Count matching: %s/%s | %s replicates | %s positions',resource,arm,nrep,sum(quotas.values()))
            rng = rng_for(cfg['seed'], 'count-matching:' + arm + ':' + resource)
            arr = {name: np.zeros((nrep,len(windows)), dtype=np.int64) for name in selected}
            for rep in range(nrep):
                p = sample_points(ps[resource], quotas, rng)
                for name in selected:
                    arr[name][rep,:] = counts(nearest(queries[name], p), windows)
            for name, cc in arr.items():
                n = len(queries[name])
                raw = counts(nearest(queries[name], ps[resource]), windows)
                f = cc / n
                lo, hi = np.quantile(f, [0.025,0.975], axis=0, method='linear')
                for k, w in enumerate(windows):
                    summaries.append(dict(resource=resource, arm=arm, query_set=name, flank_bp=w,
                        n_queries=n, n_matched_unique_positions=sum(quotas.values()), replicates=nrep,
                        raw_fraction_recovered=int(raw[k])/n, matched_mean_fraction=float(f[:,k].mean()),
                        matched_q025=float(lo[k]), matched_q975=float(hi[k]),
                        note='Monte Carlo resampling spread, not a confidence interval. A and C quotas are computed separately; do not infer coordinate-only effects from this comparison.'))
                for rep in range(nrep):
                    replicates.append(dict(resource=resource, arm=arm, query_set=name, replicate=rep+1,
                        n_queries=n, **{f'n_recovered_{w}_bp':int(cc[rep,k]) for k,w in enumerate(windows)}))
    return dict(summary=summaries, chromosome_quotas=quota_rows, replicate_counts=replicates)


def background_catalog(gene_rows, candidate_ids, queries, lengths, bgcfg):
    """Fix an annotation universe; no tested-marker universe is inferred."""
    seen, genes, exclusions = set(), [], []
    for r in gene_rows:
        gid = text(r.get('gene_id'))
        if not gid or gid in seen:
            raise DataError('Empty/duplicate ID in genome-wide gene catalog.')
        seen.add(gid)
        if not valid_interval(r, lengths) or not text(r.get('biotype')):
            exclusions.append(dict(gene_id=gid,reason='NONPRIMARY_INVALID_INTERVAL_OR_MISSING_BIOTYPE'))
            continue
        genes.append(dict(query_id=gid, gene_id=gid, chromosome=text(r['chromosome']),
                          start_bp=integer(r['start_bp']),end_bp=integer(r['end_bp']),biotype=text(r['biotype']),
                          gene_length_bp=integer(r['end_bp'])-integer(r['start_bp'])+1,
                          excluded_candidate=gid in candidate_ids))
    genes.sort(key=lambda r:r['gene_id'])
    byid = {r['gene_id']:r for r in genes}
    target = queries[bgcfg['query_set']]
    for r in target.rows:
        g = byid.get(r['query_id'])
        if g is None or any(g[k] != r[k] for k in ('chromosome','start_bp','end_bp','biotype')):
            raise DataError('Target gene is not identical to the annotation catalog.')
    density = gene_neighborhood_counts(genes, bgcfg['local_gene_density_radius_bp'])
    pool_by_type = defaultdict(list)
    for i,g in enumerate(genes):
        g['local_annotated_gene_count'] = int(density[i])
        if not g['excluded_candidate']:
            pool_by_type[g['biotype']].append(i)
    types = sorted(set(r['biotype'] for r in target.rows))
    cuts, cut_rows = {}, []
    for bt in types:
        ix = pool_by_type.get(bt, [])
        if not ix:
            raise DataError('No noncandidate background genes for biotype: ' + bt)
        lc = bin_cuts([genes[i]['gene_length_bp'] for i in ix], bgcfg['length_bins'])
        dc = bin_cuts([genes[i]['local_annotated_gene_count'] for i in ix], bgcfg['local_gene_density_bins'])
        cuts[bt] = (lc,dc)
        for feature,values in [('gene_length_bp',lc),('local_annotated_gene_count',dc)]:
            for j,v in enumerate(values,1):
                cut_rows.append(dict(biotype=bt,feature=feature,cut_number=j,cut_value=float(v),
                                     construction='Quantiles of primary-chromosome noncandidate genes of the same biotype; ties share bins.'))
    for g in genes:
        if g['biotype'] in cuts:
            lc,dc = cuts[g['biotype']]
            g['length_bin'] = bin_of(g['gene_length_bp'],lc)
            g['gene_density_bin'] = bin_of(g['local_annotated_gene_count'],dc)
        else:
            g['length_bin'],g['gene_density_bin'] = None,None
    return genes, exclusions, cut_rows


def stratum(g, model):
    s = (g['chromosome'],g['biotype'],g['length_bin'])
    if model == 'chromosome_biotype_length_gene_density':
        s += (g['gene_density_bin'],)
    elif model != 'chromosome_biotype_length':
        raise DataError('Unknown annotation matching model: ' + model)
    return s


def background_diagnostics(gene_rows, all_candidate_ids, queries, lengths, pointsets, cfg, log):
    bg = cfg['background']
    if not bg['enabled']:
        return {'status':[dict(status='DISABLED_BY_CONFIGURATION')]}
    genes, exclusions, cut_rows = background_catalog(gene_rows, all_candidate_ids, queries, lengths, bg)
    q = Queries('annotation_gene_catalog', genes)
    target_ids = set(queries[bg['query_set']].ids)
    target_idx = np.array([i for i,g in enumerate(genes) if g['gene_id'] in target_ids],dtype=np.int64)
    if len(target_idx) != len(target_ids):
        raise DataError('Background target count differs from fixed candidate denominator.')
    windows, nrep, n = cfg['windows_bp'], bg['replicates'], len(target_ids)
    # Same matched gene set is used for all panels and arms in each replicate.
    keys = [(r,arm) for r in cfg['comparison_resources'] for arm in ('A','B','C')]
    all_d = {k:nearest(q,pointsets[k[0]][k[1]]) for k in keys}
    observed = {k:counts(all_d[k][target_idx],windows) for k in keys}
    status, strata_rows, summaries, reps, balance = [], [], [], [], []
    for model in bg['models']:
        required = Counter(stratum(genes[i],model) for i in target_idx)
        pools = defaultdict(list)
        for i,g in enumerate(genes):
            if not g['excluded_candidate'] and g['length_bin'] is not None:
                pools[stratum(g,model)].append(i)
        pools = {s:np.array(ix,dtype=np.int64) for s,ix in pools.items()}
        deficits = []
        for s,nreq in sorted(required.items()):
            avail = len(pools.get(s,[]))
            strata_rows.append(dict(model=model,chromosome=s[0],biotype=s[1],length_bin=s[2],
                                    gene_density_bin=s[3] if len(s)>3 else None,
                                    candidate_genes=nreq,available_noncandidate_genes=avail,
                                    adequate=avail>=nreq))
            if avail<nreq: deficits.append((s,nreq,avail))
        if deficits:
            status.append(dict(model=model,status='BLOCKED_INSUFFICIENT_STRATA',n_target_genes=n,
                               deficient_strata=len(deficits),replicates=0,
                               note='No widening of strata, sampling with replacement, or candidate omission was applied.'))
            continue
        status.append(dict(model=model,status='COMPLETED_DESCRIPTIVE_BACKGROUND',n_target_genes=n,
                           deficient_strata=0,replicates=nrep,note='No p-values, FDR or significance calls. This is not the tested-GWAS-marker universe.'))
        log.info('Annotated-gene reference contrast: %s | %s draws of %s genes',model,nrep,n)
        rng = rng_for(cfg['seed'],'annotation-background:'+model)
        arr = {k:np.zeros((nrep,len(windows)),dtype=np.int64) for k in keys}
        length_values=np.array([g['gene_length_bp'] for g in genes],dtype=float)
        density_values=np.array([g['local_annotated_gene_count'] for g in genes],dtype=float)
        for rep in range(nrep):
            ix = sample_strata(pools,required,rng)
            if len(ix)!=n or any(genes[i]['excluded_candidate'] for i in ix):
                raise DataError('Background draw violates candidate exclusion/size.')
            balance.append(dict(model=model,replicate=rep+1,n_genes=n,
                sample_gene_id_sha256=hashlib.sha256('\n'.join(sorted(genes[i]['gene_id'] for i in ix)).encode()).hexdigest(),
                candidate_mean_length_bp=float(length_values[target_idx].mean()),background_mean_length_bp=float(length_values[ix].mean()),
                candidate_mean_local_gene_count=float(density_values[target_idx].mean()),background_mean_local_gene_count=float(density_values[ix].mean())))
            for k in keys:
                arr[k][rep,:]=counts(all_d[k][ix],windows)
        for k,cc in arr.items():
            f=cc/n;lo,hi=np.quantile(f,[0.025,0.975],axis=0,method='linear')
            for j,w in enumerate(windows):
                summaries.append(dict(model=model,resource=k[0],arm=k[1],query_set=bg['query_set'],flank_bp=w,n_candidates=n,replicates=nrep,
                    observed_fraction=int(observed[k][j])/n,background_mean_fraction=float(f[:,j].mean()),
                    background_q025=float(lo[j]),background_q975=float(hi[j]),
                    observed_minus_background_mean=int(observed[k][j])/n-float(f[:,j].mean()),
                    note='Descriptive matched-annotation reference distribution; not an enrichment significance test, LD validation, or flavor-effect replication.'))
            for rep in range(nrep):
                reps.append(dict(model=model,resource=k[0],arm=k[1],replicate=rep+1,n_genes=n,
                                  **{f'n_recovered_{w}_bp':int(cc[rep,j]) for j,w in enumerate(windows)}))
    distance_table=[]
    for i,g in enumerate(genes):
        row=dict(g)
        for k in keys:
            d=all_d[k][i]
            row[f'{k[0]}_{k[1]}_distance_bp']=int(d) if np.isfinite(d) else None
        distance_table.append(row)
    return dict(status=status,summary=summaries,stratum_requirements=strata_rows,bin_cutpoints=cut_rows,
                replicate_counts=reps,draw_balance=balance,background_gene_distances=distance_table,catalog_exclusions=exclusions)
