"""Audit source labels and positional recovery.

Raw source labels remain authoritative. A separate provisional sensitivity flags
filter-like values in the method column; it never silently rewrites the source.
Exact positional representation is not proof of variant or allele identity.
"""
from bisect import bisect_left
from collections import Counter, defaultdict
import json
import math
import re

from xlsx_read import Book


RESOURCES = ('Diversity', 'Nacional_G7', 'Nacional_MAF5')
KNOWN_FILTERS = {'G7', 'MAF5', 'SSR'}
EVENT_FILE = 'source_data/02_association_units.xlsx'
COHORT_FILE = 'alignment_coverage/02_fixed_marker_cohorts.xlsx'
COVERAGE_FILE = 'alignment_coverage/04_paired_coverage_effects.xlsx'


def text(value):
    return '' if value is None else str(value).strip()


def integer(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(re.sub(r'[\s,\u00a0]', '', str(value)))
    except ValueError:
        return None
    return int(result) if math.isfinite(result) and result.is_integer() else None


def chromosome(value):
    result = integer(value)
    return str(result) if result is not None else text(value)


def flag(value):
    if isinstance(value, bool):
        return value
    if value in (0, 1):
        return bool(value)
    raise ValueError('Expected explicit Boolean, received ' + repr(value))


def event_tuple(row):
    return (chromosome(row['chromosome']), integer(row['peak_bp']),
            text(row['trait']), text(row['method_reported']),
            text(row['marker_filter_reported']))


def peak_key(row, position='peak_bp', chrom='chromosome'):
    return chromosome(row.get(chrom)), integer(row.get(position))


def sort_point(point):
    return int(point[0]), point[1]


def read_raw_master(path):
    """Read the named source sheet directly, preserving literal strings such as NA."""
    parsed, ignored, header = [], [], None
    required = {'Chromosome', 'Position of the association peak (bp)',
                'GWAS method', 'Sorting of marker', 'Traits'}
    with Book(path) as book:
        if 'synth_tot' not in book.sheets:
            raise ValueError('Original association workbook has no synth_tot sheet')
        for rownum, values in book.rows('synth_tot'):
            if header is None:
                clean = [text(v) for v in values]
                if required.issubset(set(clean)):
                    if len(clean) != len(set(clean)):
                        raise ValueError('Duplicate synth_tot header names')
                    header = clean
                continue
            if not any(value is not None and text(value) != '' for value in values):
                continue
            record = {name: values[i] if i < len(values) else None
                      for i, name in enumerate(header)}
            chrom = chromosome(record['Chromosome'])
            peak = integer(record['Position of the association peak (bp)'])
            trait = text(record['Traits'])
            if chrom not in {str(i) for i in range(1, 11)} or peak is None or peak <= 0 or not trait:
                ignored.append({'source_sheet': 'synth_tot', 'source_row': rownum,
                                'reason': 'INVALID_PRIMARY_CHROMOSOME_PEAK_OR_TRAIT',
                                'raw_record_json': json.dumps(record, ensure_ascii=False)})
                continue
            sorting = text(record['Sorting of marker'])
            method = text(record['GWAS method'])
            candidate = method if sorting not in KNOWN_FILTERS and method in KNOWN_FILTERS else ''
            parsed.append({'source_sheet': 'synth_tot', 'source_row': rownum,
                           'chromosome': chrom, 'peak_bp': peak, 'trait': trait,
                           'method_reported': method, 'marker_filter_reported': sorting,
                           'provisional_filter_candidate': candidate,
                           'provisional_basis': ('FILTER_LIKE_LABEL_IN_METHOD_FIELD_UNRESOLVED'
                                                 if candidate else ''),
                           'raw_record_json': json.dumps(record, ensure_ascii=False)})
    if header is None or not parsed:
        raise ValueError('No valid raw synth_tot source records parsed')
    return parsed, ignored


def nearest(points, target):
    if not points:
        return None
    j = bisect_left(points, target)
    candidates = []
    if j < len(points):
        candidates.append(abs(points[j] - target))
    if j:
        candidates.append(abs(points[j - 1] - target))
    return min(candidates)


def recovered(distance, flank):
    return distance is not None and distance <= flank


def run(ctx):
    raw, invalid = read_raw_master(ctx.source_paths['associations'])
    master = ctx.read(EVENT_FILE, 'master_events')
    old_raw = [r for r in ctx.read(EVENT_FILE, 'raw_records')
               if r['source_sheet'] == 'synth_tot' and flag(r['valid_peak'])]
    old_peaks = ctx.read(EVENT_FILE, 'unique_peaks')
    ledger = ctx.read(COHORT_FILE, 'marker_record_ledger')
    mapped = ctx.read(COHORT_FILE, 'all_panel_coordinates')
    inventory = ctx.read('source_data/04_marker_identity.xlsx', 'marker_inventory')
    frozen_distances = ctx.read(COVERAGE_FILE, 'query_distances')
    frozen_summary = ctx.read(COVERAGE_FILE, 'summary')

    raw_by_row = {int(r['source_row']): r for r in raw}
    old_by_row = {int(r['source_row']): r for r in old_raw}
    differences = []
    for rownum in sorted(set(raw_by_row) | set(old_by_row)):
        a, b = raw_by_row.get(rownum), old_by_row.get(rownum)
        if a is None or b is None or event_tuple(a) != event_tuple(b):
            differences.append({'check': 'source_row_event_fields_vs_frozen_raw_record', 'source_row': rownum,
                                'observed': json.dumps(event_tuple(a)) if a else 'MISSING',
                                'expected': json.dumps(event_tuple(b)) if b else 'MISSING'})
    raw_events = {event_tuple(r) for r in raw}
    old_events = {event_tuple(r) for r in master}
    for item in sorted(raw_events ^ old_events, key=str):
        differences.append({'check': 'distinct_event_vs_frozen_master', 'source_row': '',
                            'observed': json.dumps(item) if item in raw_events else 'MISSING',
                            'expected': json.dumps(item) if item in old_events else 'MISSING'})
    ctx.check('source_master_valid_rows', len(raw), len(old_raw))
    ctx.check('source_master_distinct_events', len(raw_events), len(master))
    ctx.check('source_master_event_field_mismatches', len(differences), 0)
    ctx.write_csv('01_source_master_records.csv', raw)
    ctx.write_csv('01_source_invalid_records.csv', invalid,
                  fieldnames=['source_sheet', 'source_row', 'reason', 'raw_record_json'])
    ctx.write_csv('01_source_frozen_comparison_differences.csv', differences,
                  fieldnames=['check', 'source_row', 'observed', 'expected'])
    anomalies = [r for r in raw if r['provisional_filter_candidate'] or
                 r['marker_filter_reported'] not in KNOWN_FILTERS or
                 r['method_reported'] in KNOWN_FILTERS]
    ctx.write_csv('01_source_filter_field_anomalies.csv', anomalies)

    event_rows = {}
    for row in raw:
        event_rows.setdefault(event_tuple(row), row)
    strict, provisional = defaultdict(set), defaultdict(set)
    event_counts, raw_counts = Counter(), Counter()
    by_peak = defaultdict(list)
    for row in raw:
        raw_counts[row['marker_filter_reported']] += 1
    for row in event_rows.values():
        point = peak_key(row)
        strict[row['marker_filter_reported']].add(point)
        proposed = row['provisional_filter_candidate'] or row['marker_filter_reported']
        provisional[proposed].add(point)
        event_counts[row['marker_filter_reported']] += 1
        by_peak[point].append(row)
    all_points = set(by_peak)
    old_point_set = {peak_key(row, 'start_bp') for row in old_peaks}
    ctx.check('source_unique_peak_positions', len(all_points), len(old_peaks))
    ctx.check('source_unique_peak_position_set_mismatches', len(all_points ^ old_point_set), 0)
    peak_ids = {peak_key(row, 'start_bp'): row['query_id'] for row in old_peaks}
    if len(peak_ids) != len(old_peaks):
        raise ValueError('Frozen peak table contains duplicate chromosome-position keys')
    filter_rows = []
    for label in sorted(set(strict) | set(provisional)):
        filter_rows.append({'raw_filter_label': label, 'raw_records': raw_counts[label],
                            'distinct_events': event_counts[label],
                            'strict_unique_peak_positions': len(strict[label]),
                            'provisional_method_field_sensitivity_peak_positions': len(provisional[label]),
                            'provisional_is_final_correction': False,
                            'note': 'Filter-specific peak sets overlap; do not sum as disjoint counts.'})
    ctx.write_csv('01_filter_denominators.csv', filter_rows)
    memberships, membership_counts = [], Counter()
    for point in sorted(all_points, key=sort_point):
        labels = sorted(label for label, peaks in strict.items() if point in peaks)
        proposed_labels = sorted(label for label, peaks in provisional.items() if point in peaks)
        signature = ';'.join(labels)
        membership_counts[signature] += 1
        memberships.append({'query_id': peak_ids.get(point, ''), 'chromosome': point[0],
                            'peak_bp': point[1], 'strict_filter_membership': signature,
                            'provisional_filter_membership': ';'.join(proposed_labels),
                            'in_strict_G7': point in strict['G7'],
                            'in_strict_MAF5': point in strict['MAF5'],
                            'in_strict_SSR': point in strict['SSR'],
                            'outside_strict_G7_MAF5_union': point not in strict['G7'] | strict['MAF5'],
                            'source_rows': ';'.join(str(r['source_row']) for r in by_peak[point]),
                            'n_distinct_events': len(by_peak[point])})
    ctx.write_csv('01_peak_source_memberships.csv', memberships)
    ctx.write_csv('01_peak_filter_overlap_groups.csv',
                  [{'strict_filter_membership': name, 'unique_peak_positions': n}
                   for name, n in sorted(membership_counts.items())])

    positions = defaultdict(set)
    record_sets = defaultdict(set)
    for row in ledger:
        resource = row['resource']
        if resource not in RESOURCES:
            continue
        for arm in 'ABC':
            if flag(row['in_' + arm]):
                c = chromosome(row['reference_chromosome'] if arm == 'C' else row['original_chromosome'])
                p = integer(row['reference_snp_bp'] if arm == 'C' else row['original_position_bp'])
                if c not in {str(i) for i in range(1, 11)} or p is None or p <= 0:
                    raise ValueError('Retained ledger coordinate is not a valid positive primary-chromosome point')
                positions[(resource, arm, c)].add(p)
                record_sets[(resource, arm)].add(row['record_key'])
    for resource in RESOURCES:
        ctx.check(resource + '_B_C_record_membership_mismatches',
                  len(record_sets[(resource, 'B')] ^ record_sets[(resource, 'C')]), 0)
        ctx.check(resource + '_B_records_outside_A',
                  len(record_sets[(resource, 'B')] - record_sets[(resource, 'A')]), 0)
    positions = {key: sorted(value) for key, value in positions.items()}
    frozen = {}
    for row in frozen_distances:
        if row['query_set'] == 'master_unique_peaks' and row['resource'] in RESOURCES:
            key = (row['resource'], chromosome(row['chromosome']), integer(row['start_bp']))
            if key in frozen:
                raise ValueError('Duplicate frozen master peak-distance row')
            frozen[key] = row
    coverage, new_distances, distance_errors = [], {}, []
    for resource in RESOURCES:
        for point in sorted(all_points, key=sort_point):
            key = (resource, point[0], point[1])
            old = frozen.get(key)
            result = {'resource': resource, 'query_id': peak_ids.get(point, ''),
                      'chromosome': point[0], 'peak_bp': point[1],
                      'strict_filter_membership': ';'.join(sorted(k for k, p in strict.items() if point in p)),
                      'scope': 'POSITIONAL_REPRESENTATION_NOT_VARIANT_IDENTITY_OR_PROXY_VALIDATION'}
            for arm in 'ABC':
                distance = nearest(positions.get((resource, arm, point[0]), []), point[1])
                new_distances[(resource, arm, point)] = distance
                result['distance_' + arm + '_bp'] = distance
                result['exact_' + arm] = recovered(distance, 0)
                result['within_50kb_' + arm] = recovered(distance, 50000)
                if old is None or distance != integer(old['distance_' + arm + '_bp']):
                    distance_errors.append({'resource': resource, 'query_id': peak_ids.get(point, ''),
                                            'arm': arm, 'observed_distance_bp': distance,
                                            'expected_distance_bp': (old.get('distance_' + arm + '_bp') if old else 'MISSING_ROW')})
            coverage.append(result)
    ctx.check('frozen_master_peak_distance_row_count', len(frozen), len(all_points) * len(RESOURCES))
    ctx.check('recomputed_master_peak_distance_mismatches', len(distance_errors), 0)
    ctx.write_csv('02_peak_positional_coverage.csv', coverage)
    ctx.write_csv('02_peak_distance_validation_differences.csv', distance_errors,
                  fieldnames=['resource', 'query_id', 'arm', 'observed_distance_bp', 'expected_distance_bp'])
    recovery = []
    frozen_counts = {(r['resource'], r['arm'], integer(r['flank_bp'])): r
                     for r in frozen_summary if r['query_set'] == 'master_unique_peaks'
                     and r['resource'] in RESOURCES and integer(r['flank_bp']) in (0, 50000)}
    for resource in RESOURCES:
        subsets = [('ALL_MASTER_PEAKS', 'ALL', all_points)]
        if resource.startswith('Nacional_'):
            label = resource.split('_', 1)[1]
            subsets.extend([('STRICT_RAW_FILTER_LABEL', label, strict[label]),
                            ('PROVISIONAL_METHOD_FIELD_SENSITIVITY', label, provisional[label])])
        for policy, label, points in subsets:
            for arm in 'ABC':
                for flank in (0, 50000):
                    n = sum(recovered(new_distances[(resource, arm, p)], flank) for p in points)
                    recovery.append({'resource': resource, 'denominator_policy': policy,
                                     'originating_filter': label, 'arm': arm, 'flank_bp': flank,
                                     'n_source_peak_positions': len(points), 'n_recovered': n,
                                     'n_not_recovered': len(points) - n,
                                     'fraction_recovered': n / len(points) if points else None,
                                     'automatic_100_percent_attainable_ceiling_assumed': False})
                    if policy == 'ALL_MASTER_PEAKS':
                        old = frozen_counts.get((resource, arm, flank))
                        ctx.check('coverage_count_' + resource + '_' + arm + '_' + str(flank),
                                  n, integer(old['n_recovered']) if old else 'MISSING')
                        ctx.check('coverage_denominator_' + resource + '_' + arm + '_' + str(flank),
                                  len(points), integer(old['n_queries']) if old else 'MISSING')
    ctx.write_csv('02_peak_recovery_summary.csv', recovery)

    # A Diversity record originally positioned at the published peak is a candidate
    # identity bridge only: mmc2 lacks focal marker IDs and source alleles.
    div_at, div_key = defaultdict(list), defaultdict(list)
    targets_by_key, ledger_by_row, mapping_by_row = defaultdict(list), {}, {}
    for row in inventory:
        if row['resource'] == 'Diversity':
            point = peak_key(row, 'position_bp')
            if point[1] is not None:
                div_at[point].append(row)
            if text(row.get('reference_key')):
                div_key[text(row['reference_key'])].append(row)
        elif row['resource'] in RESOURCES:
            targets_by_key[(row['resource'], text(row['canonical_id']))].append(row)
    for row in ledger:
        ledger_by_row[(row['resource'], int(row['original_row']))] = row
    for row in mapped:
        mapping_by_row[(row['resource'], int(row['original_row']))] = row
    misses, bridges = [], []
    for label in ('G7', 'MAF5'):
        resource = 'Nacional_' + label
        for point in sorted(strict[label], key=sort_point):
            if recovered(new_distances[(resource, 'C', point)], 0):
                continue
            divrows = div_at.get(point, [])
            keys = sorted({text(r.get('reference_key')) for r in divrows} - {''})
            target_rows = {int(r['source_row']): r for key in keys
                           for r in targets_by_key.get((resource, key), [])}
            local = [ledger_by_row[(resource, number)] for number in target_rows]
            has_accepted = any(flag(r['mapping_accepted']) for r in local)
            has_C = any(flag(r['in_C']) for r in local)
            small_linked_offset = any(flag(r['in_C']) and chromosome(r['reference_chromosome']) == point[0]
                                      and abs(integer(r['reference_snp_bp']) - point[1]) <= 1000 for r in local)
            reason = ('NO_DIVERSITY_ORIGINAL_POSITION_CANDIDATE' if not divrows else
                      'NO_CANDIDATE_REFERENCE_KEY' if not keys else
                      'NO_TARGET_MARKER_WITH_CANDIDATE_KEY' if not target_rows else
                      'CANDIDATE_TARGET_RECORDS_NOT_MAPPING_ACCEPTED' if not has_accepted else
                      'ACCEPTED_CANDIDATE_TARGET_RECORDS_NOT_IN_PAIRED_C' if not has_C else
                      'PAIRED_C_CANDIDATE_AT_OTHER_POSITION')
            misses.append({'resource': resource, 'query_id': peak_ids.get(point, ''),
                           'chromosome': point[0], 'peak_bp': point[1],
                           'denominator_policy': 'STRICT_RAW_FILTER_LABEL',
                           'diagnostic_class': reason,
                           'candidate_diversity_rows': ';'.join(str(r['source_row']) for r in divrows),
                           'candidate_reference_keys': ';'.join(keys),
                           'candidate_keys_nonunique_in_Diversity': ';'.join(k for k in keys if len(div_key[k]) != 1),
                           'target_source_rows_for_candidate_keys': ';'.join(str(k) for k in sorted(target_rows)),
                           'candidate_target_status_counts_json': json.dumps(dict(Counter(r['cohort_reason'] for r in local)), sort_keys=True),
                           'nearest_C_distance_bp': new_distances[(resource, 'C', point)],
                           'nearest_C_within_1000_bp': recovered(new_distances[(resource, 'C', point)], 1000),
                           'candidate_linked_C_within_1000_bp': small_linked_offset,
                           'focal_variant_identity_verified': False,
                           'note': 'Candidate bridge from original Diversity position and reference key; not a proven causal or source-GWAS SNP identity.'})
            for div in divrows:
                key = text(div.get('reference_key'))
                matched = targets_by_key.get((resource, key), []) if key else []
                for target in matched or [None]:
                    number = int(target['source_row']) if target else None
                    led = ledger_by_row.get((resource, number), {})
                    m = mapping_by_row.get((resource, number), {})
                    same_chr = chromosome(led.get('reference_chromosome')) == point[0]
                    p = integer(led.get('reference_snp_bp'))
                    bridges.append({'resource': resource, 'query_id': peak_ids.get(point, ''),
                                    'peak_chromosome': point[0], 'peak_bp': point[1],
                                    'diversity_source_row': div['source_row'],
                                    'diversity_marker_id': div['marker_id'], 'candidate_reference_key': key,
                                    'diversity_records_with_key': len(div_key.get(key, [])),
                                    'target_source_row': number,
                                    'target_marker_id': target['marker_id'] if target else '',
                                    'target_cohort_reason': led.get('cohort_reason', ''),
                                    'target_in_C': led.get('in_C'),
                                    'target_mapping_accepted': led.get('mapping_accepted'),
                                    'target_mapping_status': m.get('coordinate_status', ''),
                                    'reference_chromosome': led.get('reference_chromosome'),
                                    'reference_snp_bp': p,
                                    'same_chromosome_signed_offset_bp': p - point[1] if p is not None and same_chr else None,
                                    'scope': 'CANDIDATE_BRIDGE_ONLY'})
    ctx.write_csv('02_own_filter_exact_misses.csv', misses)
    ctx.write_csv('02_own_filter_miss_candidate_bridges.csv', bridges)
    ctx.summary['source_provenance'] = {
        'valid_master_source_rows': len(raw), 'distinct_master_events': len(raw_events),
        'master_unique_peaks': len(all_points), 'raw_filter_counts': filter_rows,
        'strict_G7_MAF5_union_peaks': len(strict['G7'] | strict['MAF5']),
        'outside_strict_G7_MAF5_union_peaks': len(all_points - (strict['G7'] | strict['MAF5'])),
        'source_field_anomaly_records': len(anomalies),
        'candidate_reassigned_records_for_sensitivity_only': sum(bool(r['provisional_filter_candidate']) for r in raw),
        'source_labels_overwritten': False,
    }
    ctx.summary['peak_recovery'] = {
        'positional_peak_resource_rows': len(coverage), 'summary': recovery,
        'strict_own_filter_exact_C_misses': dict(Counter(r['resource'] for r in misses)),
        'miss_candidate_diagnostic_counts': dict(Counter(r['resource'] + ':' + r['diagnostic_class'] for r in misses)),
        'proxy_LD_evaluated': False,
        'variant_identity_proven_from_coordinate_match': False,
    }
    ctx.note('Raw G7/MAF5/SSR/NA labels are preserved. Provisional method-field sensitivity is reported separately from the source-label analysis.')
    ctx.note('Exact recovery denotes an identical chromosome-position, not allele identity or independent validation. Missing-peak bridge classifications are candidate evidence because source peaks lack marker IDs and alleles.')
    ctx.note('This stage recomputes positional recovery from frozen marker ledgers; it does not rerun alignment or compute adjusted LD/proxy tagging.')
