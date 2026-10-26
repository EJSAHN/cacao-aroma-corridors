"""Audit existing coordinate evidence without changing alignment or source loci.

Exact positional matches do not establish SNP/allele identity or trait effects.
All alignment reasons below are the recorded operational filter decisions.
"""
from collections import Counter, defaultdict
import re


NACIONAL = ('Nacional_G7', 'Nacional_MAF5')


def text(value):
    return '' if value is None else str(value).strip()


def integer(value):
    if value is None or value == '':
        return None
    if isinstance(value, bool):
        raise ValueError('Boolean is not a coordinate')
    number = int(value)
    if float(value) != number:
        raise ValueError('Nonintegral source coordinate: ' + str(value))
    return number


def flag(value):
    if value is True or value == 1:
        return True
    if value is False or value == 0:
        return False
    raise ValueError('Expected explicit Boolean, got ' + repr(value))


def point(chromosome, position):
    return text(chromosome), integer(position)


def unique(rows, key, label):
    result = {}
    for row in rows:
        k = key(row)
        if k in result:
            raise ValueError('Duplicate ' + label + ': ' + repr(k))
        result[k] = row
    return result


def grouped(rows, field):
    result = defaultdict(list)
    for row in rows:
        k = text(row.get(field))
        if k:
            result[k].append(row)
    return result


def convention(chrom, pos, refchrom, start):
    same = text(chrom) == text(refchrom)
    pos, start = integer(pos), integer(start)
    if not same:
        return 'CHROMOSOME_DIFFERS'
    if pos is not None and start is not None and pos == start:
        return 'EQUALS_1BASED_TAG_START'
    if pos is not None and start is not None and pos == start - 1:
        return 'EQUALS_0BASED_TAG_START'
    return 'OTHER_SAME_CHROMOSOME'


def recorded_convention(row):
    if row.get('same_chromosome') is not True:
        return 'CHROMOSOME_DIFFERS'
    if row.get('source_equals_tag_start_1based') is True:
        return 'EQUALS_1BASED_TAG_START'
    if row.get('source_equals_tag_start_0based') is True:
        return 'EQUALS_0BASED_TAG_START'
    return 'OTHER_SAME_CHROMOSOME'


def run(ctx):
    inventory = ctx.read('source_data/04_marker_identity.xlsx', 'marker_inventory')
    ledger = ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx', 'marker_record_ledger')
    coordinates = ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx', 'all_panel_coordinates')
    decisions = ctx.read('alignment_coverage/01_realignment_evidence.xlsx', 'marker_decisions')
    allele_decisions = ctx.read('alignment_coverage/01_realignment_evidence.xlsx', 'allele_decisions')
    peaks = ctx.read('source_data/02_association_units.xlsx', 'unique_peaks')
    coverage = ctx.read('alignment_coverage/04_paired_coverage_effects.xlsx', 'summary')
    exact = ctx.read('reference_mapping/01_exact_tag_mapping.xlsx', 'marker_mapping')
    stored_conventions = ctx.read('reference_mapping/03_coordinate_conventions.xlsx', 'offset_strand_summary')
    sequences = ctx.read('source_evidence/01_sequence_identity.xlsx', 'sequence_checks')
    old = unique(inventory, lambda r: (r['resource'], integer(r['source_row'])), 'inventory record')
    coord = unique(coordinates, lambda r: (r['resource'], integer(r['original_row'])), 'coordinate record')
    cohort = unique(ledger, lambda r: (r['resource'], integer(r['original_row'])), 'cohort record')
    mapped = unique(decisions, lambda r: r['query_id'], 'alignment query')
    decision_by_source = unique(decisions, lambda r: integer(r['source_row']), 'Diversity alignment source row')
    alleles = unique(allele_decisions, lambda r: (r['query_id'], r['allele']), 'alignment allele')
    seq = unique(sequences, lambda r: integer(r['source_row']), 'sequence check source row')
    bykey = grouped(decisions, 'reference_key')
    exact_bykey = grouped(exact, 'reference_key')
    panel_key_counts = Counter((r['resource'], text(r['canonical_id'])) for r in inventory)
    ctx.check('coordinate_inventory_record_keys', set(coord) == set(old), True)
    ctx.check('coordinate_decision_source_rows', set(decision_by_source) ==
              {row for resource, row in old if resource == 'Diversity'}, True)

    # R1 C1: identify query positions lost, then show every contributing marker.
    div_a = defaultdict(list)
    div_b = set()
    for r in ledger:
        if r['resource'] != 'Diversity':
            continue
        p = point(r['original_chromosome'], r['original_position_bp'])
        if flag(r['in_A']):
            div_a[p].append(r)
        if flag(r['in_B']):
            div_b.add(p)
    peak_keys = unique(peaks, lambda r: point(r['chromosome'], r['start_bp']), 'master peak position')
    lost = sorted((set(peak_keys) & set(div_a)) - div_b, key=lambda p: (int(p[0]), p[1]))
    lost_rows, lost_alleles = [], []
    for p in lost:
        for r in div_a[p]:
            source_row = integer(r['original_row'])
            c = coord[('Diversity', source_row)]
            d = decision_by_source[source_row]
            o = old[('Diversity', source_row)]
            if c.get('mapped_query_id') != d['query_id']:
                raise ValueError('Lost-peak Diversity row/query identity drift')
            if flag(r['in_B']):
                raise ValueError('Lost query has a retained marker at the same original position')
            result = dict(query_id=peak_keys[p]['query_id'], chromosome=p[0], peak_bp=p[1],
                          marker_record_key=r['record_key'], diversity_source_row=source_row,
                          diversity_marker_id=r['marker_id'], source_reference_key=d['reference_key'],
                          mapping_query_id=d['query_id'], cohort_reason=r['cohort_reason'],
                          placement_status=d['placement_status'], ref_allele_status=d['ref_allele_status'],
                          alt_allele_status=d['alt_allele_status'], mapping_accepted=flag(r['mapping_accepted']),
                          mapped_chromosome=d.get('mapped_chromosome'), mapped_snp_bp=d.get('mapped_snp_bp'),
                          mapped_strand=d.get('mapped_strand'), ref_cigar=d.get('ref_cigar'),
                          alt_cigar=d.get('alt_cigar'), source_alleles=o.get('alleles_reported'),
                          source_assembly=o.get('assembly_reported'), source_remark=o.get('remark'),
                          interpretation='Recorded alignment-filter exclusion; coordinate colocation does not establish association SNP identity.')
            lost_rows.append(result)
            for allele in ('R', 'A'):
                a = alleles.get((d['query_id'], allele))
                if a is not None:
                    lost_alleles.append(dict(peak_query_id=peak_keys[p]['query_id'],
                                             peak_chromosome=p[0], peak_bp=p[1],
                                             diversity_source_row=source_row,
                                             diversity_marker_id=r['marker_id'], **a))
    declared = {(r['resource'], r['arm']): int(r['n_recovered']) for r in coverage
                if r['query_set'] == 'master_unique_peaks' and integer(r['flank_bp']) == 0}
    ctx.check('coordinate_diversity_A_exact_positions', len(set(peak_keys) & set(div_a)), declared[('Diversity', 'A')])
    ctx.check('coordinate_diversity_B_exact_positions', len(set(peak_keys) & div_b), declared[('Diversity', 'B')])
    ctx.check('coordinate_diversity_lost_A_to_B_peak_count', len(lost), declared[('Diversity', 'A')] - declared[('Diversity', 'B')])
    ctx.write_csv('coordinate_diversity_lost_peaks_markers.csv', lost_rows)
    ctx.write_csv('coordinate_diversity_lost_peaks_alleles.csv', lost_alleles)
    reason_summary = []
    for status in sorted({r['placement_status'] for r in lost_rows}):
        rr = [r for r in lost_rows if r['placement_status'] == status]
        reason_summary.append(dict(placement_status=status, marker_records=len(rr),
                                   distinct_peak_positions=len({r['query_id'] for r in rr})))
    ctx.write_csv('coordinate_diversity_lost_reason_summary.csv', reason_summary)

    # R1 C2: reconstruct the declared unique-reference-key bridge from all rows.
    link_rows = []
    for o in inventory:
        resource = o['resource']
        if resource not in NACIONAL:
            continue
        key = text(o['canonical_id'])
        rownum = integer(o['source_row'])
        c, l = coord[(resource, rownum)], cohort[(resource, rownum)]
        candidates = bykey.get(key, [])
        accepted = flag(c['accepted_for_coordinate_comparison'])
        paired = flag(l['in_C'])
        if not candidates:
            category = 'NO_REFERENCE_KEY_LINK'
        elif len(candidates) != 1:
            category = 'AMBIGUOUS_REFERENCE_KEY_LINK'
        elif accepted:
            category = 'ACCEPTED_PAIRED' if paired else 'ACCEPTED_UNPAIRED_ORIGINAL_COORDINATE'
        elif panel_key_counts[(resource, key)] != 1:
            category = 'UNIQUE_LINK_PANEL_KEY_NOT_UNIQUE'
        else:
            category = 'UNIQUE_LINK_ALIGNMENT_NOT_ACCEPTED'
        row = dict(resource=resource, source_row=rownum, marker_id=o['marker_id'],
                   canonical_key=key, reference_key_candidate_records=len(candidates),
                   panel_key_records=panel_key_counts[(resource, key)], linkage_category=category,
                   coordinate_status=c['coordinate_status'], in_A=flag(l['in_A']),
                   in_B=flag(l['in_B']), in_C=paired, accepted_coordinate=accepted,
                   original_chromosome=o.get('chromosome'), original_position_bp=o.get('position_bp'),
                   corrected_chromosome=c.get('reference_chromosome'), corrected_snp_bp=c.get('reference_snp_bp'),
                   comparison='NOT_EVALUABLE_WITHOUT_ACCEPTED_MAPPING_AND_UNIQUE_LINK',
                   sequence_identity_basis=c.get('identity_basis'))
        if len(candidates) == 1:
            d = candidates[0]
            source_row = integer(d['source_row'])
            od = old[('Diversity', source_row)]
            row.update(diversity_source_row=source_row, diversity_marker_id=od['marker_id'],
                       diversity_reference_key=d['reference_key'], mapping_query_id=d['query_id'],
                       diversity_original_chromosome=od.get('chromosome'),
                       diversity_original_snp_bp=od.get('position_bp'),
                       mapping_placement_status=d['placement_status'],
                       mapped_strand=d.get('mapped_strand'), ref_cigar=d.get('ref_cigar'),
                       alt_cigar=d.get('alt_cigar'), diversity_source_assembly=od.get('assembly_reported'),
                       diversity_source_remark=od.get('remark'),
                       sequence_supported_reference_key=seq[source_row].get('sequence_supported_reference_key'))
            if text(c.get('mapped_query_id')) != text(d['query_id']):
                raise ValueError('Unique-reference-key reconstruction disagrees with saved mapping query')
            if accepted:
                op = point(od.get('chromosome'), od.get('position_bp'))
                cp = point(c.get('reference_chromosome'), c.get('reference_snp_bp'))
                if op[0] == '' or op[1] is None:
                    comparison = 'LINKED_DIVERSITY_ORIGINAL_COORDINATE_MISSING'
                elif op == cp:
                    comparison = 'SAME_AS_LINKED_DIVERSITY_ORIGINAL'
                elif op[0] != cp[0]:
                    comparison = 'CHROMOSOME_CHANGED_FROM_LINKED_DIVERSITY_ORIGINAL'
                else:
                    comparison = 'POSITION_CHANGED_FROM_LINKED_DIVERSITY_ORIGINAL'
                row['comparison'] = comparison
                row['corrected_minus_diversity_original_bp'] = cp[1] - op[1] if op[0] == cp[0] and op[1] is not None else None
        elif accepted:
            raise ValueError('Accepted Nacional mapping lacks a unique reference-key link')
        link_rows.append(row)
    ctx.write_csv('coordinate_identifier_link_all_nacional.csv', link_rows)
    ctx.write_csv('coordinate_identifier_link_changed_accepted.csv', [r for r in link_rows
                  if r['accepted_coordinate'] and r['comparison'] != 'SAME_AS_LINKED_DIVERSITY_ORIGINAL'])
    link_summary, comparison_summary = [], []
    linkage_categories = ('NO_REFERENCE_KEY_LINK', 'AMBIGUOUS_REFERENCE_KEY_LINK',
                          'UNIQUE_LINK_PANEL_KEY_NOT_UNIQUE', 'UNIQUE_LINK_ALIGNMENT_NOT_ACCEPTED',
                          'ACCEPTED_UNPAIRED_ORIGINAL_COORDINATE', 'ACCEPTED_PAIRED')
    comparison_categories = ('SAME_AS_LINKED_DIVERSITY_ORIGINAL',
                             'POSITION_CHANGED_FROM_LINKED_DIVERSITY_ORIGINAL',
                             'CHROMOSOME_CHANGED_FROM_LINKED_DIVERSITY_ORIGINAL',
                             'LINKED_DIVERSITY_ORIGINAL_COORDINATE_MISSING')
    for resource in NACIONAL:
        rr = [r for r in link_rows if r['resource'] == resource]
        counts = Counter(r['linkage_category'] for r in rr)
        link_summary.extend(dict(resource=resource, linkage_category=k, marker_records=counts[k]) for k in linkage_categories)
        ctx.check('coordinate_link_partition_' + resource, sum(counts.values()), len(rr))
        ctx.check('coordinate_paired_link_count_' + resource, sum(r['in_C'] for r in rr),
                  sum(flag(r['in_C']) for r in ledger if r['resource'] == resource))
        for scope, subset in [('ALL_ACCEPTED_COORDINATES', [r for r in rr if r['accepted_coordinate']]),
                              ('PAIRED_C', [r for r in rr if r['in_C']])]:
            counts = Counter(r['comparison'] for r in subset)
            for name in comparison_categories:
                n = counts[name]
                comparison_summary.append(dict(resource=resource, scope=scope, comparison=name,
                                               records=n, denominator_records=len(subset),
                                               fraction=n / len(subset) if subset else None))
    ctx.write_csv('coordinate_identifier_link_summary.csv', link_summary)
    ctx.write_csv('coordinate_join_vs_alignment_summary.csv', comparison_summary)

    # R1 C3: diagnostic exact subset reconstructed independently of its aggregate.
    convention_rows = []

    def add_convention(scope, resource, o, d, source_row, refchrom, snp, tagstart, tagend, strand):
        m = re.search(r'--(\d+)$', text(o['canonical_id']))
        suffix = int(m.group(1)) if m else None
        sr = seq[integer(d['source_row'])]
        p, snp = integer(o.get('position_bp')), integer(snp)
        same = text(o.get('chromosome')) == text(refchrom)
        convention_rows.append(dict(scope=scope, resource=resource, original_row=source_row,
                                    marker_id=o['marker_id'], original_chromosome=o.get('chromosome'),
                                    original_position_bp=p, reference_chromosome=refchrom,
                                    reference_snp_bp=snp, reference_tag_start_bp=tagstart,
                                    reference_tag_end_bp=tagend, mapped_strand=strand,
                                    source_strand_reported=o.get('strand_reported'),
                                    source_minus_reference_snp_bp=p-snp if same and p is not None else None,
                                    convention_class=convention(o.get('chromosome'), p, refchrom, tagstart),
                                    original_position_missing=p is None,
                                    mapping_query_id=d['query_id'], source_reference_key=d['reference_key'],
                                    marker_terminal_index=suffix, linked_five_flank_length=sr.get('five_flank_length'),
                                    linked_three_flank_length=sr.get('three_flank_length'),
                                    terminal_equals_five_length=(suffix == integer(sr.get('five_flank_length'))) if suffix is not None else None,
                                    terminal_equals_three_length=(suffix == integer(sr.get('three_flank_length'))) if suffix is not None else None))

    for o in inventory:
        resource = o['resource']
        if resource not in NACIONAL:
            continue
        key, source_row = text(o['canonical_id']), integer(o['source_row'])
        candidates = exact_bykey.get(key, [])
        if len(candidates) == 1:
            d = candidates[0]
            if d['mapping_status'] == 'UNIQUE_EXACT_PLACEMENT' and flag(d['primary_chromosome']) and panel_key_counts[(resource, key)] == 1:
                add_convention('EXACT_ACCEPTED_DIAGNOSTIC', resource, o, d, source_row,
                               d['mapped_chromosome'], d['mapped_snp_bp'], d['tag_start_bp'], d['tag_end_bp'], d['mapped_strand'])
        c = coord[(resource, source_row)]
        if flag(c['accepted_for_coordinate_comparison']):
            d = mapped[c['mapped_query_id']]
            add_convention('GAPPED_ALL_ACCEPTED', resource, o, d, source_row,
                           c['reference_chromosome'], c['reference_snp_bp'], c['reference_tag_start_bp'],
                           c['reference_tag_end_bp'], c['source_tag_strand_on_reference'])
            if flag(cohort[(resource, source_row)]['in_C']):
                add_convention('GAPPED_PAIRED_C', resource, o, d, source_row,
                               c['reference_chromosome'], c['reference_snp_bp'], c['reference_tag_start_bp'],
                               c['reference_tag_end_bp'], c['source_tag_strand_on_reference'])
    observed = Counter((r['resource'], r['mapped_strand'], r['convention_class']) for r in convention_rows
                       if r['scope'] == 'EXACT_ACCEPTED_DIAGNOSTIC')
    recorded = Counter()
    for r in stored_conventions:
        if r['resource'] in NACIONAL:
            recorded[(r['resource'], r['mapped_strand'], recorded_convention(r))] += integer(r['records'])
    ctx.check('coordinate_exact_convention_cross_tab_reconstruction', dict(observed) == dict(recorded), True)
    cross = Counter((r['scope'], r['resource'], r['mapped_strand'], r['convention_class']) for r in convention_rows)
    cross_rows = [dict(scope=s, resource=res, mapped_strand=strand, convention_class=cc, records=n)
                  for (s,res,strand,cc),n in sorted(cross.items(), key=lambda item: str(item[0]))]
    suffix_counts = Counter((r['scope'], r['resource'], r['mapped_strand'], r['convention_class'],
                             r['terminal_equals_five_length'], r['terminal_equals_three_length'],
                             text(r['source_strand_reported'])) for r in convention_rows)
    suffix_rows = [dict(scope=s,resource=res,mapped_strand=st,convention_class=c,
                        terminal_equals_five_length=f,terminal_equals_three_length=t,
                        source_strand_reported=ss,records=n)
                   for (s,res,st,c,f,t,ss),n in sorted(suffix_counts.items(),key=lambda item:str(item[0]))]
    ctx.write_csv('coordinate_convention_record_details.csv', convention_rows)
    ctx.write_csv('coordinate_convention_by_mapped_strand.csv', cross_rows)
    ctx.write_csv('coordinate_convention_suffix_source_strand.csv', suffix_rows)
    deterministic = []
    for resource in NACIONAL:
        rr = [r for r in convention_rows if r['resource'] == resource and r['scope'] == 'EXACT_ACCEPTED_DIAGNOSTIC']
        classes = defaultdict(set)
        for r in rr:
            if r['convention_class'].startswith('EQUALS_'):
                classes[r['mapped_strand']].add(r['convention_class'])
        deterministic.append(dict(resource=resource,
                                  both_start_equalities_observed_with_same_mapped_strand=any(len(v)>1 for v in classes.values()),
                                  diagnostic_records=len(rr)))
    ctx.summary['coordinate_audit'] = dict(
        lost_diversity_peak_positions=len(lost), lost_diversity_marker_records=len(lost_rows),
        lost_diversity_filter_reasons=reason_summary, identifier_link_categories=link_summary,
        join_vs_alignment=comparison_summary, exact_strand_diagnostic=deterministic,
        alignment_rerun=False, automatic_coordinate_rule_applied_to_unlinked_markers=False,
        association_variant_identity_established=False)
    ctx.note('Coordinate audit reconstructs the existing unique-reference-key linkage and fixed paired cohorts; it does not remap tags or change source peak positions.')
    ctx.note('A Diversity marker at a published peak position is coordinate evidence, not proof that it was the discovery SNP. Lost-peak reasons are recorded algorithmic filters, not diagnosed biological causes.')
    ctx.note('Tag-start equality classes and mapped strand are observations on retained markers; no deterministic rule is extrapolated to unlinked markers.')
