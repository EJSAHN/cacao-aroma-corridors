"""Tag length, unique-position quota and marker composition checks.

Reads fixed evidence tables for tag lengths, quotas and marker composition.
"""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import math
import re

PANELS = ('Diversity', 'Nacional_G7', 'Nacional_MAF5')
NACIONAL = PANELS[1:]


def integer(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        v = float(value)
        return int(v) if math.isfinite(v) and v.is_integer() else None
    except (ValueError, TypeError, OverflowError):
        return None


def truth(value):
    if isinstance(value, bool):
        return value
    if value in (0, 1):
        return bool(value)
    if isinstance(value, str) and value.lower() in ('true', 'false'):
        return value.lower() == 'true'
    raise ValueError('Expected explicit Boolean: ' + repr(value))


def chrom(value):
    n = integer(value)
    return str(n) if n is not None and 1 <= n <= 10 else None


def quantile(values, p):
    if not values:
        return None
    pos = (len(values)-1)*p
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo] + (values[hi]-values[lo])*(pos-lo)


def tag_audit(ctx):
    seq = ctx.read('source_evidence/01_sequence_identity.xlsx', 'sequence_checks')
    inputs = ctx.read('alignment_coverage/01_realignment_evidence.xlsx', 'input_tag_audit')
    decisions = ctx.read('alignment_coverage/01_realignment_evidence.xlsx', 'marker_decisions')
    minimum = integer(ctx.baseline_config['alignment']['minimum_tag_length'])
    if minimum is None or minimum <= 0:
        raise ValueError('Invalid configured minimum tag length')
    eligible = {r['query_id'] for r in inputs if r['input_status']=='ELIGIBLE'}
    retained = {r['query_id'] for r in decisions if r['placement_status']=='RETAINED_PRIMARY_CHROMOSOME'}
    ctx.check('methods_unique_sequence_source_rows', len({r['source_row'] for r in seq}), len(seq))
    ctx.check('methods_unique_alignment_input_ids', len({r['query_id'] for r in inputs}), len(inputs))
    source_ids = {'diversity_row_'+str(integer(r['source_row'])) for r in seq}
    ctx.check('methods_alignment_input_ids_match_sequence_rows', {r['query_id'] for r in inputs}==source_ids, True)
    populations = defaultdict(list)
    exceptions, mismatches = [], []
    for r in seq:
        q = 'diversity_row_'+str(integer(r['source_row']))
        five, three = str(r.get('five_flank') or ''), str(r.get('three_flank') or '')
        span = len(five)+1+len(three)
        length = integer(r.get('tag_length'))
        valid = (re.fullmatch('[ACGT]*', five) is not None and re.fullmatch('[ACGT]*', three) is not None
                 and re.fullmatch('[ACGT]>[ACGT]', str(r.get('variation',''))) is not None)
        if (valid and length != span) or (not valid and length is not None):
            mismatches.append(q)
        populations['all_source_literal_flank_spans'].append(span)
        if length is not None:
            populations['source_reconstructable_tags'].append(length)
        for name, subset in [('alignment_input_eligible_tags',eligible), ('retained_primary_chromosome_tags',retained)]:
            if q in subset:
                if length is None:
                    raise ValueError('Eligible/retained tag lacks a valid length: '+q)
                populations[name].append(length)
        if length is None or q not in eligible:
            exceptions.append(dict(source_row=r['source_row'], query_id=q, marker_id=r['marker_id'],
                                   reference_key=r['reference_key'], variation=r['variation'],
                                   literal_flank_span_bp=span, supported_tag_length_bp=length,
                                   reference_lookup_status=r['ref_lookup_status'], alignment_input_eligible=q in eligible))
    ctx.check('methods_reconstructed_tag_length_mismatches', len(mismatches), 0)
    ctx.check('methods_eligible_tag_length_count', len(populations['alignment_input_eligible_tags']), len(eligible))
    ctx.check('methods_retained_tag_length_count', len(populations['retained_primary_chromosome_tags']), len(retained))
    summary, distribution = [], []
    for name, lengths in populations.items():
        lengths.sort()
        counts = Counter(lengths)
        summary.append(dict(population=name, source_marker_records=len(seq), records=len(lengths),
                            min_bp=lengths[0] if lengths else None, q25_bp=quantile(lengths,.25),
                            median_bp=quantile(lengths,.5), q75_bp=quantile(lengths,.75),
                            max_bp=lengths[-1] if lengths else None,
                            mean_bp=sum(lengths)/len(lengths) if lengths else None,
                            configured_minimum_tag_length=minimum, records_below_minimum=sum(x<minimum for x in lengths),
                            records_at_69bp=counts[69], quantile_definition='Linear interpolation at (n-1)*p'))
        distribution.extend(dict(population=name, tag_length_bp=k, records=v, fraction=v/len(lengths)) for k,v in sorted(counts.items()))
    ctx.write_csv('methods_tag_length_summary.csv', summary)
    ctx.write_csv('methods_tag_length_distribution.csv', distribution)
    ctx.write_csv('methods_tag_exceptions.csv', exceptions)
    ctx.note('Tag lengths distinguish literal source spans (which may contain unsupported bases), reconstructable tags, alignment-eligible tags and retained primary placements. The constant score floor is not a length-normalized error rate; this audit does not change the acceptance rule.')
    return summary


def ledger_index(ctx):
    ledger = ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx','marker_record_ledger')
    ctx.check('methods_unique_marker_record_keys',len({r['record_key'] for r in ledger}),len(ledger))
    positions, records, lookup = defaultdict(set), defaultdict(list), {}
    for r in ledger:
        resource = r['resource']
        if resource not in PANELS:
            continue
        lookup[(resource,integer(r['original_row']))] = r
        for arm in ('A','B','C'):
            if not truth(r['in_'+arm]):
                continue
            cf,pf = ('reference_chromosome','reference_snp_bp') if arm=='C' else ('original_chromosome','original_position_bp')
            c,p = chrom(r[cf]),integer(r[pf])
            if c is None or p is None or p<1:
                raise ValueError('Included marker has invalid coordinate: '+str(r['record_key'])+' '+arm)
            key = (resource,arm,c)
            positions[key].add(p)
            records[key].append((p,r))
    for key in list(positions):
        positions[key] = sorted(positions[key])
        records[key].sort(key=lambda x:(x[0],x[1]['record_key']))
    return positions,records,lookup


def quota_audit(ctx, positions, records):
    counts = ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx','chromosome_counts')
    quotas = ctx.read('alignment_coverage/05_count_matched_sensitivity.xlsx','chromosome_quotas')
    cs = {(r['resource'],r['arm'],chrom(r['chromosome'])):r for r in counts if r['resource'] in PANELS}
    qs = {(r['resource'],r['arm'],chrom(r['chromosome'])):r for r in quotas}
    expected = {(r,a,str(c)) for r in PANELS for a in ('A','C') for c in range(1,11)}
    ctx.check('methods_frozen_quota_keys',set(qs)==expected,True)
    detail, sums = [], Counter()
    for arm in ('A','B','C'):
        for c in map(str,range(1,11)):
            quota = min(len(positions.get((r,arm,c),[])) for r in PANELS) if arm in ('A','C') else None
            if quota is not None:
                sums[arm] += quota
            for resource in PANELS:
                key = (resource,arm,c)
                unique, total = len(positions.get(key,[])),len(records.get(key,[]))
                frozen = cs[key]
                ctx.check('methods_record_count_'+'_'.join(key),total,integer(frozen['marker_records']))
                ctx.check('methods_unique_count_'+'_'.join(key),unique,integer(frozen['unique_positions']))
                if quota is not None:
                    ctx.check('methods_available_quota_positions_'+'_'.join(key),unique,integer(qs[key]['available_unique_positions']))
                    ctx.check('methods_minimum_quota_'+'_'.join(key),quota,integer(qs[key]['quota']))
                detail.append(dict(resource=resource,arm=arm,chromosome=c,marker_records=total,
                                   unique_positions=unique,duplicate_position_records=total-unique,
                                   chromosome_quota=quota,quota_set_by_resource=unique==quota if quota is not None else None))
    summary=[]
    for resource in PANELS:
        for arm in ('A','B','C'):
            total=sum(len(records.get((resource,arm,str(c)),[])) for c in range(1,11))
            unique=sum(len(positions.get((resource,arm,str(c)),[])) for c in range(1,11))
            summary.append(dict(resource=resource,arm=arm,marker_records=total,unique_chromosome_positions=unique,
                                duplicate_position_records=total-unique,matched_positions_per_draw=sums[arm] if arm in ('A','C') else None))
    ctx.write_csv('methods_chromosome_counts_and_quotas.csv',detail)
    ctx.write_csv('methods_count_summary.csv',summary)
    ctx.note('Quotas use unique chromosome-position pairs separately in A and C. B is shown only for reconciliation; no B subsampling quota is defined. Duplicate original tag-origin coordinates can collapse multiple marker records.')
    return summary


def composition(ctx, lookup):
    inventory=ctx.read('source_data/04_marker_identity.xlsx','marker_inventory')
    selected={r:[x for x in inventory if x['resource']==r] for r in NACIONAL}
    ids={r:{str(x['marker_id']) for x in selected[r] if x.get('marker_id') not in (None,'')} for r in NACIONAL}
    canonical={r:{str(x['canonical_id']) for x in selected[r] if x.get('canonical_id') not in (None,'')} for r in NACIONAL}
    summary=[]
    for label,sets in [('literal_marker_id',ids),('canonical_marker_id',canonical)]:
        a,b=sets[NACIONAL[0]],sets[NACIONAL[1]]
        summary.append(dict(identifier_definition=label,G7_source_records=len(selected[NACIONAL[0]]),
                            MAF5_source_records=len(selected[NACIONAL[1]]),G7_unique_ids=len(a),MAF5_unique_ids=len(b),
                            shared_unique_ids=len(a&b),G7_only_unique_ids=len(a-b),MAF5_only_unique_ids=len(b-a),
                            G7_is_subset_of_MAF5=a<=b,MAF5_is_subset_of_G7=b<=a,
                            interpretation='Identifier overlap does not establish sequence or allele equivalence.'))
    overlaps=ctx.read('source_data/04_marker_identity.xlsx','id_overlap')
    baseline=[r for r in overlaps if {r['resource_a'],r['resource_b']}==set(NACIONAL)]
    if len(baseline)!=1:
        raise ValueError('Expected one frozen Nacional ID-overlap comparison')
    ctx.check('methods_literal_ID_overlap',len(ids[NACIONAL[0]] & ids[NACIONAL[1]]),integer(baseline[0]['shared_raw_ids']))
    ctx.check('methods_canonical_ID_overlap',len(canonical[NACIONAL[0]] & canonical[NACIONAL[1]]),integer(baseline[0]['shared_canonical_ids']))
    detail=[]
    for resource in NACIONAL:
        other=NACIONAL[1] if resource==NACIONAL[0] else NACIONAL[0]
        for r in selected[resource]:
            row=integer(r['source_row']); ledger=lookup.get((resource,row))
            if ledger is None:
                raise ValueError('Nacional source marker missing from cohort ledger: '+repr((resource,row)))
            marker,canon=str(r['marker_id']),str(r['canonical_id'])
            detail.append(dict(resource=resource,source_row=row,marker_id=marker,canonical_id=canon,
                               reference_key=r.get('reference_key'),literal_id_in_other_export=marker in ids[other],
                               canonical_id_in_other_export=canon in canonical[other],
                               in_A=truth(ledger['in_A']),in_B=truth(ledger['in_B']),in_C=truth(ledger['in_C']),
                               original_chromosome=ledger['original_chromosome'],original_position_bp=ledger['original_position_bp'],
                               corrected_chromosome=ledger['reference_chromosome'],corrected_snp_bp=ledger['reference_snp_bp'],
                               cohort_reason=ledger['cohort_reason']))
    ctx.write_csv('methods_marker_set_overlap.csv',summary)
    ctx.write_csv('methods_nacional_marker_composition.csv',detail)
    ctx.note('Shared sample IDs do not imply nested marker sets. Literal and canonical marker-ID overlap are reported separately; neither is independent sequence or allele validation.')
    return summary,ids,canonical,inventory


def gene_distance(ps,start,end):
    i=bisect_left(ps,start)
    if i<len(ps) and ps[i]<=end:
        return 0
    ds=([start-ps[i-1]] if i else [])+([ps[i]-end] if i<len(ps) else [])
    return min(ds) if ds else None


def monoterpene(ctx,positions,records,ids,canonical,inventory):
    genes=ctx.read('alignment_coverage/03_fixed_query_sets.xlsx','unique_candidate_genes')
    members=ctx.read('alignment_coverage/03_fixed_query_sets.xlsx','pathway_membership')
    gene_ids={r['gene_id'] for r in members if r['pathway']=='Monoterpene'}
    selected=[r for r in genes if r['gene_id'] in gene_ids]
    ctx.check('methods_monoterpene_gene_denominator',len(gene_ids),11)
    ctx.check('methods_monoterpene_coordinate_count',len(selected),len(gene_ids))
    ctx.check('methods_monoterpene_unique_coordinate_ids',len({r['gene_id'] for r in selected}),len(selected))
    meta={(r['resource'],integer(r['source_row'])):r for r in inventory if r['resource'] in PANELS}
    coverage,hits=[],[]
    for gene in sorted(selected,key=lambda r:r['gene_id']):
        c,start,end=chrom(gene['chromosome']),integer(gene['start_bp']),integer(gene['end_bp'])
        if c is None or start is None or end is None or start<1 or end<start:
            raise ValueError('Invalid monoterpene gene interval: '+repr(gene))
        for resource in PANELS:
            other=NACIONAL[1] if resource==NACIONAL[0] else NACIONAL[0] if resource==NACIONAL[1] else None
            for arm in ('A','B','C'):
                key=(resource,arm,c);ps=positions.get(key,[])
                distance=gene_distance(ps,start,end)
                hitpos=set(ps[bisect_left(ps,start-50000):bisect_right(ps,end+50000)])
                matching=[(p,r) for p,r in records.get(key,[]) if p in hitpos]
                exclusive=[]
                canonical_exclusive=[]
                for p,r in matching:
                    sr=integer(r['original_row']);source=meta[(resource,sr)]
                    shared=str(r['marker_id']) in ids[other] if other else None
                    cshared=str(source['canonical_id']) in canonical[other] if other else None
                    if other and not shared:
                        exclusive.append(r['marker_id'])
                    if other and not cshared:
                        canonical_exclusive.append(str(source['canonical_id']))
                    dist=0 if start<=p<=end else start-p if p<start else p-end
                    hits.append(dict(gene_id=gene['gene_id'],resource=resource,arm=arm,gene_chromosome=c,
                                     gene_start_bp=start,gene_end_bp=end,marker_record_key=r['record_key'],
                                     marker_id=r['marker_id'],source_row=sr,canonical_id=source['canonical_id'],
                                     source_tag_key=r['source_tag_key'],marker_position_bp=p,distance_to_gene_bp=dist,
                                     literal_id_in_other_Nacional_export=shared,canonical_id_in_other_Nacional_export=cshared))
                coverage.append(dict(gene_id=gene['gene_id'],chromosome=c,start_bp=start,end_bp=end,
                                     resource=resource,arm=arm,nearest_marker_distance_bp=distance,
                                     recovered_50kb=distance is not None and distance<=50000,
                                     marker_records_within_50kb=len(matching),unique_positions_within_50kb=len(hitpos),
                                     hit_literal_ids_absent_from_other_Nacional_export=';'.join(sorted(set(exclusive))),
                                     hit_canonical_ids_absent_from_other_Nacional_export=';'.join(sorted(set(canonical_exclusive)))))
    frozen=ctx.read('alignment_coverage/04_paired_coverage_effects.xlsx','summary')
    fs={(r['resource'],r['arm']):r for r in frozen if r['resource'] in PANELS and r['query_set']=='pathway:Monoterpene' and integer(r['flank_bp'])==50000}
    summary=[]
    for resource in PANELS:
        for arm in ('A','B','C'):
            rr=[r for r in coverage if r['resource']==resource and r['arm']==arm]
            n=sum(r['recovered_50kb'] for r in rr); f=fs[(resource,arm)]
            ctx.check('methods_monoterpene_denominator_'+resource+'_'+arm,len(rr),integer(f['n_queries']))
            ctx.check('methods_monoterpene_recovered_'+resource+'_'+arm,n,integer(f['n_recovered']))
            summary.append(dict(resource=resource,arm=arm,flank_bp=50000,n_genes=len(rr),n_recovered=n,fraction_recovered=n/len(rr)))
    cv={(r['gene_id'],r['resource'],r['arm']):r for r in coverage}
    differences=[]
    for gid in sorted(gene_ids):
        for arm in ('A','B','C'):
            g,m=cv[(gid,NACIONAL[0],arm)],cv[(gid,NACIONAL[1],arm)]
            differences.append(dict(gene_id=gid,arm=arm,G7_covered_50kb=g['recovered_50kb'],MAF5_covered_50kb=m['recovered_50kb'],
                                    G7_nearest_distance_bp=g['nearest_marker_distance_bp'],MAF5_nearest_distance_bp=m['nearest_marker_distance_bp'],
                                    G7_hit_canonical_ids_absent_from_MAF5=g['hit_canonical_ids_absent_from_other_Nacional_export'],
                                    MAF5_hit_canonical_ids_absent_from_G7=m['hit_canonical_ids_absent_from_other_Nacional_export']))
    ctx.write_csv('methods_monoterpene_coverage.csv',coverage)
    ctx.write_csv('methods_monoterpene_covering_marker_records.csv',hits)
    ctx.write_csv('methods_monoterpene_summary.csv',summary)
    ctx.write_csv('methods_monoterpene_G7_MAF5_comparison.csv',differences)
    ctx.note('Monoterpene coverage uses the eleven explicitly assigned genes and extends each gene boundary by 50 kb. Every covering marker is listed. These are physical-coverage/composition results, not LD, functional validation, or transfer of aroma effects.')
    ctx.note('Monoterpene between-export absence is summarized using canonical marker IDs. Literal ID absence is shown separately because naming-format differences can make literal intersections empty.')
    return summary


def run(ctx):
    tags=tag_audit(ctx)
    positions,records,lookup=ledger_index(ctx)
    counts=quota_audit(ctx,positions,records)
    overlap,ids,canonical,inventory=composition(ctx,lookup)
    mono=monoterpene(ctx,positions,records,ids,canonical,inventory)
    ctx.summary['methods_audit']=dict(tag_lengths=tags,marker_count_summary=counts,
                                     marker_set_overlap=overlap,monoterpene_coverage=mono,
                                     scope='Tag lengths, unique-position quotas and physical marker composition.')
