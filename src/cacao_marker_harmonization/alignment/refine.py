"""Resource projection, paired coverage and prespecified protein-coding sensitivity."""
from __future__ import annotations
from collections import Counter,defaultdict
from copy import deepcopy
from .common import DataError,text,integer
from ..coverage.numerics import Queries
from ..coverage.analysis import background_diagnostics


def project(inventory,decisions,cfg):
    byrow={r['source_row']:r for r in decisions};bykey=defaultdict(list)
    if len(byrow)!=len(decisions):raise DataError('Duplicate new mapping decision rows.')
    for r in decisions:
        if r['reference_key']:bykey[r['reference_key']].append(r)
    # Count keys over all source records, not just retained records.
    counts=Counter((r['resource'],text(r.get('canonical_id'))) for r in inventory)
    out=[]
    for old in inventory:
        resource=text(old['resource']);key=text(old.get('canonical_id'))
        row=dict(resource=resource,original_row=old['source_row'],original_marker_id=old['marker_id'],original_marker_key=key,
                 original_chromosome=old.get('chromosome'),original_position_bp=old.get('position_bp'),
                 accepted_for_coordinate_comparison=False,coordinate_status='NO_SEQUENCE_IDENTITY_BRIDGE',identity_basis='',
                 genotype_allele_polarity='NOT_ESTABLISHED')
        if resource=='Diversity':
            hit=byrow.get(integer(old['source_row']));possible=[hit] if hit else []
            if hit and hit['source_marker_id']!=old['marker_id']:raise DataError('Diversity source-row identity drift.')
            row['identity_basis']='own_flank_variation_tag_biallelic_alignment'
        elif resource in ('Nacional_G7','Nacional_MAF5'):
            possible=bykey.get(key,[]);row['identity_basis']='unique_reference_sequence_key_not_independent_sequence_validation'
        else:out.append(row);continue
        if len(possible)!=1:row['coordinate_status']='MISSING_OR_AMBIGUOUS_REFERENCE_KEY'
        else:
            m=possible[0];row['mapped_query_id']=m['query_id'];row['source_tag_key']=m['reference_key'];row['coordinate_status']=m['placement_status']
            if m.get('accepted_for_mapping'):
                row.update(reference_sequence_id=m['mapped_sequence_id'],reference_chromosome=m['mapped_chromosome'],reference_snp_bp=m['mapped_snp_bp'],
                    reference_tag_start_bp=m['tag_start_bp'],reference_tag_end_bp=m['tag_end_bp'],source_tag_strand_on_reference=m['mapped_strand'],
                    assembly_base=m['assembly_base'],other_source_base=m['other_source_base'],ref_cigar=m['ref_cigar'],alt_cigar=m['alt_cigar'],
                    accepted_for_coordinate_comparison=counts[resource,key]==1)
                row['coordinate_status']='RETAINED_DIRECT_ALIGNMENT' if resource=='Diversity' else 'RETAINED_PROJECTED_ALIGNMENT'
                if counts[resource,key]!=1:row['coordinate_status']='DUPLICATE_PANEL_MARKER_KEY_REVIEW'
        out.append(row)
    return out


def add_coding_query(queries):
    rows=queries['unique_candidate_genes'].rows
    coding=[r for r in rows if r['biotype']=='protein_coding']
    if not coding:raise DataError('No annotation-verified protein-coding candidate genes.')
    queries=dict(queries);name='protein_coding_candidate_genes';queries[name]=Queries(name,coding)
    ledger=[dict(gene_id=r['query_id'],biotype=r['biotype'],in_all_candidate_analysis=True,
                 in_protein_coding_sensitivity=r['biotype']=='protein_coding',
                 reason='INCLUDED_IN_BOTH_SETS' if r['biotype']=='protein_coding' else 'RETAINED_IN_ALL_CANDIDATES_EXCLUDED_ONLY_FROM_CODING_SENSITIVITY') for r in rows]
    return queries,ledger


def mapping_transitions(old_rows,new_rows):
    old={r['query_id']:r for r in old_rows};out=[]
    for r in new_rows:
        p=old[r['query_id']]
        prev=p['mapping_status']=='UNIQUE_EXACT_PLACEMENT' and bool(p.get('primary_chromosome'))
        now=bool(r.get('accepted_for_mapping'))
        same=prev and now and (p['mapped_sequence_id'],p['mapped_snp_bp'],p['mapped_strand'])==(r['mapped_sequence_id'],r['mapped_snp_bp'],r['mapped_strand'])
        status=('RETAINED_SAME_FOCAL_LOCATION' if same else 'RETAINED_DIFFERENT_FOCAL_LOCATION') if prev and now else 'NEWLY_RETAINED' if now else 'EXACT_ONLY_RETAINED_NOW_EXCLUDED' if prev else 'NOT_RETAINED_IN_EITHER_POLICY'
        out.append(dict(query_id=r['query_id'],reference_key=r['reference_key'],previous_status=p['mapping_status'],new_status=r['placement_status'],
                        transition=status,previous_snp_bp=p.get('mapped_snp_bp'),previous_sequence_id=p.get('mapped_sequence_id'),
                        new_sequence_id=r.get('mapped_sequence_id'),new_snp_bp=r.get('mapped_snp_bp')))
    return out


def background_sets(genes,all_ids,queries,lengths,points,cfg,log):
    combined=defaultdict(list)
    for name in cfg['background']['query_sets']:
        if name not in queries:raise DataError('Missing background query set: '+name)
        conf=deepcopy(cfg);conf['background']['query_set']=name
        result=background_diagnostics(genes,all_ids,queries,lengths,points,conf,log)
        for sheet,rows in result.items():
            if sheet=='background_gene_distances':
                # Reference distances do not change when the candidate target is
                # restricted. Store one copy, retaining every candidate exclusion.
                if not combined[sheet]:combined[sheet].extend(rows)
            elif sheet=='catalog_exclusions':
                if not combined[sheet]:combined[sheet].extend(rows)
            else:combined[sheet].extend(dict(target_set=name,**r) for r in rows)
    return dict(combined)


def compare_previous(old_summary,new_summary):
    def k(r):return r['resource'],r['arm'],r['query_set'],integer(r['flank_bp'])
    old={k(r):r for r in old_summary};out=[]
    if len(old)!=len(old_summary):raise DataError('Duplicate keys in exact-only coverage summary.')
    for r in new_summary:
        p=old.get(k(r))
        if p is None:continue
        if integer(p['n_queries'])!=integer(r['n_queries']):raise DataError('Common query denominator changed across mapping policies.')
        if r['arm']=='A' and (integer(p['n_recovered'])!=integer(r['n_recovered']) or integer(p['n_marker_records'])!=integer(r['n_marker_records'])):
            raise DataError('Original A arm changed relative to the exact-only policy.')
        out.append(dict(resource=r['resource'],arm=r['arm'],query_set=r['query_set'],flank_bp=r['flank_bp'],n_queries=r['n_queries'],
                        exact_only_marker_records=p['n_marker_records'],realignment_marker_records=r['n_marker_records'],
                        exact_only_recovered=p['n_recovered'],realignment_recovered=r['n_recovered'],
                        fraction_difference=r['fraction_recovered']-p['fraction_recovered'],
                        interpretation='Mapping-policy/subset sensitivity. Only B-C within one policy isolates coordinate correction.'))
    return out
