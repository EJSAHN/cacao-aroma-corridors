"""Local sequence evidence for identifier fields and position conventions.

Matching against an exported tag is not a reference-genome alignment. None of
these diagnostics grants permission to replace source coordinates.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import re
from .common import dart_suffix, text, fraction


def dna(value, allow_empty=False):
    s=re.sub(r'\s+', '', text(value)).upper()
    return s if (s or allow_empty) and re.fullmatch('[ACGT]*',s) else None


def reverse_complement(s):
    return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def reconstruct(row):
    left=dna(row.get('five_flank'),True); right=dna(row.get('three_flank'),True)
    m=re.fullmatch(r'\s*([ACGT])\s*>\s*([ACGT])\s*',text(row.get('variation')).upper())
    if left is None or right is None or m is None or m[1]==m[2] or not (left or right):
        return None
    return dict(left=left,right=right,reference_base=m[1],alternate_base=m[2],
        ref_tag=left+m[1]+right, alt_tag=left+m[2]+right)


def compare_tag(parts, refrows):
    if not refrows:return 'missing_reference_name'
    if len(refrows)!=1:return 'duplicate_reference_name'
    seq=dna(refrows[0]['sequence'])
    if parts is None:return 'unsupported_flank_or_variation'
    if seq is None:return 'unsupported_reference_sequence'
    if seq==parts['ref_tag']:return 'exact_ref_tag'
    if seq==parts['alt_tag']:return 'exact_alt_tag'
    if seq==reverse_complement(parts['ref_tag']):return 'reverse_complement_ref_tag'
    if seq==reverse_complement(parts['alt_tag']):return 'reverse_complement_alt_tag'
    return 'sequence_mismatch'


def sequence_evidence(rows, refs):
    detail=[];summary=Counter();anomalies=[]
    for row in rows:
        p=reconstruct(row)
        refmatch=compare_tag(p,refs.get(row['reference_key'],[]))
        idmatch=compare_tag(p,refs.get(row['marker_key'],[]))
        summary[('reference_key',refmatch)]+=1;summary[('marker_key',idmatch)]+=1
        suffix=dart_suffix(row['reference_key'])
        r={**row,'marker_vs_reference_name_equal':row['marker_key']==row['reference_key'],
           'ref_lookup_status':refmatch,'marker_lookup_status':idmatch,
           'ref_sequence_rows':';'.join(str(x['source_row']) for x in refs.get(row['reference_key'],[])),
           'reference_sequence':refs[row['reference_key']][0]['sequence'] if len(refs.get(row['reference_key'],[]))==1 else '',
           'five_flank_length':len(p['left']) if p else None,'three_flank_length':len(p['right']) if p else None,
           'tag_length':len(p['ref_tag']) if p else None,'id_terminal_index':suffix,
           'suffix_equals_five_length':bool(p and suffix is not None and suffix==len(p['left'])),
           'suffix_equals_three_length':bool(p and suffix is not None and suffix==len(p['right'])),
           'sequence_supported_reference_key':refmatch=='exact_ref_tag'}
        detail.append(r)
    for name,rr in sorted(refs.items()):
        if len(rr)!=1 or any(dna(r['sequence']) is None for r in rr):
            anomalies.extend(dict(reference_key=name,**r,issue='duplicate_name' if len(rr)>1 else 'invalid_sequence') for r in rr)
    summaryrows=[dict(key_type=k,status=s,records=n,fraction=fraction(n,len(rows))) for (k,s),n in sorted(summary.items())]
    return detail,summaryrows,anomalies


def make_index(rows, field):
    result=defaultdict(list)
    for r in rows:
        k=text(r.get(field))
        if k:result[k].append(r)
    return result


def coordinate_relation(a,b,chromosomes):
    c1=a.get('chromosome');c2=b.get('chromosome');p1=a.get('position_bp');p2=b.get('position_bp')
    if c1 not in chromosomes or c2 not in chromosomes or p1 is None or p2 is None or p1<=0 or p2<=0:
        return 'unavailable'
    if c1!=c2:return 'different_chromosomes'
    return 'exact_coordinate' if p1==p2 else 'same_chromosome_offset'


MODEL_DESCRIPTIONS={
    'same_position':'P_nacional = P_diversity',
    'five_flank_start_1based':'P_nacional = P_diversity - length(5-flank)',
    'five_flank_start_0based':'P_nacional = P_diversity - length(5-flank) - 1',
    'suffix_forward_0':'Index equals length(5-flank); P_nacional = P_diversity - index',
    'suffix_reverse_0':'Index equals length(3-flank); P_nacional = P_diversity - (tag_length - index)',
    'suffix_forward_1':'Index equals length(5-flank)+1; P_nacional = P_diversity - (index-1)',
    'suffix_reverse_1':'Index equals length(3-flank)+1; P_nacional = P_diversity - (tag_length-index+1)',
}

def models(row):
    """All models are predeclared hypotheses, not an automatically selected mapping."""
    p=row.get('position_bp');left=row.get('five_flank_length');right=row.get('three_flank_length')
    idx=row.get('id_terminal_index');length=row.get('tag_length')
    result={k:None for k in MODEL_DESCRIPTIONS}
    if p is None or p<=0:return result
    result['same_position']=p
    if left is not None:
        result['five_flank_start_1based']=p-left
        result['five_flank_start_0based']=p-left-1
    if idx is not None and length is not None:
        if idx==left:result['suffix_forward_0']=p-idx
        if idx==right:result['suffix_reverse_0']=p-(length-idx)
        if idx==left+1:result['suffix_forward_1']=p-(idx-1)
        if idx==right+1:result['suffix_reverse_1']=p-(length-idx+1)
    return result


def coordinate_hypotheses(diversity, matrices, cfg):
    details=[];summary=[];models_summary=[];ambiguities=[];candidates=[]
    chromosomes=set(cfg['chromosomes'])
    for name in cfg['nacional_resources']:
        target=matrices[name]
        ti=make_index(target.markers,'canonical_id')
        for ktype in ('marker_key','reference_key'):
            di=make_index(diversity,ktype)
            counts=Counter();mc=Counter()
            for k in sorted(set(di)&set(ti)):
                if len(di[k])!=1 or len(ti[k])!=1:
                    counts['ambiguous_key']+=1
                    ambiguities.append(dict(resource=name,key_type=ktype,key=k,diversity_rows=len(di[k]),target_rows=len(ti[k])))
                    continue
                a=di[k][0];b=ti[k][0];rel=coordinate_relation(a,b,chromosomes);counts[rel]+=1
                delta=b['position_bp']-a['position_bp'] if rel in ('same_chromosome_offset','exact_coordinate') else None
                predictions=models(a);matches=[]
                seq_supported=a['ref_lookup_status']=='exact_ref_tag' and ktype=='reference_key'
                for model,pred in predictions.items():
                    evaluable=seq_supported and delta is not None and pred is not None and pred>0
                    if evaluable:
                        mc[(model,'eligible')]+=1
                        if pred==b['position_bp']:
                            mc[(model,'matches')]+=1;matches.append(model)
                r=dict(resource=name,join_field=ktype,join_key=k,diversity_row=a['source_row'],target_row=b['source_row'],
                    diversity_marker_id=a['marker_id'],diversity_reference_name=a['reference_sequence_name'],target_marker_id=b['marker_id'],
                    diversity_chromosome=a['chromosome'],diversity_position_bp=a['position_bp'],
                    target_chromosome=b['chromosome'],target_position_bp=b['position_bp'],
                    coordinate_relation=rel,offset_target_minus_diversity_bp=delta,
                    five_flank_length=a['five_flank_length'],three_flank_length=a['three_flank_length'],tag_length=a['tag_length'],
                    id_terminal_index=a['id_terminal_index'],ref_lookup_status=a['ref_lookup_status'],marker_lookup_status=a['marker_lookup_status'],
                    suffix_equals_five_length=a['suffix_equals_five_length'],suffix_equals_three_length=a['suffix_equals_three_length'],
                    matched_hypotheses=';'.join(matches),source_variation=a['variation'],target_alleles=b['alleles_reported'],
                    target_strand=b['strand_reported'],source_version=a['assembly_reported'],legacy_remark=a['remark'])
                details.append(r)
                if ktype=='reference_key':
                    candidates.append(dict(resource=name,marker_id=b['marker_id'],diversity_source_row=a['source_row'],
                        proposed_identity_key=k,identity_evidence=a['ref_lookup_status'],
                        old_chromosome=b['chromosome'],old_position_bp=b['position_bp'],
                        diversity_chromosome=a['chromosome'],diversity_position_bp=a['position_bp'],
                        offset_bp=delta,coordinate_relation=rel,matched_hypotheses=';'.join(matches),
                        action='REVIEW_ONLY_NOT_APPLIED',
                        evidence_scope='Within-export tag/offset consistency only; assembly, strand and SNP coordinate require independent alignment.'))
            summary.append(dict(resource=name,join_field=ktype,shared_keys=sum(counts.values()),**{s:counts[s] for s in ('exact_coordinate','same_chromosome_offset','different_chromosomes','unavailable','ambiguous_key')}))
            if ktype=='reference_key':
                for model,description in MODEL_DESCRIPTIONS.items():
                    n=mc[(model,'eligible')];m=mc[(model,'matches')]
                    models_summary.append(dict(resource=name,hypothesis=model,definition=description,
                        eligible_same_chromosome_sequence_supported_pairs=n,exact_model_matches=m,match_fraction=fraction(m,n),
                        approved_coordinate_conversion=False))
    return details,summary,models_summary,ambiguities,candidates
