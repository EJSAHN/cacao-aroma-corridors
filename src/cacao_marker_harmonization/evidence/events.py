"""Evidence-preserving association crosswalks; no silent row merging or unit conversion."""
from __future__ import annotations
from collections import Counter, defaultdict
from itertools import permutations
import re
import unicodedata
from .common import DataError, close_number, stable_id


def lexical_trait(s):
    """Only normalize whitespace around a terminal UR/R suffix; preserve trait identity."""
    s=unicodedata.normalize('NFC',s).strip()
    return re.sub(r'\s*_\s*(UR|R)$',lambda m:'_'+m[1],s)


def event_tuple(r, trait=None):
    return (r['chromosome'],r['peak_bp'],r['trait'] if trait is None else trait,
            r['method_reported'],r['marker_filter_reported'])


def interval_tuple(r):
    return (r['source_start_bp'],r['source_end_bp'],r['source_block_number'])


def interval_state(r, chromosomes):
    s,e,_=interval_tuple(r);p=r['peak_bp']
    if r['chromosome'] not in chromosomes or p is None or p<=0:return 'invalid_peak'
    if s is None and e is None:return 'no_boundaries'
    if s is None or e is None:return 'partial_boundaries'
    if s<=0 or e<=0:return 'nonpositive_boundary'
    if s>e:return 'reversed_boundaries'
    if not s<=p<=e:return 'peak_outside_interval'
    return 'valid_as_labeled'


def scalar_comparison(a,b,factors,rtol):
    if a is None or b is None:return 'unavailable',None
    if close_number(a,b,rtol):return 'equal',1.0
    if b==0:return 'different_zero_denominator',None
    for factor in factors:
        if factor!=1 and close_number(a,b*factor,rtol):
            return 'scale_factor_candidate',factor
    return 'different',a/b


def validated_aliases(cfg):
    aliases={}
    for item in cfg.get('trait_alias_candidates',[]):
        a=item.get('source','');b=item.get('target','');reason=item.get('reason','')
        if not a or not b or not reason or a==b or a in aliases:
            raise DataError('Invalid or duplicate candidate trait alias.')
        if item.get('apply_to_primary',False):
            raise DataError('This evidence stage does not apply aliases to primary analysis.')
        aliases[a]=(b,reason)
    if any(b in aliases for b,_ in aliases.values()):
        raise DataError('Chained/cyclic candidate aliases are not permitted.')
    return aliases


def reconcile(records,cfg):
    chromosomes=set(cfg['chromosomes']);rtol=cfg.get('numeric_comparison_rtol',1e-9)
    factors=cfg.get('reported_variance_scale_candidates',[1.0,100.0,0.01])
    master=[r for r in records if r['in_master']]
    valid_peak=lambda r:r['chromosome'] in chromosomes and r['peak_bp'] is not None and r['peak_bp']>0 and bool(r['trait'])
    master=[r for r in master if valid_peak(r)]
    if not master:raise DataError('Association master contains no valid peaks.')
    exact=defaultdict(list);lex=defaultdict(list)
    for m in master:
        exact[event_tuple(m)].append(m);lex[event_tuple(m,lexical_trait(m['trait']))].append(m)
    aliases=validated_aliases(cfg)
    cross=[];payload=[];intervals=[];notes=[]
    summary=defaultdict(Counter)
    valid_master_intervals=defaultdict(set)
    for m in master:
        if interval_state(m,chromosomes)=='valid_as_labeled':valid_master_intervals[event_tuple(m)].add(interval_tuple(m))
    for r in records:
        candidate=lexical_trait(r['trait']);basis='';alias_reason=''
        options=exact.get(event_tuple(r),[])
        if options:basis='exact_event_key'
        if not options:
            options=lex.get(event_tuple(r,candidate),[])
            if options:basis='terminal_suffix_whitespace_candidate'
        if not options and candidate in aliases:
            candidate,alias_reason=aliases[candidate]
            options=lex.get(event_tuple(r,candidate),[])
            if options:basis='configured_trait_alias_candidate'
        state='master' if r['in_master'] else 'matched_unique_master_event' if len(options)==1 else 'ambiguous_master_event' if len(options)>1 else 'not_matched'
        m=options[0] if len(options)==1 else None
        ps,pf=scalar_comparison(r['p_value_reported'],m['p_value_reported'] if m else None,[1],rtol)
        vs,vf=scalar_comparison(r['explained_variance_reported'],m['explained_variance_reported'] if m else None,factors,rtol)
        raw=interval_tuple(r);target=interval_tuple(m) if m else None
        exact_interval=raw==target if m else None
        transforms=[]
        if target and all(v is not None for v in raw) and raw!=target:
            for perm in permutations(range(3)):
                if tuple(raw[i] for i in perm)==target:transforms.append(','.join(str(i+1) for i in perm))
        relation='same_labeled_fields' if exact_interval else 'field_permutation_matches_master' if transforms else 'different_or_unavailable'
        source_state=interval_state(r,chromosomes)
        review='NO_PRIMARY_CHANGE'
        if not r['in_master'] and m and ps=='equal':
            review='DUPLICATE_EVENT_CANDIDATE_WITH_PAYLOAD_REVIEW' if relation!='same_labeled_fields' or vs!='equal' else 'DUPLICATE_EVENT_CANDIDATE'
        row=dict(source_sheet=r['source_sheet'],source_row=r['source_row'],trait_original=r['trait'],trait_candidate=candidate,
            chromosome=r['chromosome'],peak_bp=r['peak_bp'],method_reported=r['method_reported'],marker_filter_reported=r['marker_filter_reported'],
            source_event_key=r['event_key'],match_status=state,matching_basis=basis,
            master_row=m['source_row'] if m else None,master_trait=m['trait'] if m else '',master_matches=len(options),
            p_comparison=ps,variance_comparison=vs,source_to_master_variance_factor=vf,
            interval_relation=relation,field_permutation=';'.join(transforms),interval_state=source_state,
            alias_reason=alias_reason,action=review)
        cross.append(row)
        summary[r['source_sheet']][(state,basis,ps,vs,relation)]+=1
        if not r['in_master'] and (basis!='exact_event_key' or relation!='same_labeled_fields' or vs!='equal' or ps!='equal'):
            payload.append({**row,'source_start_bp':raw[0],'source_end_bp':raw[1],'source_block_number':raw[2],
                'master_start_bp':target[0] if target else None,'master_end_bp':target[1] if target else None,'master_block_number':target[2] if target else None,
                'source_p_value':r['p_value_reported'],'master_p_value':m['p_value_reported'] if m else None,
                'source_explained_variance':r['explained_variance_reported'],
                'master_explained_variance':m['explained_variance_reported'] if m else None})
        if source_state!='valid_as_labeled':
            geometric=set()
            if valid_peak(r) and all(v is not None for v in raw):
                for perm in permutations(raw):
                    if 0<perm[0]<=r['peak_bp']<=perm[1]:geometric.add(perm)
            intervals.append(dict(source_sheet=r['source_sheet'],source_row=r['source_row'],chromosome=r['chromosome'],peak_bp=r['peak_bp'],
                trait=r['trait'],state=source_state,source_start_bp=raw[0],source_end_bp=raw[1],source_block_number=raw[2],
                master_row=m['source_row'] if m else None,master_start_bp=target[0] if target else None,master_end_bp=target[1] if target else None,
                geometry_only_permutations='; '.join(str(x) for x in sorted(geometric)),
                action='REVIEW_ONLY_NO_BOUNDARIES_IMPUTED'))
    for k,group in exact.items():
        if len(group)>1:
            notes.append(dict(issue='duplicate_master_event_key',event_key=stable_id('e_',k),source_rows=';'.join(str(r['source_row']) for r in group)))
    countrows=[dict(source_sheet=sheet,match_status=s[0],matching_basis=s[1],p_comparison=s[2],variance_comparison=s[3],interval_relation=s[4],records=n)
               for sheet,co in sorted(summary.items()) for s,n in sorted(co.items())]
    counts=[dict(metric='source_records',count=len(records)),dict(metric='master_records',count=len(master)),
            dict(metric='master_unique_exact_events',count=len(exact)),
            dict(metric='master_unique_peak_positions',count=len({(r['chromosome'],r['peak_bp']) for r in master})),
            dict(metric='all_sheet_unique_peak_positions',count=len({(r['chromosome'],r['peak_bp']) for r in records if valid_peak(r)})),
            dict(metric='extra_exact_event_records',count=sum(not r['in_master'] and event_tuple(r) not in exact for r in records)),
            dict(metric='extra_records_matched_by_candidates',count=sum(not r['source_sheet']==cfg['association_master_sheet'] and r['matching_basis'] in ('terminal_suffix_whitespace_candidate','configured_trait_alias_candidate') and r['master_matches']==1 for r in cross)),
            dict(metric='master_boundaries_not_valid_as_labeled',count=sum(interval_state(r,chromosomes)!='valid_as_labeled' for r in master))]
    return dict(counts=counts,crosswalk=cross,payload_differences=payload,interval_review=intervals,match_summary=countrows,master_ambiguities=notes)
