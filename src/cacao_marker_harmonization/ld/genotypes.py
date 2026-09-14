"""Strict genotype decoding, fixed QC and shared-sample sensitivity comparisons."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
import re
import numpy as np
from .common import DataError, text, integer
from .readers import Matrix


@dataclass
class Panel:
    name: str
    samples: list[str]
    markers: list[dict]
    dosage: np.ndarray
    groups: dict[str,np.ndarray]
    metadata: dict


def decode_matrix(matrix: Matrix, cfg: dict):
    e = cfg['encoding']; codes = e['symbol_to_dosage']; missing = set(e['missing_symbols'])
    dosage = np.full((len(matrix.markers),len(matrix.samples)), np.nan, dtype=np.float32)
    states = []; counts = Counter(); invalid = []
    for i, (m, calls) in enumerate(zip(matrix.markers,matrix.calls)):
        declared = set(re.findall('[ACGT]', text(m['alleles_reported']).upper()))
        status = 'OK' if declared == set(e['declared_allele_set']) else 'UNSUPPORTED_DECLARED_ALLELES'
        unknown = sorted(set(calls) - set(codes) - missing)
        if unknown: status = 'UNSUPPORTED_SYMBOLS'
        for s, token in zip(matrix.samples,calls):
            counts[token] += 1
            if token not in codes and token not in missing:
                invalid.append(dict(resource=matrix.name, source_row=m['source_row'], marker_id=m['marker_id'],sample_id=s,symbol=token))
        if status == 'OK':
            dosage[i] = [codes.get(t,np.nan) for t in calls]
        states.append(status)
    return dosage, np.array(states,dtype=object), [dict(resource=matrix.name,symbol=k,count=int(v),meaning='dosage '+str(codes[k]) if k in codes else 'missing' if k in missing else 'unsupported') for k,v in sorted(counts.items())], invalid


def row_qc(d: np.ndarray, q: dict):
    n = np.isfinite(d).sum(axis=1)
    rate = n/max(1,d.shape[1])
    total = np.nansum(d,axis=1,dtype=np.float64)
    p = np.divide(total,2*n,out=np.full(len(d),np.nan),where=n>0)
    maf = np.minimum(p,1-p)
    sumsquare = np.nansum(d.astype(np.float64)**2,axis=1)
    variance = n*sumsquare-total**2
    keep = (rate >= q['marker_min_call_rate']) & (maf >= q['min_maf']) & (variance > 0)
    return dict(called=n,call_rate=rate,maf=maf,variable=variance>0,keep=keep)


def raw_concordance(a: Matrix, b: Matrix, cfg: dict):
    """Compare source symbols, not named biological allele polarity; no calls are changed."""
    def lookup(m):
        d = {r['canonical_id']:i for i,r in enumerate(m.markers)}
        if len(d)!=len(m.markers): raise DataError('Canonical source marker IDs are not unique: '+m.name)
        return d
    aa,bb=lookup(a),lookup(b)
    ss=sorted(set(a.samples)&set(b.samples));mm=sorted(set(aa)&set(bb))
    ai={s:i for i,s in enumerate(a.samples)};bi={s:i for i,s in enumerate(b.samples)}
    missing=set(cfg['encoding']['missing_symbols']); flags=set(cfg['sensitivity_sample_ids'])
    totals=Counter();different=Counter();miss=Counter();details=[];marker_stats=[]
    for mid in mm:
        ma,mb=a.markers[aa[mid]],b.markers[bb[mid]];good=bad=unavailable=0
        for s in ss:
            x,y=a.calls[aa[mid]][ai[s]],b.calls[bb[mid]][bi[s]]
            if x in missing or y in missing: miss[s]+=1;unavailable+=1;continue
            totals[s]+=1;good+=1
            if x!=y:
                different[s]+=1;bad+=1
                details.append(dict(marker_key=mid,sample_id=s,resource_a=a.name,symbol_a=x,resource_b=b.name,symbol_b=y,
                    source_row_a=ma['source_row'],source_row_b=mb['source_row'],sample_in_declared_sensitivity=s in flags))
        marker_stats.append(dict(marker_key=mid,shared_nonmissing_calls=good,discordant_calls=bad,missing_pairs=unavailable,
                                 original_chromosome_a=ma['chromosome'],original_chromosome_b=mb['chromosome'],
                                 original_position_a=ma['position_bp'],original_position_b=mb['position_bp']))
    rows=[dict(sample_id=s,shared_nonmissing_calls=totals[s],discordant_calls=different[s],missing_pairs=miss[s],
               concordant_fraction=(totals[s]-different[s])/totals[s] if totals[s] else None,
               selected_for_sensitivity=s in flags) for s in ss]
    stats=dict(resource_a=a.name,resource_b=b.name,samples_a=len(a.samples),samples_b=len(b.samples),
               shared_samples=len(ss),shared_markers=len(mm),shared_nonmissing_calls=sum(totals.values()),
               discordant_calls=sum(different.values()),samples_with_discordance=';'.join(s for s in ss if different[s]),
               discordance_outside_declared_samples=sum(v for s,v in different.items() if s not in flags))
    membership=[dict(sample_id=s,in_resource_a=s in ai,in_resource_b=s in bi,
                     selected_for_sensitivity=s in flags) for s in sorted(set(a.samples)|set(b.samples))]
    return dict(summary=[stats],sample_concordance=rows,marker_concordance=marker_stats,discordant_cells=details,sample_overlap=membership)


def prepare_panel(matrix: Matrix, ledger: list[dict], cfg: dict):
    if matrix.errors: raise DataError('Source Excel error cells detected: '+matrix.name)
    n=len(matrix.markers)
    d,decoding,symbols,invalid=decode_matrix(matrix,cfg)
    ll={integer(r['original_row']):r for r in ledger if r['resource']==matrix.name}
    raw_by_row={m['source_row']:m for m in matrix.markers}
    if len(ll)!=sum(r['resource']==matrix.name for r in ledger) or set(ll)!=set(raw_by_row):
        raise DataError('Source rows and upstream marker ledger differ: '+matrix.name)
    keys=[m['canonical_id'] for m in matrix.markers]
    if len(set(keys))!=n:raise DataError('Ambiguous canonical genotype identifiers: '+matrix.name)
    selected=set(cfg['sensitivity_sample_ids'])
    if not selected.issubset(set(matrix.samples)):
        raise DataError('Configured sensitivity sample missing from '+matrix.name+': '+str(sorted(selected-set(matrix.samples))))
    mapped=np.zeros(n,dtype=bool); coord_counts=Counter(); combined=[]
    for i,m in enumerate(matrix.markers):
        r=ll[m['source_row']]
        if (m['marker_id']!=r['marker_id'] or text(m['chromosome'])!=text(r['original_chromosome'])
                or m['position_bp']!=integer(r['original_position_bp'])):
            raise DataError('Source marker/coordinate does not match alignment ledger: '+matrix.name+':'+str(m['source_row']))
        ok=r['in_C']
        if ok and (text(r['reference_chromosome']) not in cfg['chromosomes'] or integer(r['reference_snp_bp']) is None or integer(r['reference_snp_bp'])<1):
            raise DataError('Accepted marker has invalid corrected coordinates.')
        if ok: coord_counts[(text(r['reference_chromosome']),integer(r['reference_snp_bp']))]+=1
        mapped[i]=ok
        combined.append(dict(resource=matrix.name,marker_key=m['canonical_id'],marker_id=m['marker_id'],record_key=r['record_key'],
            original_row=m['source_row'],original_chromosome=m['chromosome'],original_position_bp=m['position_bp'],
            chromosome=text(r['reference_chromosome']),position_bp=integer(r['reference_snp_bp']),
            source_tag_key=text(r['source_tag_key']),identity_basis=text(r['identity_basis']),source_matrix_index=i))
    dec=decoding=='OK'; valid=mapped&dec
    if not valid.any():raise DataError('No mapped decodable marker in '+matrix.name)
    initial_rate=np.isfinite(d[valid]).mean(axis=0)
    keep_samples=initial_rate>=cfg['qc']['sample_min_call_rate']
    groups={'all_qc_samples':np.flatnonzero(keep_samples),
            'omit_flagged_samples':np.array([i for i,s in enumerate(matrix.samples) if keep_samples[i] and s not in selected],dtype=np.int32)}
    if any(len(ids)<cfg['qc']['min_pair_samples'] for ids in groups.values()):
        raise DataError('Too few samples after fixed sensitivity/sample QC: '+matrix.name)
    stat={k:row_qc(d[:,v],cfg['qc']) for k,v in groups.items()}
    duplicate=np.array([mapped[i] and coord_counts[(r['chromosome'],r['position_bp'])]>1 for i,r in enumerate(combined)])
    keep=valid&~duplicate
    for q in stat.values():keep &= q['keep']
    sq=[]
    for i,s in enumerate(matrix.samples):
        calls=d[valid,i]; nc=int(np.isfinite(calls).sum())
        sq.append(dict(resource=matrix.name,sample_id=s,initial_call_rate=float(initial_rate[i]),
            genotype_qc_pass=bool(keep_samples[i]),selected_for_sensitivity=s in selected,
            in_all_qc_samples=bool(keep_samples[i]),in_omit_flagged_samples=bool(keep_samples[i] and s not in selected),
            marker_denominator=int(valid.sum()),heterozygous_call_fraction=float(np.sum(calls==1)/nc) if nc else None))
    mq=[]
    for i,r in enumerate(combined):
        reasons=[]
        if not mapped[i]:reasons.append('NOT_IN_REFINED_FIXED_COHORT')
        if not dec[i]:reasons.append(str(decoding[i]))
        if duplicate[i]:reasons.append('DUPLICATE_CORRECTED_POSITION_ALL_EXCLUDED')
        row=dict(r,decoding_status=str(decoding[i]),mapped_fixed_cohort=bool(mapped[i]),duplicate_corrected_position=bool(duplicate[i]))
        for g,q in stat.items():
            row[g+'_call_rate']=float(q['call_rate'][i]);row[g+'_maf']=float(q['maf'][i]) if np.isfinite(q['maf'][i]) else None
            row[g+'_variable']=bool(q['variable'][i])
            if not q['keep'][i]:reasons.append('QC_FAIL_'+g.upper())
        row['retained_for_fixed_ld']=bool(keep[i]);row['exclusion_reasons']=';'.join(reasons)
        mq.append(row)
    rank={c:i for i,c in enumerate(cfg['chromosomes'])}
    ids=sorted(np.flatnonzero(keep), key=lambda i:(rank[combined[i]['chromosome']],combined[i]['position_bp'],combined[i]['marker_key']))
    if len(ids)<cfg['qc']['min_markers']:raise DataError('Too few fixed markers after QC: '+matrix.name)
    markers=[combined[i] for i in ids]
    panel=Panel(matrix.name,matrix.samples,markers,d[ids].copy(),groups,
                dict(source_markers=n,mapped_markers=int(mapped.sum()),fixed_markers=len(ids),
                     samples_all=len(groups['all_qc_samples']),samples_after_omission=len(groups['omit_flagged_samples']),
                     engine=matrix.engine,zero_filled_missing_calls=False))
    return panel,dict(marker_qc=mq,sample_qc=sq,symbol_counts=symbols,unsupported_calls=invalid)


def prepare_common(a: Panel,b: Panel,cfg:dict):
    ai={s:i for i,s in enumerate(a.samples)};bi={s:i for i,s in enumerate(b.samples)}
    asa={a.samples[i] for i in a.groups['all_qc_samples']};bsb={b.samples[i] for i in b.groups['all_qc_samples']}
    shared=sorted(asa&bsb);omit=sorted(set(shared)-set(cfg['sensitivity_sample_ids']))
    if len(omit)<cfg['qc']['min_pair_samples']:raise DataError('Insufficient shared samples for the fixed comparison.')
    ga={'shared_all':np.array([ai[s] for s in shared]),'shared_omit_flagged':np.array([ai[s] for s in omit])}
    gb={'shared_all':np.array([bi[s] for s in shared]),'shared_omit_flagged':np.array([bi[s] for s in omit])}
    ma={r['marker_key']:i for i,r in enumerate(a.markers)};mb={r['marker_key']:i for i,r in enumerate(b.markers)}
    audit=[];ia=[];ib=[]
    for k in sorted(set(ma)&set(mb)):
        x,y=a.markers[ma[k]],b.markers[mb[k]]
        ok=(x['chromosome'],x['position_bp'],x['source_tag_key'])==(y['chromosome'],y['position_bp'],y['source_tag_key'])
        audit.append(dict(marker_key=k,coordinates_and_tag_agree=ok,chromosome_a=x['chromosome'],position_a=x['position_bp'],
                          chromosome_b=y['chromosome'],position_b=y['position_bp'],retained=False,
                          reason='PENDING_SHARED_SAMPLE_QC' if ok else 'COORDINATE_OR_IDENTITY_CONFLICT'))
        if ok:ia.append(ma[k]);ib.append(mb[k])
    ia=np.array(ia,dtype=np.int32);ib=np.array(ib,dtype=np.int32)
    if not len(ia):raise DataError('No unambiguous common markers.')
    keep=np.ones(len(ia),dtype=bool)
    for p,inds,g in [(a,ia,ga),(b,ib,gb)]:
        for sid in g.values():keep &= row_qc(p.dosage[inds][:,sid],cfg['qc'])['keep']
    kept_keys={a.markers[i]['marker_key'] for i in ia[keep]}
    for r in audit:
        if r['coordinates_and_tag_agree']:
            r['retained']=r['marker_key'] in kept_keys
            r['reason']='FIXED_SHARED_COMPARISON' if r['retained'] else 'SHARED_SAMPLE_QC_FAIL'
    ia,ib=ia[keep],ib[keep]
    order=sorted(range(len(ia)),key=lambda j:(cfg['chromosomes'].index(a.markers[ia[j]]['chromosome']),a.markers[ia[j]]['position_bp'],a.markers[ia[j]]['marker_key']))
    ia,ib=ia[order],ib[order]
    if len(ia)<cfg['qc']['min_markers']:raise DataError('Too few shared QC markers.')
    return Panel(a.name,a.samples,[a.markers[i] for i in ia],a.dosage[ia],ga,{}),Panel(b.name,b.samples,[b.markers[i] for i in ib],b.dosage[ib],gb,{}),audit,shared,omit
