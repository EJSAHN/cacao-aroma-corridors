"""CODE2 genotype audit. Core QC is adapted from submitted v3.0.0 ld/genotypes.py.
No genotype is repaired; A/M/C dosage interpretation remains conditional.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
import re
import numpy as np
from xlsx_read import DataError
from math import isfinite
def text(v):
    return "" if v is None or (isinstance(v,float) and not isfinite(v)) else str(v).strip()
def integer(v):
    if v is None or isinstance(v,bool): return None
    try:
        n=float(v)
        return int(n) if isfinite(n) and n.is_integer() else None
    except (ValueError,TypeError,OverflowError): return None



@dataclass
class Panel:
    name: str
    samples: list[str]
    markers: list[dict]
    dosage: np.ndarray
    groups: dict[str,np.ndarray]
    metadata: dict


def decode_matrix(matrix, cfg: dict):
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


def raw_concordance(a, b, cfg: dict):
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


def prepare_panel(matrix, ledger: list[dict], cfg: dict):
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

@dataclass
class Matrix:
    name: str
    samples: list
    markers: list
    calls: list
    sheet: str
    engine: str
    errors: list


def canonical_marker(value):
    s = text(value)
    m = re.fullmatch(r'(\d+)\|([FRfr])\|(\d+)--(\d+)', s)
    if not m: m = re.fullmatch(r'(\d+)\.([FRfr])\.(\d+)\.(\d+)', s)
    return '%s|%s|%s--%s' % (m[1],m[2].upper(),m[3],m[4]) if m else s


def read_matrix(path, name):
    """Read literal cells, using the submitted restricted BIFF8 parser, no edits."""
    from biff8_read import BiffBook
    book = BiffBook(path)
    if len(book.sheets) != 1: raise DataError('Expected one frozen matrix sheet: '+name)
    sheet = next(iter(book.sheets)); it = iter(book.rows(sheet))
    rownum, header = next(it)
    header = [text(x) for x in header]
    if rownum != 1 or header[:5] != ['marker_id','alleles','chrom','pos','strand']:
        raise DataError('Unexpected frozen source headers: '+name)
    samples = header[5:]
    if not samples or any(not s for s in samples) or len(set(samples)) != len(samples):
        raise DataError('Missing or duplicate sample names: '+name)
    markers, calls = [], []
    for rownum, row in it:
        if len(row) > len(header) and any(text(x) for x in row[len(header):]):
            raise DataError('Source row exceeds header: '+name)
        row += [None] * max(0,len(header)-len(row))
        marker = text(row[0])
        if not marker:
            if any(text(x) for x in row): raise DataError('Nonempty row without marker: '+name)
            continue
        chrom = text(row[2]); ic = integer(chrom)
        chrom = str(ic) if ic is not None else chrom
        markers.append(dict(resource=name, source_file=path.name, source_sheet=sheet, source_row=rownum,
                            marker_id=marker, canonical_id=canonical_marker(marker), chromosome=chrom,
                            position_bp=integer(row[3]), chromosome_raw=text(row[2]), position_raw=text(row[3]),
                            alleles_reported=text(row[1]), strand_reported=text(row[4]), record_index=len(markers)))
        calls.append(tuple(text(x).upper() for x in row[5:]))
    if len({r['marker_id'] for r in markers}) != len(markers): raise DataError('Duplicate marker ID: '+name)
    return Matrix(name,samples,markers,calls,sheet,'restricted_biff8',[])


def _near(a, b):
    import math
    if isinstance(a,np.generic): a = a.item()
    if isinstance(b,np.generic): b = b.item()
    if a is None or b is None: return a is None and b is None
    if isinstance(a,bool) or isinstance(b,bool): return a == b
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
    return a == b


def compare_frozen(ctx, label, actual, expected, keys):
    """Compare by immutable keys, not worksheet order; no silent duplicate keys."""
    aa = {tuple(r[k] for k in keys):r for r in actual}
    ee = {tuple(r[k] for k in keys):r for r in expected}
    ctx.check(label+':unique_fresh_keys',len(aa),len(actual))
    ctx.check(label+':unique_frozen_keys',len(ee),len(expected))
    ctx.check(label+':keys',set(aa)==set(ee),True)
    differences = []
    for k in sorted(set(aa)&set(ee)):
        for f,v in ee[k].items():
            if not _near(aa[k].get(f),v):
                differences.append(dict(table=label,key='|'.join(map(str,k)),field=f,fresh=aa[k].get(f),frozen=v))
    ctx.check(label+':different_cells',len(differences),0)
    return differences


def sample_symbol_qc(matrix, indices, scope, missing, flags):
    """Literal M fraction is a conditional heterozygosity estimate only."""
    rows=[]
    for j,sample in enumerate(matrix.samples):
        counts=Counter(matrix.calls[i][j] for i in indices)
        total=sum(counts.values()); nm=sum(v for k,v in counts.items() if k in missing); called=total-nm
        rows.append(dict(resource=matrix.name,scope=scope,sample_id=sample,marker_denominator=total,
                         called=called,missing=nm,call_rate=called/total if total else None,
                         A_calls=counts['A'],M_calls=counts['M'],C_calls=counts['C'],
                         M_fraction_of_called=counts['M']/called if called else None,
                         conditional_heterozygous_fraction=counts['M']/called if called else None,
                         interpretation='Conditional on A/M/C = homozygote/heterozygote/homozygote; not independently verified',
                         selected_for_exclusion_sensitivity=sample in flags))
    return rows


def compare_panel_cells(a,b,scope):
    """Both state differences and missingness differences are retained explicitly."""
    ai={s:i for i,s in enumerate(a.samples)}; bi={s:i for i,s in enumerate(b.samples)}
    ma={r['marker_key']:i for i,r in enumerate(a.markers)}; mb={r['marker_key']:i for i,r in enumerate(b.markers)}
    samples=sorted(set(ai)&set(bi)); markers=sorted(set(ma)&set(mb)); rows=[]
    for sample in samples:
        x=a.dosage[[ma[k] for k in markers],ai[sample]]
        y=b.dosage[[mb[k] for k in markers],bi[sample]]
        both=np.isfinite(x)&np.isfinite(y); one=np.isfinite(x)^np.isfinite(y)
        bothmissing=~np.isfinite(x)&~np.isfinite(y)
        different=both&(x!=y)
        rows.append(dict(scope=scope,sample_id=sample,markers=len(markers),paired_nonmissing=int(both.sum()),
                         different_nonmissing_calls=int(different.sum()),one_missing_calls=int(one.sum()),
                         both_missing_calls=int(bothmissing.sum()),
                         exact_equal_including_missing=not bool(np.any(different|one))))
    return rows


def run(ctx):
    cfg=ctx.baseline_config['dosage_ld']; names=cfg['resources']
    flags=set(cfg['sensitivity_sample_ids']); missing=set(cfg['encoding']['missing_symbols'])
    ledger=ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx','marker_record_ledger')
    matrices={}; panels={}; qc={k:[] for k in ('marker_qc','sample_qc','symbol_counts','unsupported_calls')}
    literal=[]; declaration=[]; coding=[]; marker_literal=[]
    for name in names:
        print('[GENOTYPES] Reading literal source cells:',name,flush=True)
        m=read_matrix(ctx.source_paths[name],name); matrices[name]=m
        p, rr=prepare_panel(m,ledger,cfg); panels[name]=p
        marker_literal.extend(marker_symbol_qc(m,missing))
        for k in qc: qc[k].extend(rr[k])
        literal.extend(sample_symbol_qc(m,range(len(m.markers)),'all_source_markers',missing,flags))
        fixed_indices=[r['source_matrix_index'] for r in p.markers]
        literal.extend(sample_symbol_qc(m,fixed_indices,'fixed_LD_markers',missing,flags))
        counts=Counter(r['alleles_reported'] for r in m.markers)
        for label,count in sorted(counts.items()):
            declaration.append(dict(resource=name,source_file=ctx.source_paths[name].name,source_sheet=m.sheet,
                                    source_column='alleles',literal_declaration=label,marker_rows=count))
        coding.append(dict(resource=name,source_markers=len(m.markers),source_samples=len(m.samples),
                           literal_allele_declarations=';'.join(sorted(counts)),
                           source_symbols=';'.join(sorted({x for row in m.calls for x in row})),
                           declared_dosage='A=0;M=1;C=2',
                           coding_status='FORMAT_COMPATIBLE_CONDITIONAL_DOSAGE',
                           exact_exporter_codebook_confirmed=False,
                           biological_allele_polarity_established=False,
                           claim_supported='Export cells declare A/C and use A/M/C states; literal counts/QC reproducible',
                           claim_not_supported='Cell declarations alone do not prove biological heterozygote coding or sequencing-depth cause'))
        ctx.check('genotype:'+name+':unsupported_calls',len(rr['unsupported_calls']),0)
        ctx.check('genotype:'+name+':decoding_failures',sum(r['decoding_status']!='OK' for r in rr['marker_qc']),0)
    common_raw_ids=set(r['canonical_id'] for r in matrices[names[0]].markers)&set(r['canonical_id'] for r in matrices[names[1]].markers)
    for m in matrices.values():
        indices=[i for i,r in enumerate(m.markers) if r['canonical_id'] in common_raw_ids]
        literal.extend(sample_symbol_qc(m,indices,'common_raw_marker_ids',missing,flags))
    concord=raw_concordance(matrices[names[0]],matrices[names[1]],cfg)
    transitions=Counter((r['sample_id'],r['symbol_a'],r['symbol_b']) for r in concord['discordant_cells'])
    transition_rows=[dict(sample_id=k[0],symbol_G7=k[1],symbol_MAF5=k[2],discordant_cells=v) for k,v in sorted(transitions.items())]
    differences=[]
    frozen='dosage_ld/01_genotype_coordinate_qc.xlsx'
    for table, keys in (('marker_qc',['resource','original_row']),('sample_qc',['resource','sample_id']),
                        ('symbol_counts',['resource','symbol'])):
        differences.extend(compare_frozen(ctx,'frozen_'+table,qc[table],ctx.read(frozen,table),keys))
    for table, data, keys in (('raw_sample_concordance',concord['sample_concordance'],['sample_id']),
                              ('raw_discordant_cells',concord['discordant_cells'],['marker_key','sample_id'])):
        differences.extend(compare_frozen(ctx,'frozen_'+table,data,ctx.read(frozen,table),keys))
    a,b,common_audit,shared,omitted=prepare_common(panels[names[0]],panels[names[1]],cfg)
    common_panels={a.name:a,b.name:b}
    common_fixed_ids={r['marker_key'] for r in a.markers}
    for m in matrices.values():
        indices=[i for i,r in enumerate(m.markers) if r['canonical_id'] in common_fixed_ids]
        literal.extend(sample_symbol_qc(m,indices,'common_fixed_QC_marker_ids',missing,flags))
    differences.extend(compare_frozen(ctx,'frozen_shared_marker_audit',common_audit,
                        ctx.read('dosage_ld/04_shared_collection_comparison.xlsx','shared_marker_audit'),['marker_key']))
    cell_rows=compare_panel_cells(a,b,'fixed_common_QC_markers')
    totaldiff=sum(r['different_nonmissing_calls'] for r in cell_rows)
    afterdiff=sum(r['different_nonmissing_calls']+r['one_missing_calls'] for r in cell_rows if r['sample_id'] not in flags)
    ctx.check('raw_common_marker_count',concord['summary'][0]['shared_markers'],4444)
    ctx.check('raw_discordant_calls',concord['summary'][0]['discordant_calls'],535)
    ctx.check('raw_discordance_outside_declared_sensitivity',concord['summary'][0]['discordance_outside_declared_samples'],0)
    ctx.check('fixed_common_QC_discordant_calls',totaldiff,503)
    ctx.check('fixed_common_QC_state_and_missingness_differences_after_omission',afterdiff,0)
    for k,rows in qc.items():ctx.write_csv('genotype_'+k+'.csv',rows)
    for k,rows in concord.items():ctx.write_csv('genotype_raw_'+k+'.csv',rows)
    ctx.write_csv('genotype_literal_sample_QC.csv',literal)
    ctx.write_csv('genotype_literal_marker_QC.csv',marker_literal)
    ctx.write_csv('genotype_source_allele_declarations.csv',declaration)
    ctx.write_csv('genotype_coding_evidence.csv',coding)
    ctx.write_csv('genotype_discordance_transitions.csv',transition_rows)
    ctx.write_csv('genotype_frozen_comparison_differences.csv',differences)
    ctx.write_csv('genotype_common_marker_audit.csv',common_audit)
    ctx.write_csv('genotype_common_sample_cells.csv',cell_rows)
    ctx.summary['genotype_audit']=dict(
        coding_status='FORMAT_COMPATIBLE_CONDITIONAL_DOSAGE',exact_exporter_codebook_confirmed=False,source_cells_edited=False,
        primary_includes_flagged_samples=True,sensitivity_sample_ids=sorted(flags),
        raw_concordance=concord['summary'][0],fixed_common_markers=len(a.markers),
        fixed_common_nonmissing_discordances=totaldiff,
        fixed_common_state_and_missingness_differences_after_omission=afterdiff,
        frozen_comparison_differences=len(differences),
        panels={name:p.metadata for name,p in panels.items()},
        sequencing_depth_available=False,
        interpretation='Post-omission equality is algebraic agreement of shared exported cells, not independent biological replication.')
    ctx.note('A/M/C state counts, literal A/C allele declarations, and conditional M fractions are reported separately; no new codebook or biological coding validation was inferred.')
    ctx.note('SNA604 and L17H53 remain in primary analyses. Omission is sensitivity only. Missingness and conditional heterozygosity do not identify the cause of discordance; read depth is unavailable.')
    ctx.genotype_matrices=matrices;ctx.genotype_panels=panels;ctx.common_genotype_panels=common_panels
    return {'matrices':matrices,'panels':panels,'common_panels':common_panels}


def marker_symbol_qc(matrix, missing):
    """Descriptive frequencies conditional on A/M/C dosage; no HWE test."""
    rows=[]
    for marker,calls in zip(matrix.markers,matrix.calls):
        counts=Counter(calls); n=len(calls); nm=sum(v for k,v in counts.items() if k in missing); called=n-nm
        p=(counts['M']+2*counts['C'])/(2*called) if called else None
        h=counts['M']/called if called else None
        rows.append(dict(resource=matrix.name,source_row=marker['source_row'],marker_id=marker['marker_id'],
                         canonical_id=marker['canonical_id'],alleles_reported=marker['alleles_reported'],
                         sample_denominator=n,called=called,missing=nm,call_rate=called/n if n else None,
                         A_calls=counts['A'],M_calls=counts['M'],C_calls=counts['C'],
                         unsupported_calls=sum(v for k,v in counts.items() if k not in missing and k not in {'A','M','C'}),
                         M_fraction_of_called=h,conditional_heterozygous_fraction=h,
                         conditional_C_state_allele_frequency=p,
                         conditional_two_p_one_minus_p=2*p*(1-p) if p is not None else None,
                         minimum_ACM_class_count=min(counts['A'],counts['M'],counts['C']),
                         conditional_interpretation=True,HWE_test_performed=False))
    return rows
