"""Shared-sample cross-export coordinate-anchor correlation diagnostics.
No exported cells are merged; a position match is not established variant identity.
"""
from collections import Counter, defaultdict
import re
import numpy as np
from genotype_audit import Panel, row_qc

SCOPE='CONDITIONAL_A_M_C_DOSAGE;SHARED_INDIVIDUALS;COORDINATE_ANCHOR_NOT_VARIANT_IDENTITY'
POLICIES=('ALL_LOCAL_PAIRS','EXCLUDE_SAME_SOURCE_TAG_FAMILY')

def tag_family(value):
    m=re.fullmatch(r'(.*)--([0-9]+)',str(value or ''))
    return m.group(1) if m else None

def prepare_shared(ctx):
    cfg=ctx.baseline_config['dosage_ld']; names=list(cfg['resources']); originals=ctx.genotype_panels
    shared=sorted(set.intersection(*[{p.samples[int(i)] for i in p.groups['all_qc_samples']} for p in (originals[name] for name in names)]))
    flags=set(cfg['sensitivity_sample_ids'])
    groups={'shared_all':np.arange(len(shared),dtype=np.int64),'shared_omit_flagged':np.array([i for i,s in enumerate(shared) if s not in flags],dtype=np.int64)}
    if min(map(len,groups.values()))<cfg['qc']['min_pair_samples']:raise ValueError('Insufficient shared samples')
    panels,audit,sample_rows={},[],[]
    for name in names:
        original=originals[name]; lookup={s:i for i,s in enumerate(original.samples)}
        source_indices=np.array([lookup[s] for s in shared],dtype=np.int64); dosage=original.dosage[:,source_indices]
        metrics={g:row_qc(dosage[:,ii],cfg['qc']) for g,ii in groups.items()}
        keep=np.logical_and.reduce([m['keep'] for m in metrics.values()])
        for i,m in enumerate(original.markers):
            r=dict(resource=name,marker_key=m['marker_key'],chromosome=str(m['chromosome']),position_bp=int(m['position_bp']),retained_shared_sample_QC=bool(keep[i]),source_fixed_panel_markers=len(original.markers))
            for g,q in metrics.items():
                r[g+'_call_rate']=float(q['call_rate'][i]);r[g+'_maf']=float(q['maf'][i]) if np.isfinite(q['maf'][i]) else None;r[g+'_pass']=bool(q['keep'][i])
            audit.append(r)
        markers=[m for m,yes in zip(original.markers,keep) if yes]
        if len(markers)<cfg['qc']['min_markers']:raise ValueError('Insufficient shared-QC markers: '+name)
        coordinates=[(str(m['chromosome']),int(m['position_bp'])) for m in markers];keys=[m['marker_key'] for m in markers]
        if len(set(coordinates))!=len(coordinates) or len(set(keys))!=len(keys):raise ValueError('Ambiguous cross-export anchor: '+name)
        panels[name]=Panel(name,shared,markers,dosage[keep].copy(),groups,dict(source_fixed_markers=len(original.markers),shared_QC_markers=len(markers),shared_samples=len(shared)))
        for i,s in enumerate(shared):sample_rows.append(dict(resource=name,sample_id=s,shared_column_index=i,original_column_index=int(source_indices[i]),omitted_in_sensitivity=s in flags))
    ctx.write_csv('cross_export_shared_marker_QC.csv',audit);ctx.write_csv('cross_export_shared_sample_order.csv',sample_rows)
    ctx.check('cross_export_shared_samples',len(shared),152);ctx.check('cross_export_omission_samples',len(groups['shared_omit_flagged']),150)
    ctx.cross_panels=panels
    return panels

def physical_coverage(ledger):
    positions,exact=defaultdict(list),Counter()
    for r in ledger:
        if r['in_C'] and r['resource'] in ('Diversity','Nacional_G7','Nacional_MAF5'):
            key=(r['resource'],str(r['reference_chromosome']));bp=int(r['reference_snp_bp']);positions[key].append(bp);exact[key+(bp,)]+=1
    return {k:np.array(sorted(set(v)),dtype=np.int64) for k,v in positions.items()},exact

def geometry(peaks,panels,ledger):
    names=list(panels)
    if len(names)!=2:raise ValueError('Exactly two Nacional exports are required')
    anchors={name:{(str(m['chromosome']),int(m['position_bp'])):i for i,m in enumerate(p.markers)} for name,p in panels.items()}
    by_chrom={}
    for name,p in panels.items():
        by=defaultdict(list)
        for i,m in enumerate(p.markers):by[str(m['chromosome'])].append(i)
        by_chrom[name]={c:sorted(ii,key=lambda i:(int(p.markers[i]['position_bp']),p.markers[i]['marker_key'])) for c,ii in by.items()}
    positions,exact=physical_coverage(ledger);records,pairs=[],[]
    for peak in peaks:
        chrom,bp=str(peak['chromosome']),int(peak['start_bp'])
        if int(peak['end_bp'])!=bp:raise ValueError('Expected point peak')
        for target in names+['Diversity']:
            coord=(chrom,bp);source=None
            if target in panels:
                other=next(n for n in names if n!=target)
                source=target if coord in anchors[target] else other if coord in anchors[other] else None
            index=anchors[source][coord] if source else None;marker=panels[source].markers[index] if source else None
            pp=positions.get((target,chrom),np.empty(0,dtype=np.int64));i=int(np.searchsorted(pp,bp))
            nearest=min([abs(int(pp[j])-bp) for j in (i-1,i) if 0<=j<len(pp)],default=None)
            count50=int(np.searchsorted(pp,bp+50000,side='right')-np.searchsorted(pp,bp-50000,side='left'));nexact=exact[(target,chrom,bp)]
            if target=='Diversity':selection='NO_ACCESSION_GENOTYPES'
            elif source==target:selection='TARGET_PANEL_QC_ANCHOR'
            elif source:selection='OTHER_EXPORT_QC_ANCHOR_TARGET_PHYSICALLY_ABSENT' if not nexact else 'OTHER_EXPORT_QC_ANCHOR_TARGET_QC_EXCLUDED'
            else:selection='NO_QC_ANCHOR_IN_EITHER_NACIONAL_EXPORT'
            record=dict(query_id=peak['query_id'],chromosome=chrom,position_bp=bp,n_events=peak['n_events'],traits=peak['traits'],target_resource=target,
                target_exact_C_records=nexact,target_physical_status='EXACT_COORDINATE' if nexact else 'NEARBY_ONLY' if count50 else 'NO_MARKER_WITHIN_50KB',
                nearest_target_C_distance_bp=nearest,target_unique_C_positions_within_50kb=count50,target_shared_QC_anchor=coord in anchors.get(target,{}),
                coordinate_anchor_source=source or '',coordinate_anchor_marker_key=marker['marker_key'] if marker else '',anchor_selection=selection,
                original_variant_identity='NOT_ESTABLISHED_BY_COORDINATE',analysis_scope=SCOPE)
            ri=len(records);records.append(record)
            if not source:continue
            for j in by_chrom[target].get(chrom,[]):
                partner=panels[target].markers[j];distance=abs(int(partner['position_bp'])-bp)
                if not 0<distance<=50000 or partner['marker_key']==marker['marker_key']:continue
                a,b=tag_family(marker['source_tag_key']),tag_family(partner['source_tag_key'])
                if a is None or b is None:raise ValueError('Unparseable source tag family')
                pairs.append(dict(record_index=ri,query_id=peak['query_id'],chromosome=chrom,source_resource=source,target_resource=target,
                    source_index=index,target_index=j,coordinate_anchor_marker_key=marker['marker_key'],partner_marker_key=partner['marker_key'],
                    anchor_position_bp=bp,partner_position_bp=int(partner['position_bp']),distance_bp=distance,same_source_tag_family=a==b,cross_export_pair=source!=target))
    return records,pairs

def classification(record,eligible,valid,best):
    if record['target_resource']=='Diversity':return 'NOT_ASSESSABLE_NO_ACCESSION_GENOTYPES'
    if not record['coordinate_anchor_source']:return 'NOT_ASSESSABLE_NO_QC_ANCHOR_IN_EITHER_EXPORT'
    if not eligible:return 'NO_ELIGIBLE_TARGET_NEIGHBOR_WITHIN_50KB'
    if not valid:return 'NOT_ASSESSABLE_NO_VALID_CORRELATION'
    high=best>=.8
    if record['target_shared_QC_anchor']:return 'TARGET_EXACT_QC_ANCHOR_HIGH_NEIGHBOR' if high else 'TARGET_EXACT_QC_ANCHOR_NEIGHBORS_BELOW_R2_0_8'
    if record['target_exact_C_records']:return 'TARGET_EXACT_QC_EXCLUDED_DONOR_ANCHOR_HIGH_NEIGHBOR' if high else 'TARGET_EXACT_QC_EXCLUDED_DONOR_NEIGHBORS_BELOW_R2_0_8'
    return 'TARGET_ABSENT_DONOR_COORDINATE_ANCHOR_HIGH_NEIGHBOR' if high else 'TARGET_ABSENT_DONOR_NEIGHBORS_BELOW_R2_0_8'

def raw_pair_correlations(xrows,yrows,minimum=50):
    """Direct integer-sum Pearson formula for complete calls."""
    x,y=np.asarray(xrows,dtype=float),np.asarray(yrows,dtype=float)
    if x.ndim!=2 or x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Complete, aligned pair matrices required')
    if not np.isin(x,[0,1,2]).all() or not np.isin(y,[0,1,2]).all():raise ValueError('Invalid dosage')
    n=x.shape[1];sx=x.sum(axis=1);sy=y.sum(axis=1)
    cov=n*(x*y).sum(axis=1)-sx*sy;vx=n*(x*x).sum(axis=1)-sx*sx;vy=n*(y*y).sum(axis=1)-sy*sy
    valid=(n>=minimum)&(vx>0)&(vy>0);r2=np.full(len(x),np.nan)
    r2[valid]=cov[valid]**2/(vx[valid]*vy[valid])
    if np.any(r2[valid]>1+1e-12):raise ValueError('Raw correlation out of range')
    return dict(r2=np.clip(r2,0,1),n=np.full(len(x),n,np.int32),status=np.where(n<minimum,1,np.where(valid,0,2)).astype(np.uint8))

def best_partner(valid,pairs,r2,tolerance=1e-12):
    if not valid:return None
    maximum=max(float(r2[i]) for i in valid)
    tied=[i for i in valid if maximum-float(r2[i])<=tolerance]
    return min(tied,key=lambda i:(pairs[i]['partner_position_bp'],pairs[i]['partner_marker_key']))

def run(ctx):
    from adjusted_math import pair_correlations
    from model_setup import fit_models, METHODS
    cfg=ctx.baseline_config['dosage_ld'];panels=prepare_shared(ctx)
    peaks=ctx.read('source_data/02_association_units.xlsx','unique_peaks')
    ledger=ctx.read('alignment_coverage/02_fixed_marker_cohorts.xlsx','marker_record_ledger')
    records,pairs=geometry(peaks,panels,ledger)
    ctx.check('cross_peak_unique_queries',len(peaks),1415);ctx.check('cross_peak_three_panel_inventory',len(records),3*len(peaks))
    ctx.write_csv('cross_peak_inventory.csv',records)
    pair_by_chrom=defaultdict(list);pair_by_record=defaultdict(list)
    for i,p in enumerate(pairs):pair_by_chrom[p['chromosome']].append(i);pair_by_record[p['record_index']].append(i)
    first=next(iter(panels.values()));groups=first.groups;shared=first.samples
    for p in panels.values():
        if p.samples!=shared:raise ValueError('Cross-panel sample order changed')
    common=ctx.common_genotype_panels
    basis=[[(m['marker_key'],str(m['chromosome']),int(m['position_bp']),m['source_tag_key']) for m in p.markers] for p in common.values()]
    ctx.check('cross_model_common_marker_basis_matches',all(x==basis[0] for x in basis),True);ctx.check('cross_model_common_marker_count',len(basis[0]),4285)
    backgrounds={}
    for source,p in common.items():
        by_sample={s:i for i,s in enumerate(p.samples)};backgrounds[source]=p.dosage[:,[by_sample[s] for s in shared]]
    group_outputs={};pair_out=[];assessment=[];summary=Counter();training=[];background_audit=[];projection_audit=[];eigen_audit=[]
    for covariance_source,background in backgrounds.items():
        for group,indices in groups.items():
            print('[CROSS PEAK]',covariance_source,group,'LOCO models and coordinate anchors',flush=True)
            values={method:dict(r2=np.full(len(pairs),np.nan),n=np.zeros(len(pairs),np.int32),status=np.ones(len(pairs),np.uint8)) for method in METHODS}
            covered=np.zeros(len(pairs),np.int32)
            for chromosome in cfg['chromosomes']:
                fitted=fit_models(background[:,indices],common[covariance_source].markers,chromosome)
                selected=fitted['background_indices']
                leakage=sum(str(common[covariance_source].markers[int(i)]['chromosome'])==str(chromosome) for i in selected)
                ctx.check('cross_LOCO_leakage:'+covariance_source+':'+group+':'+str(chromosome),leakage,0)
                if leakage:raise ValueError('Held-out chromosome entered covariance background')
                context=dict(covariance_source=covariance_source,sample_group=group,held_out_chromosome=str(chromosome))
                background_audit.append(dict(covariance_source=covariance_source,sample_group=group,**fitted['background_audit']))
                for rank,index in enumerate(fitted['background_indices'],1):
                    m=common[covariance_source].markers[int(index)];training.append(dict(**context,selection_order=rank,marker_key=m['marker_key'],chromosome=str(m['chromosome']),position_bp=int(m['position_bp']),source_tag_key=m['source_tag_key']))
                for rank,value in enumerate(fitted['eigenvalues'],1):eigen_audit.append(dict(**context,eigen_rank=rank,eigenvalue=float(value)))
                ii=np.array(pair_by_chrom.get(str(chromosome),[]),dtype=np.int64)
                covered[ii]+=1
                if len(ii):
                    x=np.stack([panels[pairs[int(i)]['source_resource']].dosage[pairs[int(i)]['source_index'],indices] for i in ii])
                    y=np.stack([panels[pairs[int(i)]['target_resource']].dosage[pairs[int(i)]['target_index'],indices] for i in ii])
                for method,model in fitted['models'].items():
                    projection_audit.append(dict(**context,method=method,**model['diagnostics']))
                    if len(ii):
                        v=pair_correlations(x,y,model['T'],cfg['qc']['min_pair_samples'])
                        if method=='RAW':
                            direct=raw_pair_correlations(x,y,cfg['qc']['min_pair_samples'])
                            ctx.check('cross_projected_raw:'+covariance_source+':'+group+':'+str(chromosome),bool(np.allclose(v['r2'],direct['r2'],rtol=1e-10,atol=1e-12,equal_nan=True)),True)
                            v=direct
                        for key in values[method]:values[method][key][ii]=v[key]
                if len(ii):
                    jj=np.linspace(0,len(ii)-1,min(7,len(ii)),dtype=int);expected=np.array([np.corrcoef(x[j],y[j])[0,1]**2 for j in jj])
                    ctx.check('cross_raw_Pearson:'+covariance_source+':'+group+':'+str(chromosome),bool(np.allclose(values['RAW']['r2'][ii[jj]],expected,rtol=1e-10,atol=1e-12,equal_nan=True)),True)
            ctx.check('cross_each_pair_evaluated_once:'+covariance_source+':'+group,bool(np.all(covered==1)),True)
            for method,v in values.items():
                ctx.check('cross_pair_sample_count:'+covariance_source+':'+group+':'+method,bool(np.all(v['n']==len(indices))),True)
            group_outputs[(covariance_source,group)]=values
            for i,p in enumerate(pairs):
                row={k:v for k,v in p.items() if k not in ('record_index','source_index','target_index')};row.update(covariance_source=covariance_source,sample_group=group,analysis_scope=SCOPE)
                for method,v in values.items():
                    row[method+'__r2']=float(v['r2'][i]) if np.isfinite(v['r2'][i]) else None;row[method+'__n']=int(v['n'][i]);row[method+'__status_code']=int(v['status'][i])
                pair_out.append(row)
            for ri,record in enumerate(records):
                if record['target_resource']=='Diversity':continue
                all_ids=pair_by_record.get(ri,[])
                for policy in POLICIES:
                    ii=[i for i in all_ids if policy=='ALL_LOCAL_PAIRS' or not pairs[i]['same_source_tag_family']]
                    row=dict(record,covariance_source=covariance_source,sample_group=group,pair_policy=policy,shared_sample_count=len(indices),eligible_target_neighbors=len(ii))
                    for method,v in values.items():
                        valid=[i for i in ii if np.isfinite(v['r2'][i])]
                        best=best_partner(valid,pairs,v['r2'])
                        value=float(v['r2'][best]) if best is not None else None;label=classification(record,len(ii),len(valid),value)
                        row[method+'__valid_pairs']=len(valid);row[method+'__best_r2']=value;row[method+'__best_partner_key']=pairs[best]['partner_marker_key'] if best is not None else ''
                        row[method+'__best_partner_distance_bp']=pairs[best]['distance_bp'] if best is not None else None;row[method+'__best_pair_samples']=int(v['n'][best]) if best is not None else None;row[method+'__classification']=label
                        summary[(record['target_resource'],covariance_source,group,policy,method,label)]+=1
                    assessment.append(row)
    sources=list(backgrounds)
    if len(sources)!=2:raise ValueError('Two covariance-source sensitivities required')
    for group in groups:
        a,b=(group_outputs[(source,group)] for source in sources)
        ctx.check('cross_raw_covariance_source_invariance:'+group,bool(np.allclose(a['RAW']['r2'],b['RAW']['r2'],rtol=0,atol=1e-12,equal_nan=True)),True)
        if group=='shared_omit_flagged':
            for method in METHODS:ctx.check('cross_omitted_covariance_source_agreement:'+method,bool(np.allclose(a[method]['r2'],b[method]['r2'],rtol=1e-8,atol=1e-10,equal_nan=True)),True)
    for record in records:
        if record['target_resource']!='Diversity':continue
        row=dict(record,covariance_source='NOT_APPLICABLE',sample_group='NOT_APPLICABLE',pair_policy='NOT_APPLICABLE',shared_sample_count=None,eligible_target_neighbors=None)
        for method in METHODS:
            for suffix in ('valid_pairs','best_r2','best_partner_distance_bp','best_pair_samples'):row[method+'__'+suffix]=None
            row[method+'__best_partner_key']='';row[method+'__classification']='NOT_ASSESSABLE_NO_ACCESSION_GENOTYPES'
        assessment.append(row)
    summaries=[dict(target_resource=k[0],covariance_source=k[1],sample_group=k[2],pair_policy=k[3],method=k[4],classification=k[5],peak_count=count,peak_denominator=len(peaks),analysis_scope=SCOPE) for k,count in sorted(summary.items())]
    ctx.check('cross_peak_wide_assessment_rows',len(assessment),len(peaks)*(2*len(groups)*len(sources)*len(POLICIES)+1))
    denom=Counter()
    for key,count in summary.items():denom[key[:5]]+=count
    ctx.check('cross_peak_every_method_denominator',set(denom.values())=={len(peaks)},True)
    noanchor=[r for r in assessment if not r['coordinate_anchor_source']]
    ctx.check('cross_missing_anchor_correlations_blank',all(r[m+'__best_r2'] is None for r in noanchor for m in METHODS),True)
    for name,rows in [('cross_peak_adjusted_assessment.csv',assessment),('cross_peak_adjusted_summary.csv',summaries),('cross_peak_pair_correlations.csv',pair_out),('cross_model_training_markers.csv',training),('cross_model_background_audit.csv',background_audit),('cross_model_projection_audit.csv',projection_audit),('cross_model_eigenvalues.csv',eigen_audit)]:ctx.write_csv(name,rows)
    counts=Counter((r['target_resource'],r['anchor_selection']) for r in records)
    ctx.summary['cross_export_peak_audit']=dict(analysis_scope=SCOPE,shared_samples=len(shared),omitted_samples=len(groups['shared_omit_flagged']),shared_QC_markers={n:len(p.markers) for n,p in panels.items()},common_covariance_marker_basis=len(basis[0]),directed_anchor_neighbor_pairs=len(pairs),wide_pair_rows=len(pair_out),wide_peak_rows=len(assessment),source_selection_counts=[dict(target_resource=k[0],selection=k[1],peaks=v) for k,v in sorted(counts.items())],covariance_sources=list(backgrounds),methods=list(METHODS),source_genotypes_merged=False,identity_status='POSITIONAL_ANCHOR_ONLY;NOT_VALIDATED_DISCOVERY_VARIANT_OR_BREEDING_PROXY')
    ctx.note('Cross-export analysis aligns shared individuals and rechecks fixed-marker QC in both shared groups. Target-panel QC anchors take precedence; only unavailable target anchors are borrowed from the other export. No exported cells are merged or repaired.')
    ctx.note('Cross-export RAW values are repeated under covariance-source labels only to align comparisons, not as independent estimates. G7/MAF5 covariance backgrounds use the same marker IDs and sample order, preserving both discordant samples in the inclusive analysis.')
    ctx.note('Best partners within an absolute r2 difference of 1e-12 are tied and selected by target position then canonical key. RAW values use the direct integer-sum Pearson formula after checking agreement with the projection implementation.')
    ctx.note('Donor-export coordinate-anchor correlation does not establish original peak-variant identity. Target physical absence, target QC exclusion, no anchor in either export, and non-estimable correlation remain distinct.')
