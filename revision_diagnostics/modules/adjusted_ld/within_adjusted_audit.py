"""Within-panel conditional adjusted dosage r2 on fixed quality-controlled cohorts.

Models are estimated using chromosomes other than the tested chromosome.
RAW uses complete-case arithmetic after a separate projection-equivalence check.
Alleles, coordinates and QC cohorts are held fixed.
"""
from collections import defaultdict
import numpy as np
from adjusted_math import metric_pairvalues
from model_setup import fit_models,METHODS
from ld_audit import local_pairs,correlation,statistics,compare,equal,family,scalar

SCOPE='CONDITIONAL_A_M_C_DOSAGE;LOCO_SENSITIVITY;OBSERVED_MARKERS_ONLY'
POLICIES=('ALL_LOCAL_PAIRS','EXCLUDE_SAME_SOURCE_TAG_FAMILY','EXCLUDE_IDENTICAL_ALIGNED_TAG_SPAN')


def neighbor_rows(panel,pairs,dist,values,include,window=50000,tie_tolerance=0.0):
    n=len(panel.markers);local=np.flatnonzero((dist>0)&(dist<=window)&include)
    eligible=np.bincount(pairs[local].ravel(),minlength=n)
    valid=local[np.isfinite(values['r2'][local])]
    evaluable=np.bincount(pairs[valid].ravel(),minlength=n)
    best=np.full(n,np.nan);partner=np.full(n,-1,int);distance=np.full(n,-1,int)
    for z in valid:
        i,j=map(int,pairs[z]);r=float(values['r2'][z]);d=int(dist[z])
        for k,other in ((i,j),(j,i)):
            better=partner[k]<0 or r>best[k]+tie_tolerance
            tied=partner[k]>=0 and abs(r-best[k])<=tie_tolerance and other<partner[k]
            if better or tied:best[k],partner[k],distance[k]=r,other,d
    rows=[]
    for i,m in enumerate(panel.markers):
        j=int(partner[i])
        rows.append(dict(marker_key=m['marker_key'],chromosome=str(m['chromosome']),position_bp=int(m['position_bp']),window_bp=window,
            eligible_other_markers=int(eligible[i]),valid_other_markers=int(evaluable[i]),
            best_partner_key=panel.markers[j]['marker_key'] if j>=0 else '',
            best_partner_distance_bp=int(distance[i]) if j>=0 else None,best_r2=scalar(best[i]),
            status='EVALUABLE' if j>=0 else 'NO_VALID_PAIR' if eligible[i] else 'NO_OTHER_MARKER_IN_WINDOW'))
    summary=dict(window_bp=window,typed_markers=n,markers_with_other_marker=int(np.sum(eligible>0)),
                 markers_with_valid_pair=int(np.sum(evaluable>0)),markers_without_any_neighbor=int(np.sum(eligible==0)),
                 markers_with_neighbors_but_no_valid_pair=int(np.sum((eligible>0)&(evaluable==0))))
    for threshold in (.2,.5,.8):
        key=str(threshold).replace('.','_');nt=int(np.sum(best>=threshold));ne=int(np.sum(evaluable>0))
        summary['markers_best_r2_ge_'+key]=nt
        summary['fraction_all_markers_ge_'+key]=nt/n if n else None
        summary['fraction_evaluable_ge_'+key]=nt/ne if ne else None
    return rows,summary,best


def comparison_statistics(raw,adjusted,mask):
    ok=mask&np.isfinite(raw)&np.isfinite(adjusted);a=raw[ok];b=adjusted[ok]
    return dict(inventory_count=int(mask.sum()),valid_both=int(ok.sum()),
        raw_only_evaluable=int(np.sum(mask&np.isfinite(raw)&~np.isfinite(adjusted))),
        adjusted_only_evaluable=int(np.sum(mask&~np.isfinite(raw)&np.isfinite(adjusted))),
        high_raw=int(np.sum(a>=.8)),high_adjusted=int(np.sum(b>=.8)),
        retained_high=int(np.sum((a>=.8)&(b>=.8))),lost_high=int(np.sum((a>=.8)&(b<.8))),
        gained_high=int(np.sum((a<.8)&(b>=.8))),
        mean_raw_fixed_evaluable=float(a.mean()) if len(a) else None,
        mean_adjusted_fixed_evaluable=float(b.mean()) if len(b) else None,
        mean_delta_adjusted_minus_raw=float((b-a).mean()) if len(a) else None,
        median_absolute_delta=float(np.median(np.abs(b-a))) if len(a) else None)


def _tag_flags(panel,pairs,span_index):
    families=[family(m['source_tag_key']) for m in panel.markers]
    spans=[next(iter(span_index[str(m['source_tag_key'])])) if len(span_index[str(m['source_tag_key'])])==1 else None for m in panel.markers]
    if any(f is None for f in families) or any(s is None for s in spans):
        raise ValueError('Every fixed marker must have a source tag family and unique aligned span')
    same_family=np.array([families[i]==families[j] for i,j in pairs],bool)
    same_span=np.array([spans[i]==spans[j] for i,j in pairs],bool)
    return same_family,same_span


def run(ctx):
    cfg=ctx.baseline_config['dosage_ld'];lc=cfg['ld'];thresholds=lc['r2_thresholds']
    frozen_bins=ctx.read('dosage_ld/02_ld_distance_profiles.xlsx','distance_bins')
    frozen_sum=ctx.read('dosage_ld/05_observed_snp_neighbors.xlsx','summary')
    frozen_neighbor=ctx.read('dosage_ld/05_observed_snp_neighbors.xlsx','marker_neighbors')
    frozen_inv=ctx.read('dosage_ld/02_ld_distance_profiles.xlsx','pair_inventories')
    fb={(r['resource'],r['sample_group'],str(r['chromosome']),r['lower_bp_exclusive'],r['upper_bp_inclusive']):r for r in frozen_bins if r['coordinate_basis']=='corrected_snp'}
    fs={(r['resource'],r['sample_group']):r for r in frozen_sum if r['window_bp']==50000}
    fn={(r['resource'],r['sample_group'],r['marker_key']):r for r in frozen_neighbor if r['window_bp']==50000}
    fi={r['resource']:r for r in frozen_inv}
    decisions=ctx.read('alignment_coverage/01_realignment_evidence.xlsx','marker_decisions');spans=defaultdict(set)
    for r in decisions:
        if r['accepted_for_mapping']:
            spans[str(r['reference_key'])].add((str(r['mapped_chromosome']),int(r['tag_start_bp']),int(r['tag_end_bp']),str(r['mapped_strand'])))
    profiles=[];all_neighbors=[];summaries=[];transition_rows=[];pair_records=[];differences=[]
    background_rows=[];training_rows=[];model_rows=[];eigen_rows=[];inventory=[]
    for resource,panel in ctx.genotype_panels.items():
        print('[ADJUSTED WITHIN]',resource,'fixed markers',len(panel.markers),flush=True)
        pairs=local_pairs(panel.markers,int(lc['max_distance_bp']),int(lc['max_union_pairs']))
        chrom=np.array([str(m['chromosome']) for m in panel.markers]);pos=np.array([int(m['position_bp']) for m in panel.markers])
        pair_chrom=chrom[pairs[:,0]];dist=np.abs(pos[pairs[:,1]]-pos[pairs[:,0]])
        ctx.check('within corrected pair census:'+resource,len(pairs),int(fi[resource]['corrected_local_pairs']))
        same_family,same_span=_tag_flags(panel,pairs,spans)
        masks={POLICIES[0]:np.ones(len(pairs),bool),POLICIES[1]:~same_family,POLICIES[2]:~same_span}
        for group,sample_ids in panel.groups.items():
            d=np.asarray(panel.dosage[:,sample_ids],float)
            if not np.isfinite(d).all():raise ValueError('Within-panel calls must be complete; no imputation')
            values={method:dict(r2=np.full(len(pairs),np.nan),n=np.zeros(len(pairs),np.int32),status=np.ones(len(pairs),np.uint8)) for method in METHODS}
            covered=np.zeros(len(pairs),np.int8)
            for chromosome in cfg['chromosomes']:
                ids=np.flatnonzero(chrom==str(chromosome));which=np.flatnonzero(pair_chrom==str(chromosome))
                fit=fit_models(d,panel.markers,str(chromosome))
                background_rows.append(dict(resource=resource,sample_group=group,**fit['background_audit'],analysis_scope=SCOPE))
                selected=fit['background_indices']
                ctx.check('within LOCO leakage:'+resource+':'+group+':'+str(chromosome),sum(str(panel.markers[i]['chromosome'])==str(chromosome) for i in selected),0)
                for order,i in enumerate(selected):
                    m=panel.markers[i]
                    training_rows.append(dict(resource=resource,sample_group=group,held_out_chromosome=str(chromosome),selection_order=order+1,
                        marker_key=m['marker_key'],chromosome=str(m['chromosome']),position_bp=int(m['position_bp']),source_tag_key=m['source_tag_key']))
                for order,eig in enumerate(fit['eigenvalues']):
                    eigen_rows.append(dict(resource=resource,sample_group=group,held_out_chromosome=str(chromosome),component=order+1,eigenvalue=float(eig)))
                mapper=np.full(len(panel.markers),-1,int);mapper[ids]=np.arange(len(ids));local=mapper[pairs[which]]
                for method,model in fit['models'].items():
                    model_rows.append(dict(resource=resource,sample_group=group,held_out_chromosome=str(chromosome),method=method,**model['diagnostics']))
                    vals=metric_pairvalues(d[ids],local,model['T'],cfg['qc']['min_pair_samples'])
                    for key in vals:values[method][key][which]=vals[key]
                covered[which]+=1
            ctx.check('within every local pair evaluated once:'+resource+':'+group,bool(np.all(covered==1)),True)
            raw=correlation(d,pairs,cfg['qc']['min_pair_samples'],lc['batch_pairs'])
            valid=np.isfinite(raw['r2'])&np.isfinite(values['RAW']['r2'])
            ctx.check('within RAW projection statuses:'+resource+':'+group,bool(np.array_equal(raw['status'],values['RAW']['status'])),True)
            error=float(np.max(np.abs(raw['r2'][valid]-values['RAW']['r2'][valid]))) if valid.any() else 0.
            ctx.check('within RAW projection equivalence:'+resource+':'+group,error<1e-10,True)
            values['RAW']=raw
            inventory.append(dict(resource=resource,sample_group=group,qc_markers=len(panel.markers),samples=len(sample_ids),
                local_pairs=len(pairs),pairs_within_50kb=int(np.sum((dist>0)&(dist<=50000))),raw_projection_max_abs_error=error))
            rawbest={}
            for method in METHODS:
                v=values[method]
                for chromosome in ['ALL']+cfg['chromosomes']:
                    cmask=np.ones(len(pairs),bool) if chromosome=='ALL' else pair_chrom==str(chromosome)
                    for lo,hi in [(None,0)]+list(zip(lc['distance_edges_bp'][:-1],lc['distance_edges_bp'][1:])):
                        base=cmask&(dist==0 if lo is None else (dist>lo)&(dist<=hi))
                        for policy,include in masks.items():
                            stats=statistics(v,base&include)
                            profiles.append(dict(resource=resource,sample_group=group,method=method,chromosome=str(chromosome),
                                lower_bp_exclusive=lo,upper_bp_inclusive=hi,pair_policy=policy,**stats,analysis_scope=SCOPE))
                            if method=='RAW' and policy=='ALL_LOCAL_PAIRS':
                                key=(resource,group,str(chromosome),lo,hi)
                                differences.extend(compare(ctx,'within frozen bin:'+':'.join(map(str,key)),stats,fb[key],stats))
                for policy,include in masks.items():
                    rows,summary,best=neighbor_rows(panel,pairs,dist,v,include,tie_tolerance=0 if method=='RAW' else 1e-12)
                    fields=list(summary)
                    summaries.append(dict(resource=resource,sample_group=group,method=method,pair_policy=policy,**summary,analysis_scope=SCOPE))
                    for row in rows:all_neighbors.append(dict(resource=resource,sample_group=group,method=method,pair_policy=policy,**row,analysis_scope=SCOPE))
                    if method=='RAW':rawbest[policy]=best
                    else:
                        transition_rows.append(dict(resource=resource,sample_group=group,method=method,pair_policy=policy,
                            comparison_unit='MARKER_BEST_NEIGHBOR_50KB',**comparison_statistics(rawbest[policy],best,np.ones(len(best),bool)),analysis_scope=SCOPE))
                        transition_rows.append(dict(resource=resource,sample_group=group,method=method,pair_policy=policy,
                            comparison_unit='FIXED_MARKER_PAIRS_50KB',**comparison_statistics(raw['r2'],v['r2'],include&(dist>0)&(dist<=50000)),analysis_scope=SCOPE))
                    if method=='RAW' and policy=='ALL_LOCAL_PAIRS':
                        differences.extend(compare(ctx,'within frozen neighbor summary:'+resource+':'+group,summary,fs[(resource,group)],fields))
                        bad=[]
                        for row in rows:
                            ref=fn[(resource,group,row['marker_key'])]
                            for field in ('eligible_other_markers','valid_other_markers','best_partner_key','best_partner_distance_bp','best_r2','status'):
                                if not equal(row[field],ref[field]):bad.append(dict(comparison='within frozen marker:'+resource+':'+group,marker_key=row['marker_key'],field=field,observed=row[field],frozen=ref[field]))
                        ctx.check('within frozen neighbor rows:'+resource+':'+group,len(bad),0);differences.extend(bad)
            for z in np.flatnonzero((dist>0)&(dist<=50000)):
                i,j=map(int,pairs[z]);a,b=panel.markers[i],panel.markers[j]
                row=dict(resource=resource,sample_group=group,marker_a=a['marker_key'],marker_b=b['marker_key'],chromosome=str(a['chromosome']),
                    position_a=int(a['position_bp']),position_b=int(b['position_bp']),distance_bp=int(dist[z]),
                    same_source_tag_family=bool(same_family[z]),same_aligned_tag_span=bool(same_span[z]),analysis_scope=SCOPE)
                for method in METHODS:
                    row[method+'_r2']=scalar(values[method]['r2'][z]);row[method+'_n']=int(values[method]['n'][z]);row[method+'_status']=int(values[method]['status'][z])
                pair_records.append(row)
    for name,rows in [('within_model_background_audit.csv',background_rows),('within_model_training_markers.csv',training_rows),
                      ('within_model_diagnostics.csv',model_rows),('within_pca_eigenvalues.csv',eigen_rows),
                      ('within_pair_inventory.csv',inventory),('within_adjusted_distance_profiles.csv',profiles),
                      ('within_adjusted_neighbors_50kb.csv',all_neighbors),('within_adjusted_neighbor_summary.csv',summaries),
                      ('within_adjustment_transitions_50kb.csv',transition_rows),('within_adjusted_pairs_50kb.csv',pair_records)]:
        ctx.write_csv(name,rows)
    ctx.write_csv('within_frozen_RAW_comparison_differences.csv',differences,['comparison','marker_key','field','observed','frozen'])
    ctx.summary['within_adjusted_audit']=dict(scope=SCOPE,methods=list(METHODS),pair_policies=list(POLICIES),
        panel_inventories=inventory,training_fits=len(background_rows),frozen_RAW_difference_cells=len(differences),
        primary_neighbor_summaries=[r for r in summaries if r['sample_group']=='all_qc_samples' and r['pair_policy']=='ALL_LOCAL_PAIRS'],
        interpretation='PC and regularized covariance sensitivity estimates; phenotype effects and biological coding remain unvalidated.')
    ctx.note('Within-panel adjusted LD uses fixed QC markers and inclusive/omission sample groups. Each chromosome uses background markers only from other chromosomes.')
    ctx.note('Kinship metrics use mean-diagonal-one dosage covariance plus explicit ridge terms. PC2/PC5 are residual correlations after adjustment for two or five genetic PCs.')
    ctx.note('Adjusted best-neighbor ties within 1e-12 select the lowest fixed marker index. Raw neighbor correlations use exact complete-case arithmetic.')
