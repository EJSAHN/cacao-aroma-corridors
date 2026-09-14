"""Distance profiles, sample sensitivities and observed-SNP neighbor diagnostics."""
from __future__ import annotations
import numpy as np
from .common import DataError
from .genotypes import Panel, prepare_common
from .numerics import pair_inventory,pairwise_ld,summary_statistics


def coords(panel:Panel):
    return tuple(np.array([r[key] for r in panel.markers],dtype=dtype)
                 for key,dtype in [('original_chromosome',str),('original_position_bp',np.int64),('chromosome',str),('position_bp',np.int64)])


def geometry(pairs,chroms,pos):
    i,j=pairs[:,0],pairs[:,1]
    return chroms[i]==chroms[j], np.abs(pos[j]-pos[i])


def panel_analysis(panel:Panel,cfg:dict,logger):
    oldc,oldp,newc,newp=coords(panel)
    pairs,flags,pop=pair_inventory(oldc,oldp,newc,newp,cfg['ld'],cfg['seed'],panel.name)
    logger.info('LD %s: %d markers, %d fixed pairs; samples %s',panel.name,len(panel.markers),len(pairs),{k:len(v) for k,v in panel.groups.items()})
    values={g:pairwise_ld(panel.dosage[:,idx],pairs,cfg['qc']['min_pair_samples'],cfg['ld']['batch_pairs']) for g,idx in panel.groups.items()}
    geometries={'corrected_snp':geometry(pairs,newc,newp),'original_position':geometry(pairs,oldc,oldp)}
    profiles=[]; windows=[]; context=[]; shifts=[]
    thresholds=cfg['ld']['r2_thresholds']; edges=cfg['ld']['distance_edges_bp']
    for mode,(same,dist) in geometries.items():
        cc=newc if mode=='corrected_snp' else oldc
        for chrom in ['ALL']+cfg['chromosomes']:
            base=same.copy() if chrom=='ALL' else same&(cc[pairs[:,0]]==chrom)
            for lo,hi in [(0,0)]+list(zip(edges[:-1],edges[1:])):
                mask=base&(dist==0) if lo==hi else base&(dist>lo)&(dist<=hi)
                for group,v in values.items():
                    profiles.append(dict(resource=panel.name,sample_group=group,coordinate_basis=mode,chromosome=chrom,
                        lower_bp_exclusive=None if lo==hi else lo,upper_bp_inclusive=hi,
                        bin_definition='ZERO_DISTANCE_SEPARATE' if lo==hi else 'LOWER_EXCLUSIVE_UPPER_INCLUSIVE',
                        **summary_statistics(v,mask,thresholds)))
        for w in cfg['ld']['report_windows_bp']:
            mask=same&(dist<=w)
            for group,v in values.items():
                windows.append(dict(resource=panel.name,sample_group=group,coordinate_basis=mode,max_distance_bp=w,
                                    **summary_statistics(v,mask,thresholds)))
    for label in ['interchromosomal_context','distant_context']:
        for group,v in values.items():
            context.append(dict(resource=panel.name,sample_group=group,pair_class=label,
                                scope='UNIFORM_SAMPLE_OF_QC_MARKER_PAIRS; NOT_A_NULL_TEST',
                                **summary_statistics(v,flags[label],thresholds)))
    for w in cfg['ld']['report_windows_bp']:
        b=geometries['original_position'][0]&(geometries['original_position'][1]<=w)
        c=geometries['corrected_snp'][0]&(geometries['corrected_snp'][1]<=w)
        shifts.append(dict(resource=panel.name,max_distance_bp=w,both=int(np.sum(b&c)),
                           original_only=int(np.sum(b&~c)),corrected_only=int(np.sum(c&~b)),
                           note='Identical genotype pair values; coordinates change only pair-window membership.'))
    return dict(panel=panel,pairs=pairs,flags=flags,population=pop,values=values,geometries=geometries,
                profiles=profiles,windows=windows,context=context,coordinate_window_changes=shifts)


def sample_sensitivity(result,cfg):
    panel=result['panel'];values=result['values']
    a,b=values['all_qc_samples'],values['omit_flagged_samples']
    same,dist=result['geometries']['corrected_snp'];edges=cfg['ld']['distance_edges_bp']
    slices=[('distance_bin',lo,hi,same&(dist>lo)&(dist<=hi)) for lo,hi in zip(edges[:-1],edges[1:])]
    slices += [('cumulative_window',0,w,same&(dist<=w)) for w in cfg['ld']['report_windows_bp']]
    slices += [(label,None,None,result['flags'][label]) for label in ['interchromosomal_context','distant_context']]
    rows=[]
    for kind,lo,hi,mask in slices:
        valid=mask&np.isfinite(a.r2)&np.isfinite(b.r2);delta=b.r2[valid]-a.r2[valid]
        r=dict(resource=panel.name,pair_class=kind,lower_bp_exclusive=lo,upper_bp_inclusive=hi,
               inventory_pairs=int(mask.sum()),valid_both=int(valid.sum()),
               valid_all_only=int(np.sum(mask&np.isfinite(a.r2)&~np.isfinite(b.r2))),
               valid_omitted_only=int(np.sum(mask&~np.isfinite(a.r2)&np.isfinite(b.r2))),
               mean_r2_all_fixed_pairs=float(a.r2[valid].mean()) if valid.any() else None,
               mean_r2_omitted_fixed_pairs=float(b.r2[valid].mean()) if valid.any() else None,
               mean_delta_r2=float(delta.mean()) if len(delta) else None,
               median_abs_delta_r2=float(np.median(np.abs(delta))) if len(delta) else None,
               q95_abs_delta_r2=float(np.quantile(np.abs(delta),0.95)) if len(delta) else None,
               max_abs_delta_r2=float(np.max(np.abs(delta))) if len(delta) else None)
        for t in cfg['ld']['r2_thresholds']:
            k=str(t).replace('.','_')
            r['crossed_up_'+k]=int(np.sum((a.r2[valid]<t)&(b.r2[valid]>=t)))
            r['crossed_down_'+k]=int(np.sum((a.r2[valid]>=t)&(b.r2[valid]<t)))
        rows.append(r)
    return rows


def neighbor_diagnostics(result,cfg):
    """For each actually typed marker, evaluate other typed markers. No self tags."""
    panel=result['panel']; pairs=result['pairs'];same,dist=result['geometries']['corrected_snp']
    n=len(panel.markers);rows=[];summaries=[]
    for w in cfg['ld']['report_windows_bp']:
        local=np.flatnonzero(same&(dist>0)&(dist<=w))
        eligible=np.bincount(pairs[local].ravel(),minlength=n)
        for group,v in result['values'].items():
            valid=local[np.isfinite(v.r2[local])]
            nvalid=np.bincount(pairs[valid].ravel(),minlength=n)
            best=np.full(n,-np.inf);partner=np.full(n,-1,dtype=np.int64);bdist=np.zeros(n,dtype=np.int64)
            # Stable iteration over sorted index pairs makes ties deterministic.
            for z in valid:
                i,j=map(int,pairs[z]);r=float(v.r2[z]);d=int(dist[z])
                for k,other in ((i,j),(j,i)):
                    if r>best[k] or (r==best[k] and (partner[k]<0 or other<partner[k])):
                        best[k]=r;partner[k]=other;bdist[k]=d
            for i,m in enumerate(panel.markers):
                j=int(partner[i])
                row=dict(resource=panel.name,sample_group=group,marker_key=m['marker_key'],chromosome=m['chromosome'],
                         position_bp=m['position_bp'],window_bp=w,eligible_other_markers=int(eligible[i]),
                         valid_other_markers=int(nvalid[i]),best_partner_key=panel.markers[j]['marker_key'] if j>=0 else '',
                         best_partner_distance_bp=int(bdist[i]) if j>=0 else None,
                         best_r2=float(best[i]) if j>=0 else None,
                         status='EVALUABLE' if j>=0 else 'NO_VALID_PAIR' if eligible[i] else 'NO_OTHER_MARKER_IN_WINDOW')
                rows.append(row)
            ev=nvalid>0
            summary=dict(resource=panel.name,sample_group=group,window_bp=w,typed_markers=n,
                         markers_with_other_marker=int(np.sum(eligible>0)),markers_with_valid_pair=int(ev.sum()),
                         markers_without_any_neighbor=int(np.sum(eligible==0)),
                         markers_with_neighbors_but_no_valid_pair=int(np.sum((eligible>0)&~ev)),
                         interpretation='OTHER_OBSERVED_MARKERS_ONLY; NO_GENE_OR_CAUSAL_VARIANT_TAGGING_CLAIM')
            for t in cfg['ld']['r2_thresholds']:
                k=str(t).replace('.','_');nt=int(np.sum(ev&(best>=t)))
                summary['markers_best_r2_ge_'+k]=nt
                summary['fraction_all_markers_ge_'+k]=nt/n if n else None
                summary['fraction_evaluable_ge_'+k]=nt/int(ev.sum()) if ev.any() else None
            summaries.append(summary)
    return rows,summaries


def common_comparison(a:Panel,b:Panel,cfg,logger):
    a,b,audit,shared,omit=prepare_common(a,b,cfg)
    ca,pa=coords(a)[2:];cb,pb=coords(b)[2:]
    if not np.array_equal(ca,cb) or not np.array_equal(pa,pb):raise DataError('Shared final coordinate alignment failed.')
    pairs,flags,pop=pair_inventory(ca,pa,ca,pa,cfg['ld'],cfg['seed'],'shared_fixed_pair_comparison')
    logger.info('Shared comparison: %d fixed markers, %d/%d samples, %d pairs',len(a.markers),len(shared),len(omit),len(pairs))
    vals={}
    for p in [a,b]:
        vals[p.name]={g:pairwise_ld(p.dosage[:,idx],pairs,cfg['qc']['min_pair_samples'],cfg['ld']['batch_pairs']) for g,idx in p.groups.items()}
    same,dist=geometry(pairs,ca,pa);rows=[];cell_summary=[]
    for group in a.groups:
        aa,bb=vals[a.name][group],vals[b.name][group]
        x,y=a.dosage[:,a.groups[group]],b.dosage[:,b.groups[group]]
        both=np.isfinite(x)&np.isfinite(y)
        cell_summary.append(dict(sample_group=group,shared_markers=len(a.markers),shared_samples=len(a.groups[group]),
                                 nonmissing_symbol_pairs=int(both.sum()),discordant_dosage_cells=int(np.sum(both&(x!=y))),
                                 differing_missingness=int(np.sum(np.isfinite(x)!=np.isfinite(y)))))
        masks=[('all_comparison_pairs',None,None,np.ones(len(pairs),dtype=bool))]
        masks += [('distance_bin',lo,hi,same&(dist>lo)&(dist<=hi)) for lo,hi in zip(cfg['ld']['distance_edges_bp'][:-1],cfg['ld']['distance_edges_bp'][1:])]
        masks += [(lab,None,None,flags[lab]) for lab in ['interchromosomal_context','distant_context']]
        for label,lo,hi,mask in masks:
            valid=mask&np.isfinite(aa.r2)&np.isfinite(bb.r2);d=bb.r2[valid]-aa.r2[valid]
            rows.append(dict(sample_group=group,pair_class=label,lower_bp_exclusive=lo,upper_bp_inclusive=hi,
                resource_a=a.name,resource_b=b.name,shared_markers=len(a.markers),shared_samples=len(a.groups[group]),
                inventory_pairs=int(mask.sum()),valid_both=int(valid.sum()),
                validity_disagreement=int(np.sum(mask&(np.isfinite(aa.r2)!=np.isfinite(bb.r2)))),
                mean_r2_a=float(aa.r2[valid].mean()) if valid.any() else None,
                mean_r2_b=float(bb.r2[valid].mean()) if valid.any() else None,
                mean_delta_b_minus_a=float(d.mean()) if len(d) else None,
                max_abs_delta_r2=float(np.max(np.abs(d))) if len(d) else None,
                pairs_disagree_above_tolerance=int(np.sum(np.abs(d)>cfg['ld']['r2_agreement_tolerance'])),
                note='SAME_COLLECTION_FILTER_COMPARISON; NOT_INDEPENDENT_REPLICATION'))
    return dict(a=a,b=b,pairs=pairs,flags=flags,population=pop,values=vals,audit=audit,summary=rows,
                genotype_comparison=cell_summary,shared_samples=shared,shared_without_flagged=omit)


def scalar(v):
    if v is None:return None
    if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,np.bool_):return bool(v)
    return v


STATUSES={0:'VALID',1:'INSUFFICIENT_COMPLETE_SAMPLES',2:'ZERO_PAIRWISE_COMPLETE_VARIANCE'}


def pair_rows(result,lo=0,hi=None):
    p=result['panel'];pairs=result['pairs'];hi=len(pairs) if hi is None else hi
    b_same,b_dist=result['geometries']['original_position'];c_same,c_dist=result['geometries']['corrected_snp']
    for z in range(lo,hi):
        i,j=map(int,pairs[z]);a,b=p.markers[i],p.markers[j]
        row=dict(resource=p.name,marker_a=a['marker_key'],marker_b=b['marker_key'],source_row_a=a['original_row'],source_row_b=b['original_row'],
                 corrected_chr_a=a['chromosome'],corrected_pos_a=a['position_bp'],corrected_chr_b=b['chromosome'],corrected_pos_b=b['position_bp'],
                 corrected_distance_bp=int(c_dist[z]) if c_same[z] else None,
                 original_chr_a=a['original_chromosome'],original_pos_a=a['original_position_bp'],original_chr_b=b['original_chromosome'],original_pos_b=b['original_position_bp'],
                 original_distance_bp=int(b_dist[z]) if b_same[z] else None)
        row.update({k:bool(v[z]) for k,v in result['flags'].items()})
        for group,v in result['values'].items():
            prefix='all' if group=='all_qc_samples' else 'omit_flagged'
            row[prefix+'_n']=int(v.n[z]);row[prefix+'_r2']=scalar(v.r2[z]);row[prefix+'_status']=STATUSES[int(v.status[z])]
        yield row


def common_pair_rows(result,lo=0,hi=None):
    a,b=result['a'],result['b'];pairs=result['pairs'];hi=len(pairs) if hi is None else hi
    cc,pp=coords(a)[2:];same,dist=geometry(pairs,cc,pp)
    for z in range(lo,hi):
        i,j=map(int,pairs[z]);row=dict(marker_a=a.markers[i]['marker_key'],marker_b=a.markers[j]['marker_key'],
            chromosome_a=str(cc[i]),position_a=int(pp[i]),chromosome_b=str(cc[j]),position_b=int(pp[j]),
            distance_bp=int(dist[z]) if same[z] else None)
        row.update({k:bool(v[z]) for k,v in result['flags'].items()})
        for group in a.groups:
            va,vb=result['values'][a.name][group],result['values'][b.name][group]
            row[group+'_n_a']=int(va.n[z]);row[group+'_n_b']=int(vb.n[z])
            row[group+'_r2_a']=scalar(va.r2[z]);row[group+'_r2_b']=scalar(vb.r2[z])
            row[group+'_r2_b_minus_a']=scalar(vb.r2[z]-va.r2[z])
        yield row
