"""Complete dosage correlation workbooks from fixed in-memory cohorts."""
from __future__ import annotations
from .readers import read_matrix
from .genotypes import raw_concordance,prepare_panel
from .analysis import panel_analysis,sample_sensitivity,neighbor_diagnostics,common_comparison,pair_rows,common_pair_rows
from ..io.export import RowStream

def split_pair_sheets(results,cfg):
    sheets={};parts=[];chunk=cfg['ld']['export_rows_per_sheet']
    for no,r in enumerate(results,1):
        n=len(r['pairs'])
        columns=list(next(pair_rows(r,0,1))) if n else ['No_records']
        if n==0:sheets[f'panel_{no}_pairs']=[];continue
        for part,lo in enumerate(range(0,n,chunk),1):
            hi=min(lo+chunk,n);name=f'panel_{no}_pairs_{part}'
            sheets[name]=RowStream(hi-lo,lambda rr=r,l=lo,h=hi:pair_rows(rr,l,h),columns)
            parts.append(dict(sheet=name,resource=r['panel'].name,rows=hi-lo,pair_start_0based=lo,pair_end_exclusive=hi))
    return sheets,parts


def generate(source_files, ledger, cfg, save, logger):
    mats=[];panels=[];qc={k:[] for k in ['marker_qc','sample_qc','symbol_counts','unsupported_calls']}
    for name in cfg['resources']:
        logger.info('Reading original genotype matrix: %s',name)
        m=read_matrix(source_files[name],name,{})
        p,rr=prepare_panel(m,ledger,cfg);mats.append(m);panels.append(p)
        for k in qc:qc[k].extend(rr[k])
        logger.info('Fixed QC %s: %d markers, %d inclusive samples, %d sensitivity samples',name,len(p.markers),len(p.groups['all_qc_samples']),len(p.groups['omit_flagged_samples']))
    concord=raw_concordance(*mats,cfg)
    del mats
    states=[dict(resource=p.name,**p.metadata,status='CONDITIONAL_UNPHASED_DOSAGE_LD',source_sample_sensitivity=';'.join(cfg['sensitivity_sample_ids'])) for p in panels]
    save('01_genotype_coordinate_qc.xlsx',dict(panel_status=states,**qc,
         raw_sample_overlap=concord['sample_overlap'],raw_concordance=concord['summary'],
         raw_sample_concordance=concord['sample_concordance'],raw_marker_concordance=concord['marker_concordance'],raw_discordant_cells=concord['discordant_cells']))
    analyses=[panel_analysis(p,cfg,logger) for p in panels]
    save('02_ld_distance_profiles.xlsx',dict(
         distance_bins=[r for a in analyses for r in a['profiles']],
         window_summaries=[r for a in analyses for r in a['windows']],
         distant_context=[r for a in analyses for r in a['context']],
         coordinate_window_changes=[r for a in analyses for r in a['coordinate_window_changes']],
         pair_inventories=[dict(resource=a['panel'].name,**a['population']) for a in analyses],
         sample_sensitivity=[r for a in analyses for r in sample_sensitivity(a,cfg)]))
    sheets,parts=split_pair_sheets(analyses,cfg)
    save('03_pairwise_dosage_ld.xlsx',dict(sheet_index=parts,**sheets))
    common=common_comparison(*panels,cfg,logger)
    csheets={};n=len(common['pairs']);chunk=cfg['ld']['export_rows_per_sheet'];cp=[]
    if n:
        columns=list(next(common_pair_rows(common,0,1)))
        for part,lo in enumerate(range(0,n,chunk),1):
            hi=min(lo+chunk,n);nm=f'common_pairs_{part}'
            csheets[nm]=RowStream(hi-lo,lambda l=lo,h=hi:common_pair_rows(common,l,h),columns)
            cp.append(dict(sheet=nm,rows=hi-lo,pair_start_0based=lo,pair_end_exclusive=hi))
    save('04_shared_collection_comparison.xlsx',dict(
         summary=common['summary'],genotype_comparison=common['genotype_comparison'],shared_marker_audit=common['audit'],
         shared_samples=[dict(sample_id=s,in_shared_all=True,in_shared_omit_flagged=s in set(common['shared_without_flagged'])) for s in common['shared_samples']],
         pair_inventory=[common['population']],sheet_index=cp,**csheets))
    tagrows=[];tagsum=[]
    for a in analyses:
        rr,ss=neighbor_diagnostics(a,cfg);tagrows.extend(rr);tagsum.extend(ss)
    save('05_observed_snp_neighbors.xlsx',dict(summary=tagsum,marker_neighbors=tagrows))
    return states, concord['summary']
