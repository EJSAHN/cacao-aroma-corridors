"""Raw-input analysis with explicit source, coordinate, coverage and dosage evidence."""
from __future__ import annotations
from collections import Counter
from datetime import datetime
import gc
import hashlib
import importlib.metadata
import json
import logging
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import zipfile
import numpy as np
from . import __version__
from .errors import DataError
from .sources.common import sha256
from .sources.sources import locate, load_associations, load_candidates, load_resource
from .sources.audit import resource_summary, sample_overlap, identity_audit, genotype_symbols, gene_interval_links
from .evidence.readers import read_diversity, read_matrix, read_associations
from .evidence.sequences import sequence_evidence, coordinate_hypotheses
from .evidence.events import reconcile
from .evidence.genotypes import symbol_diagnostics, compare_matrices
from .reference.reference import acquire_asset, iter_fasta, validate_sequence_header, validate_fasta_inventory, read_gff
from .reference.matcher import build_queries as build_tags, ExactTagMatcher, summarize_mappings
from .reference.harmonize import project_markers, candidate_checks, peak_links
from .coverage.inputs import build_cohorts, build_queries
from .coverage.analysis import representation_analysis, count_matching, cohort_summaries
from .alignment.refine import add_coding_query, project, mapping_transitions, background_sets, compare_previous
from .alignment.tags import source_tags, write_fasta
from .alignment.sam import parse_sam
from .alignment.smoke import native_smoke
from .alignment.toolchain import resolve_tools, reference_cache, build_index, align_cached, Reference, free_check
from .io.export import write_workbook
from .ld.reports import generate as generate_ld


def primitive(value):
    """Match the source-workbook missing-value serialization without numeric rounding."""
    if isinstance(value, np.generic): value=value.item()
    if isinstance(value,float) and not math.isfinite(value): return None
    if isinstance(value,dict):return {k:primitive(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [primitive(v) for v in value]
    return value


def records(value):
    if hasattr(value,'to_dict'):return [primitive(x) for x in value.to_dict('records')]
    return primitive(value)


def source_paths(root: Path,cfg:dict):
    raw=root/'data_raw' if (root/'data_raw').is_dir() else root
    if not raw.is_dir():raise DataError('Source directory not found: '+str(raw))
    s=cfg['sources'];paths={}
    paths['associations']=locate(raw,s['association_file'],s['association_aliases'])
    paths['candidates']=locate(raw,s['candidate_file'],s['candidate_aliases'])
    for spec in s['resources']:
        p=locate(raw,spec['path'],spec['aliases'],spec.get('optional',False))
        if p:paths[spec['id']]=p
    for p in paths.values():
        if not p.resolve().is_relative_to(raw.resolve()):raise DataError('Source path resolves outside the selected read-only root.')
    return raw,paths


def report_writer(root,log,scope):
    root.mkdir(parents=True,exist_ok=True)
    notes=[('Software','Cacao marker harmonization '+__version__),('Scope',scope),
           ('Coordinates','1-based inclusive intervals. Flanks are measured in bp on each side. Peak exact means same position; gene exact means a gene-body hit.'),
           ('Inference','Positional representation and unadjusted genotype correlations do not establish aroma effects, optimal physical windows or breeding prediction.'),
           ('Missing','Blank numeric values indicate not available or not evaluable, not zero.'),
           ('Preservation','Input workbooks are read-only. Qualitative annotations and excluded records remain documented.')]
    def save(name,tables,extra=()):
        write_workbook(root/name,tables,notes+list(extra))
        log.info('Saved %s/%s (%.2f MiB)',root.name,name,(root/name).stat().st_size/2**20)
    return save


def parse_inputs(paths,cfg,save,log):
    s=cfg['sources'];log.info('Parsing original association records and candidate annotations')
    assoc=load_associations(paths['associations'],s);cand=load_candidates(paths['candidates'],s)
    resources=[]
    for spec in s['resources']:
        if spec['id'] in paths:
            log.info('Reading source resource: %s',spec['id']);resources.append(load_resource(paths[spec['id']],spec,s))
    ids,coord_summary,coord_details,collisions=identity_audit(resources)
    totals,patterns=genotype_symbols(resources)
    checkpoints=[dict(metric=k,count=v) for k,v in [
        ('association_raw_rows',len(assoc['raw_records'])),('master_raw_rows',len(assoc['master_raw'])),
        ('master_unique_event_keys',len(assoc['master_events'])),('master_unique_peak_positions',len(assoc['unique_peaks'])),
        ('fruity_subset_unique_peaks',len(assoc['fruity_peaks'])),('master_unique_valid_intervals',len(assoc['source_intervals'])),
        ('candidate_raw_rows',len(cand['raw_records'])),('candidate_unique_valid_intervals',len(cand['unique_genes']))]]
    save('01_resources.xlsx',dict(checkpoints=checkpoints,resource_summary=records(resource_summary(resources)),
         sample_overlap=records(sample_overlap(resources)),sample_ids=[dict(resource=r.spec['id'],sample_id=x) for r in resources for x in r.samples],
         genotype_symbols=records(totals),genotype_patterns=records(patterns),source_metadata=[x for r in resources for x in r.metadata]))
    aa={k:records(v) for k,v in assoc.items()};cc={k:records(v) for k,v in cand.items()}
    save('02_association_units.xlsx',{k:aa[k] for k in ['worksheet_audit','raw_records','master_events','extra_records','interval_audit','unique_peaks','fruity_peaks','source_intervals','union_peaks_sensitivity']})
    save('03_candidate_units.xlsx',dict(cc,source_interval_overlap=records(gene_interval_links(cand['unique_genes'],assoc['source_intervals']))))
    inventory=[x for r in resources for x in records(r.markers)]
    save('04_marker_identity.xlsx',dict(id_overlap=records(ids),coordinate_summary=records(coord_summary),matched_details=records(coord_details),ambiguous_keys=records(collisions),marker_inventory=inventory))
    return aa,cc,inventory,checkpoints


def source_evidence(paths,cfg,save,log):
    e=cfg['source_evidence'];log.info('Checking source tag identity and annotation correspondence')
    diversity,refs,meta,errors=read_diversity(paths['Diversity'],e)
    seq,summary,anomalies=sequence_evidence(diversity,refs)
    del refs,diversity
    save('01_sequence_identity.xlsx',dict(sequence_summary=summary,sequence_checks=seq,reference_anomalies=anomalies,study_metadata=meta))
    mats={name:read_matrix(paths[name],name,e) for name in [*e['nacional_resources'],'Amazonia']}
    detail,cs,ms,amb,review=coordinate_hypotheses(seq,mats,e)
    save('02_coordinate_hypotheses.xlsx',dict(pair_summary=cs,model_summary=ms,pair_details=detail,ambiguous_keys=amb,mapping_review_candidates=review))
    arr,headers,err=read_associations(paths['associations'],e);errors.extend(err)
    ev=reconcile(arr,e)
    save('03_association_correspondence.xlsx',dict(ev,source_headers=headers,original_records=arr))
    gd=symbol_diagnostics(mats,e);concord=compare_matrices(*(mats[n] for n in e['nacional_resources']),e)
    save('04_genotype_symbols.xlsx',dict(gd,**concord,source_read_errors=errors+[dict(resource=name,**row) for name,mat in mats.items() for row in mat.errors]))
    del mats,detail,ev,arr;gc.collect()
    return seq


def exact_mapping(seqrows,inventory,candidates,peaks,fasta,cfg,save,log):
    ref=cfg['reference'];ms=cfg['exact_mapping']['matcher']
    direct={int(r['source_row']):r for r in inventory if r['resource']==cfg['exact_mapping']['direct_resource']}
    if len(direct)!=len(seqrows):raise DataError('Marker-row count differs across source readers.')
    for r in seqrows:
        old=direct.get(int(r['source_row']))
        if old is None or any(str(old.get(a) or '').strip()!=str(r.get(b) or '').strip() for a,b in [('marker_id','marker_id'),('reference_key','reference_key'),('flank_5','five_flank'),('flank_3','three_flank')]):
            raise DataError('Source readers disagree on tag identity.')
    qs,audit=build_tags(seqrows,ms['minimum_tag_length'],ms['seed_k'])
    matcher=ExactTagMatcher(qs,ms['seed_k'],ms['chunk_bases'],ms['max_hits']);hits=[];inv=[]
    last=time.monotonic()
    def progress(seqid,done,total):
        nonlocal last
        now=time.monotonic()
        if now-last>=10:log.info('Exact search %s: %.1f / %.1f Mb',seqid,done/1e6,total/1e6);last=now
    for sid,header,sequence in iter_fasta(fasta,ref['maximum_contig_bases']):
        validate_sequence_header(sid,header,ref)
        inv.append(dict(sequence_id=sid,chromosome=ref['chromosome_aliases'].get(sid),length_bp=len(sequence),sequence_sha256=hashlib.sha256(sequence).hexdigest(),original_header=header))
        hits.extend(matcher.scan(sid,sequence,progress))
    validate_fasta_inventory(inv,ref)
    mapping=summarize_mappings(audit,hits,ref['chromosome_aliases']);sc=Counter(r['mapping_status'] for r in mapping)
    save('01_exact_tag_mapping.xlsx',dict(mapping_summary=[dict(status=s,records=n,fraction=n/len(mapping)) for s,n in sorted(sc.items())],
         marker_mapping=mapping,all_exact_placements=[dict(**h.as_dict(),chromosome=ref['chromosome_aliases'].get(h.sequence_id)) for h in hits]))
    panels,summary=project_markers(inventory,mapping,cfg['exact_mapping']['direct_resource'],cfg['exact_mapping']['linked_resources'])
    accepted=[r for r in panels if r['accepted_for_coordinate_comparison']]
    save('02_exact_coordinates.xlsx',dict(panel_summary=summary,all_panel_records=panels,accepted_coordinate_records=accepted))
    models=[]
    for resource in sorted({r['resource'] for r in panels}):
        rr=[r for r in accepted if r['resource']==resource]
        ct=Counter((r.get('source_tag_strand_on_reference'),r.get('source_minus_reference_snp_bp'),r.get('source_equals_tag_start_1based'),r.get('source_equals_tag_start_0based'),r.get('source_chromosome_agrees')) for r in rr)
        for (strand,offset,one,zero,same),n in sorted(ct.items(),key=lambda x:str(x[0])):
            models.append(dict(resource=resource,mapped_strand=strand,source_minus_reference_snp_bp=offset,source_equals_tag_start_1based=one,source_equals_tag_start_0based=zero,same_chromosome=same,records=n))
    save('03_coordinate_conventions.xlsx',dict(offset_strand_summary=models,changed_coordinate_records=[r for r in accepted if r.get('source_chromosome_agrees') is not True or r.get('source_minus_reference_snp_bp')!=0]))
    return inv,mapping,panels


def annotate(gff,inventory,candidates,cfg,save,assets):
    genes,headers,issues,declared=read_gff(gff,{r['sequence_id']:r['length_bp'] for r in inventory},cfg['reference']['chromosome_aliases'],cfg['reference'])
    if not declared:raise DataError('GFF3 does not declare the expected assembly; candidate intervals cannot be approved.')
    checks=candidate_checks(candidates,genes,cfg['reference']['chromosome_aliases'])
    ct=Counter(r['status'] for r in checks)
    save('04_reference_annotation.xlsx',dict(reference_assets=[a.report() for a in assets],sequence_inventory=inventory,
         candidate_summary=[dict(status=s,records=n) for s,n in sorted(ct.items())],candidate_interval_check=checks,gene_catalog=genes,annotation_headers=headers,annotation_issues=issues))
    lengths={r['chromosome']:r['length_bp'] for r in inventory if r['chromosome'] is not None}
    return genes,checks,lengths


def make_queries(assoc,cand,checks,lengths):
    source_queries={name:assoc[sheet] for name,sheet in [('master_unique_peaks','unique_peaks'),('master_fruity_subset_peaks','fruity_peaks'),('master_source_intervals','source_intervals')]}
    queries,audit,allids=build_queries(cand['unique_genes'],checks,cand['memberships'],source_queries,lengths)
    queries,coding=add_coding_query(queries)
    return queries,audit,allids,coding


def native_alignment(seqrows,exact,seq_inventory,fasta_asset,cfg,args,out,save,log):
    acfg=dict(cfg['coverage'],alignment=cfg['alignment'],tools=cfg['tools'])
    tags,audit=source_tags(seqrows,exact,cfg['alignment']['minimum_tag_length'])
    tools=resolve_tools(args.tools_root,acfg,args.bowtie2_dir,args.offline,log)
    env=dict(os.environ);env['PYTHONUNBUFFERED']='1'
    threads=args.threads or cfg['alignment']['threads']
    smoke=native_smoke(tools,args.cache_root,threads,env,acfg,log)
    (out/'native_smoke_test.json').write_text(json.dumps(smoke,indent=2)+'\n', encoding='utf-8')
    fasta,indexjson=reference_cache(fasta_asset.path,fasta_asset.report(),seq_inventory,args.cache_root,log)
    query=out/'source_allele_tags.fa';write_fasta(query,tags)
    prefix=build_index(fasta,args.cache_root,tools,threads,env,acfg,log)
    sam,identity=align_cached(query,prefix,args.cache_root,tools,threads,env,acfg,log)
    align_prov=dict(sam_filename=sam.name,sam_sha256=sha256(sam),reference_fasta_sha256=sha256(fasta),
                   alignment_cache_identity=identity,native_tool_versions={name: value.replace(str((tools['aligner'].parent/name)), name) for name,value in tools['versions'].items()},tool_sha256=tools['sha256'])
    (out/'alignment_provenance.json').write_text(json.dumps(align_prov,indent=2)+'\n', encoding='utf-8')
    shutil.copyfile(sam.parent/'alignment.log',out/'bowtie2_alignment.log')
    # Local execution paths are kept outside the scientific archive.
    (args.output_root/'local_execution_paths.json').write_text(json.dumps(dict(reference_fasta=str(fasta_asset.path),reference_cache=str(fasta),sam=str(sam),bowtie2_dir=str(tools['aligner'].parent),native_tool_versions=tools['versions']),indent=2)+'\n', encoding='utf-8')
    with Reference(fasta,indexjson) as reference:
        mapped,alleles,top,sam_records=parse_sam(sam,tags,reference,cfg['alignment'],{r['sequence_id']:r['chromosome'] for r in seq_inventory},log)
    active={r['query_id']:r for r in mapped};decisions=[]
    for a in audit:
        d=active.get(a['query_id'])
        if d is None:d=dict(query_id=a['query_id'],source_row=a['source_row'],reference_key=a['reference_key'],source_marker_id=a['marker_id'],placement_status=a['input_status'],accepted_for_mapping=False)
        decisions.append(d)
    changes=mapping_transitions(exact,decisions);ct=Counter(r['placement_status'] for r in decisions)
    save('01_realignment_evidence.xlsx',dict(mapping_summary=[dict(status=s,records=n) for s,n in sorted(ct.items())],marker_decisions=decisions,
         allele_decisions=alleles,top_two_per_allele=top,exact_policy_transitions=changes,input_tag_audit=audit),
         [('Alignment reporting','Two highest-score records per allele are listed here. The complete compressed SAM is retained in the external analysis cache.'),
          ('Search limits','Bounded end-to-end Bowtie 2 search, not proof of genome-wide uniqueness.')])
    return decisions,dict(tags=len(tags),sam_records=sam_records,mapping_statuses=dict(ct))


def evaluate_coverage(inventory,projections,queries,queryaudit,allids,coding,genes,lengths,cfg,exact_summary,save,log):
    cohorts,ledger,states=build_cohorts(inventory,projections,lengths,cfg['coverage'])
    save('02_fixed_marker_cohorts.xlsx',dict(panel_status=states,chromosome_counts=cohort_summaries(cohorts,lengths),marker_record_ledger=ledger,all_panel_coordinates=projections))
    qt={'query_inventory':[dict(query_set=n,n_queries=len(q)) for n,q in queries.items()],'candidate_set_membership':coding,'query_inclusion_audit':queryaudit}
    for n,q in queries.items():
        if not n.startswith('pathway:'):qt['protein_coding_candidates' if n=='protein_coding_candidate_genes' else n]=q.rows
    qt['pathway_membership']=[dict(gene_id=r['query_id'],pathway=n.split(':',1)[1]) for n,q in queries.items() if n.startswith('pathway:') for r in q.rows]
    save('03_fixed_query_sets.xlsx',qt)
    result,points,_=representation_analysis(cohorts,queries,cfg['coverage']['windows_bp'],log)
    result['exact_only_policy_comparison']=compare_previous(exact_summary,result['summary'])
    save('04_paired_coverage_effects.xlsx',result)
    save('05_count_matched_sensitivity.xlsx',count_matching(points,queries,cfg['coverage'],log))
    bg=background_sets(genes,allids,queries,lengths,points,cfg['coverage'],log)
    save('06_matched_gene_background.xlsx',bg)
    return ledger,states,bg['status']


def overlaps(a,b):
    a,b=a.resolve(),b.resolve()
    return a==b or a in b.parents or b in a.parents


def check_locations(args,paths):
    writes=[args.output_root,args.scratch,args.cache_root,args.tools_root,args.reference_cache]
    protected=[args.source_root,Path(__file__).resolve().parents[2]]
    for w in writes:
        for p in protected:
            if overlaps(w,p):raise DataError('Output/cache/scratch overlaps source or code: '+str(w))
        for p in paths.values():
            if overlaps(w,p):raise DataError('Write path overlaps an original input: '+str(w))
    for i,w in enumerate(writes):
        for other in writes[i+1:]:
            if overlaps(w,other):raise DataError('Output, caches, tools and scratch must be separate nonnested locations.')
    for p in (args.reference_fasta,args.reference_gff):
        if p:
            for w in (args.output_root,args.scratch,args.cache_root,args.tools_root):
                if overlaps(w,p):raise DataError('Reference input cannot overlap writable output/scratch/tool/alignment-cache paths.')


def run(args,cfg):
    started=time.monotonic();raw,paths=source_paths(args.source_root,cfg)
    check_locations(args,paths)
    for p in [args.output_root,args.scratch,args.cache_root,args.tools_root,args.reference_cache]:p.mkdir(parents=True,exist_ok=True)
    free_check(args.output_root,cfg['minimum_output_free_gib']);free_check(args.cache_root,cfg['alignment']['min_cache_free_gib'])
    out=args.output_root/('analysis_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir()
    log=logging.getLogger('cacao_marker_harmonization');log.setLevel(logging.INFO);log.propagate=False
    for h in list(log.handlers):h.close();log.removeHandler(h)
    for h in [logging.StreamHandler(sys.stdout),logging.FileHandler(out/'analysis.log',encoding='utf-8')]:
        h.setFormatter(logging.Formatter('%(asctime)s | %(levelname)s | %(message)s'));log.addHandler(h)
    oldenv={k:os.environ.get(k) for k in ['TEMP','TMP','TMPDIR']};oldtemp=tempfile.tempdir
    os.environ.update({k:str(args.scratch) for k in oldenv});tempfile.tempdir=str(args.scratch)
    source_manifest=[dict(role=k,relative_path=p.relative_to(raw).as_posix(),bytes=p.stat().st_size,sha256=sha256(p)) for k,p in paths.items()]
    protected={p:sha256(p) for p in paths.values()};protected[args.config]=sha256(args.config)
    (out/'study_config.json').write_text(json.dumps(cfg,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    try:
        log.info('Cacao marker harmonization %s | Python %s | NumPy %s',__version__,sys.version.split()[0],np.__version__)
        srcsave=report_writer(out/'source_data',log,'Source identifiers, annotation records and genotype symbols; no inferred source edits.')
        assoc,cand,inventory,checkpoints=parse_inputs(paths,cfg,srcsave,log)
        seqrows=source_evidence(paths,cfg,report_writer(out/'source_evidence',log,'Source sequence identity and diagnostic record reconciliation.'),log)
        spec=cfg['reference'];search=[args.reference_dir] if args.reference_dir else []
        search.extend([raw/'reference',raw/'reference_annotation'])
        fa=acquire_asset('fasta',spec,args.reference_cache,search,args.reference_fasta,args.offline,log)
        gff=acquire_asset('gff',spec,args.reference_cache,search,args.reference_gff,args.offline,log)
        assets=[fa,gff]
        for a in assets:protected[a.path]=a.sha256
        refsave=report_writer(out/'reference_mapping',log,'Exact full-tag mapping and independent source-gene annotation checks; original peaks retained.')
        seqinv,exact,exact_panels=exact_mapping(seqrows,inventory,cand['unique_genes'],assoc['unique_peaks'],fa.path,cfg,refsave,log)
        genes,checks,lengths=annotate(gff.path,seqinv,cand['unique_genes'],cfg,refsave,assets)
        refsave('05_association_context.xlsx',dict(peak_context=peak_links(assoc['unique_peaks'],inventory,exact_panels)))
        queries,qaudit,allids,coding=make_queries(assoc,cand,checks,lengths)
        exact_cohorts,_,_=build_cohorts(inventory,exact_panels,lengths,cfg['coverage'])
        # The exact-policy reference excludes the coding-only set, preserving its original comparison scope.
        exact_queries={k:v for k,v in queries.items() if k!='protein_coding_candidate_genes'}
        ex_result,_,_=representation_analysis(exact_cohorts,exact_queries,cfg['coverage']['windows_bp'],log)
        report_writer(out/'exact_policy_coverage',log,'Exact-match-only mapping policy sensitivity.')('01_coverage.xlsx',ex_result)
        del exact_panels,exact_cohorts,assoc;gc.collect()
        alignsave=report_writer(out/'alignment_coverage',log,'Fixed A/B/C positional representation, mapping policy, marker counts and descriptive gene backgrounds.')
        decisions,alignstats=native_alignment(seqrows,exact,seqinv,fa,cfg,args,out,alignsave,log)
        projections=project(inventory,decisions,cfg['coverage'])
        del seqrows,exact,decisions;gc.collect()
        ledger,states,bgstatus=evaluate_coverage(inventory,projections,queries,qaudit,allids,coding,genes,lengths,cfg,ex_result['summary'],alignsave,log)
        del inventory,projections,genes,queries,ex_result;gc.collect()
        ld_scope='Conditional unphased A/C/M dosage correlations; no imputation, haplotype phasing or structure correction. G7 and MAF5 are overlapping filter exports.'
        ld_states,concord=generate_ld(paths,ledger,cfg['dosage_ld'],report_writer(out/'dosage_ld',log,ld_scope),log)
        for p,h in protected.items():
            if not p.is_file() or sha256(p)!=h:raise DataError('Read-only input changed during execution: '+p.name)
        versions=[dict(package=k,version=importlib.metadata.version(k)) for k in ['numpy','pandas','openpyxl']]
        statuses=[dict(check='source_and_reference_integrity',status='PASS'),dict(check='coordinate_mapping',status='EXACT_AND_BOUNDED_ALIGNMENT_POLICIES_REPORTED'),
                  dict(check='coverage',status='FIXED_A_B_C_COMPARISON'),dict(check='gene_backgrounds',status='DESCRIPTIVE_MATCHED_ANNOTATION'),
                  dict(check='dosage_ld',status='CONDITIONAL_UNPHASED_UNADJUSTED'),dict(check='Amazonia',status='ORIGINAL_COORDINATE_CONTEXT_ONLY'),
                  dict(check='aroma_effects_and_breeding',status='NOT_TESTED')]
        report_writer(out/'provenance',log,'Input digests, declared settings, scientific restrictions and completion state.')('01_analysis_manifest.xlsx',dict(source_manifest=source_manifest,reference_assets=[a.report() for a in assets],software_versions=versions,
            checks=statuses,background_status=bgstatus,panel_status=states,ld_panel_status=ld_states,raw_concordance=concord))
        summary={'software_version':__version__,'status':'COMPLETED_WITH_DECLARED_LIMITATIONS','elapsed_seconds':round(time.monotonic()-started,3),'source_manifest':source_manifest,
                 'reference_assets':[a.report() for a in assets],'checkpoints':checkpoints,'alignment':alignstats,'coverage_panels':states,'ld_panels':ld_states,'background_status':bgstatus}
        (out/'analysis_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        lines=['CACAO MARKER HARMONIZATION','Version: '+__version__,'Status: '+summary['status'],f"Elapsed seconds: {summary['elapsed_seconds']}",'',
               'Original input and reference hashes: unchanged','Raw-input workflow: no earlier result folders required','Figures: not generated','',
               'Interpretation: A-to-B changes selection; B-to-C changes coordinates on identical records.',
               'Background draws are descriptive and do not reconstruct the original discovery universe.',
               'Dosage correlations are conditional on declared coding and are not phenotype validation.','',
               'External numerical reference comparison is a separate verification step.']
        (out/'RUN_SUMMARY.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        log.info('Completed raw-input workflow. External reference agreement is checked separately.')
    except BaseException as exc:
        (out/'FAILED.txt').write_text(type(exc).__name__+': '+str(exc)+'\n',encoding='utf-8');log.exception('Incomplete analysis; original inputs were not intentionally modified.');raise
    finally:
        for h in list(log.handlers):h.close();log.removeHandler(h)
        tempfile.tempdir=oldtemp
        for k,v in oldenv.items():
            if v is None:os.environ.pop(k,None)
            else:os.environ[k]=v
    finalize(out)
    (args.output_root/'latest_analysis_run.txt').write_text(str(out)+'\n',encoding='utf-8')
    print('[RUN FOLDER]',out,flush=True);print('[UPLOAD ZIP]',out.with_suffix('.zip'),flush=True)
    return out


def finalize(out):
    manifest={p.relative_to(out).as_posix():sha256(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_SHA256.json'}
    (out/'OUTPUT_SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    zpath=out.with_suffix('.zip');tmp=zpath.with_name(zpath.name+'.partial')
    if tmp.exists():raise DataError('Archive partial already exists: '+str(tmp))
    try:
        with zipfile.ZipFile(tmp,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for p in sorted(out.rglob('*')):
                if p.is_file():z.write(p,arcname=out.name+'/'+p.relative_to(out).as_posix())
        os.replace(tmp,zpath)
    finally:
        if tmp.exists():tmp.unlink()
