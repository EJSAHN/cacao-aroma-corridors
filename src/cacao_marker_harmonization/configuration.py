"""Explicit and validated analytical policies."""
from __future__ import annotations
import json
from pathlib import Path
from .errors import DataError
from .ld.inputs import validate_config as validate_ld

def validate_coverage(cfg):
    required=['study_id','seed','required_assembly_accession','comparison_resources',
              'original_only_resources','windows_bp','rarefaction_replicates','rarefaction_query_sets','background','minimum_output_free_gib']
    for k in required:
        if k not in cfg:raise DataError('Missing configuration key: '+k)
    for k in ['seed','rarefaction_replicates']:
        if isinstance(cfg[k],bool) or not isinstance(cfg[k],int) or cfg[k]<0 or (k!='seed' and cfg[k]==0):
            raise DataError('Invalid configuration: '+k)
    w=cfg['windows_bp']
    if not isinstance(w,list) or not w or w[0]!=0 or any(isinstance(x,bool) or not isinstance(x,int) or x<0 for x in w) or w!=sorted(set(w)):
        raise DataError('Windows must be unique increasing nonnegative integer flank lengths starting at zero.')
    rs=cfg['comparison_resources'];ctx=cfg['original_only_resources']
    if len(rs)<2 or len(rs+ctx)!=len(set(rs+ctx)) or any(not isinstance(x,str) or not x for x in rs+ctx):
        raise DataError('Invalid resource lists.')
    bg=cfg['background']
    if bg.get('formal_p_values') is not False:
        raise DataError('Formal p-values are deliberately not supported by this descriptive reference contrast.')
    for k in ['replicates','length_bins','local_gene_density_bins']:
        if isinstance(bg.get(k),bool) or not isinstance(bg.get(k),int) or bg[k]<1:raise DataError('Invalid background '+k)
    if not isinstance(bg.get('enabled'),bool):raise DataError('Background enabled must be Boolean.')
    models=bg.get('models')
    if not models or len(models)!=len(set(models)) or not set(models).issubset({'chromosome_biotype_length','chromosome_biotype_length_gene_density'}):
        raise DataError('Invalid background models.')
    if bg.get('query_set')!='unique_candidate_genes':raise DataError('The primary background anchor must remain the complete fixed candidate gene set.')
    if not isinstance(bg.get('local_gene_density_radius_bp'),int) or isinstance(bg['local_gene_density_radius_bp'],bool) or bg['local_gene_density_radius_bp']<0:
        raise DataError('Invalid annotation neighborhood radius.')
    if not cfg['rarefaction_query_sets'] or len(cfg['rarefaction_query_sets'])!=len(set(cfg['rarefaction_query_sets'])):raise DataError('Invalid rarefaction queries.')
    if not isinstance(cfg['minimum_output_free_gib'],(int,float)) or not 0<cfg['minimum_output_free_gib']<1000:raise DataError('Invalid free-space threshold.')
    return cfg


def validate_alignment(cfg):
    a=cfg.get('alignment',{})
    for k in ['seed','threads','max_reported_hits','seed_length','seed_mismatches','effort_D','effort_R','mismatch_penalty',
              'gap_open','gap_extend','max_nonfocal_edit_bases','max_gap_bases','focal_gap_guard_bases','minimum_score_gap','minimum_tag_length','timeout_seconds']:
        v=a.get(k)
        if isinstance(v,bool) or not isinstance(v,int) or v<0:raise DataError('Invalid alignment parameter: '+k)
    if not 1<=a['threads']<=8 or not 2<=a['max_reported_hits']<=20 or not 4<=a['seed_length']<=32 or a['seed_mismatches'] not in (0,1):raise DataError('Alignment configuration exceeds supported bounds.')
    if a['minimum_score_gap']<1 or a['minimum_tag_length']<20 or a['timeout_seconds']<1:raise DataError('Invalid acceptance/timeout threshold.')
    if a.get('score_min')!='L,-30,0':raise DataError('This release requires the explicit constant -30 search cutoff for score-gap censoring.')
    if a['version']!='2.5.5':raise DataError('This release is pinned to Bowtie 2 2.5.5.')
    if cfg['background'].get('query_sets')!=['unique_candidate_genes','protein_coding_candidate_genes']:raise DataError('Both full and protein-coding candidate sensitivity sets must be preserved.')
    return cfg


def load_study(path: Path) -> dict:
    cfg=json.loads(path.read_text(encoding="utf-8-sig"))
    required={'schema_version','sources','source_evidence','reference','exact_mapping','coverage','alignment','tools','dosage_ld','minimum_output_free_gib'}
    if set(cfg)!=required or cfg['schema_version']!=1:
        raise DataError('Configuration schema differs from this release.')
    validate_coverage(cfg['coverage'])
    validate_alignment(dict(cfg['coverage'],alignment=cfg['alignment']))
    validate_ld(cfg['dosage_ld'])
    if cfg['coverage']['required_assembly_accession']!=cfg['reference']['assembly_accession'] or cfg['dosage_ld']['required_assembly_accession']!=cfg['reference']['assembly_accession']:
        raise DataError('Reference accession must agree across analytical policies.')
    roles=[r['id'] for r in cfg['sources']['resources']]
    if len(roles)!=len(set(roles)):
        raise DataError('Duplicate source resource identifiers.')
    if set(cfg['reference']['chromosomes'])!=set(cfg['sources']['chromosomes']):
        raise DataError('Source and reference chromosome inventories differ.')
    if cfg['source_evidence']['apply_coordinate_corrections'] or cfg['source_evidence']['infer_amazonia_encoding']:
        raise DataError('No inferred source coding or unverified coordinate overrides are supported.')
    from .evidence.events import validated_aliases
    validated_aliases(cfg['source_evidence'])
    return cfg
