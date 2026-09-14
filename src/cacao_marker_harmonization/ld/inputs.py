"""Configuration validation for fixed-coordinate dosage correlations."""
from __future__ import annotations
from pathlib import Path
from .common import DataError, text
from ..io.tableio import TableBook


def read_sheet(path: Path, sheet: str, required=()):
    with TableBook(path) as book:
        if sheet not in book.sheet_names:
            raise DataError(f'Required worksheet missing: {path.name}:{sheet}')
        out = [r for _, r in book.table(sheet, required=required)]
        if book.errors:
            raise DataError(f'Excel error cell in {path.name}:{sheet}')
    return out


def validate_config(cfg):
    required = {'study_id',
                'required_assembly_accession','seed','resources','chromosomes',
                'sensitivity_sample_ids','sensitivity_reason','encoding','qc','ld','minimum_output_free_gib'}
    if set(cfg) != required:
        raise DataError('Configuration keys differ: ' + str(sorted(set(cfg) ^ required)))
    def whole(v): return isinstance(v, int) and not isinstance(v, bool)
    if not whole(cfg['seed']) or cfg['seed'] < 0:
        raise DataError('Seed must be a nonnegative integer.')
    for k in ['resources','chromosomes','sensitivity_sample_ids']:
        v = cfg[k]
        if not isinstance(v, list) or not v or any(not isinstance(x, str) or not x for x in v) or len(v) != len(set(v)):
            raise DataError('Configuration list must contain distinct nonempty strings: ' + k)
    if len(cfg['resources']) != 2:
        raise DataError('This release requires exactly two filter resources for the fixed shared-sample comparison.')
    if not text(cfg['sensitivity_reason']): raise DataError('Sample sensitivity needs a recorded reason.')
    e = cfg['encoding']
    if e.get('mode') != 'ACM_symbol_dosage' or set(e.get('declared_allele_set', [])) != {'A','C'}:
        raise DataError('Unsupported genotype encoding; no inference or recoding is allowed.')
    if e.get('symbol_to_dosage') not in ({'A':0,'M':1,'C':2},{'A':2,'M':1,'C':0}):
        raise DataError('A/C/M coding must keep M=1 and opposite homozygotes=0/2.')
    if set(e.get('missing_symbols', [])) & set(e['symbol_to_dosage']):
        raise DataError('Missing symbols overlap genotype states.')
    if not text(e.get('interpretation')):
        raise DataError('Genotype encoding limitations must be declared.')
    q = cfg['qc']
    for k in ['sample_min_call_rate','marker_min_call_rate']:
        if not 0 < q[k] <= 1: raise DataError('Call rate must be in (0,1].')
    if not 0 < q['min_maf'] <= 0.5: raise DataError('MAF threshold must be in (0,0.5].')
    for k in ['min_pair_samples','min_markers']:
        if not whole(q[k]) or q[k] < 3: raise DataError('Invalid QC minimum: ' + k)
    if q['duplicate_corrected_position_policy'] != 'exclude_all' or q['fixed_marker_qc'] != 'intersection_of_all_and_omitted_sample_groups':
        raise DataError('This release only implements conservative duplicate exclusion and fixed intersection QC.')
    ld = cfg['ld']
    for k in ['max_distance_bp','distant_same_chromosome_min_bp','max_union_pairs','batch_pairs','export_rows_per_sheet']:
        if not whole(ld[k]) or ld[k] <= 0: raise DataError('Invalid positive integer: ' + k)
    if not whole(ld['context_pairs_per_class']) or ld['context_pairs_per_class'] < 0:
        raise DataError('Invalid context pair count.')
    edges = ld['distance_edges_bp']; windows = ld['report_windows_bp']; thresholds = ld['r2_thresholds']
    if (not isinstance(edges,list) or len(edges) < 2 or edges[0] != 0
            or edges[-1] != ld['max_distance_bp'] or any(not whole(x) for x in edges)
            or sorted(set(edges)) != edges): raise DataError('Distance edges must increase from zero to max distance.')
    if not windows or sorted(set(windows)) != windows or any(not whole(x) or x <= 0 or x > ld['max_distance_bp'] for x in windows):
        raise DataError('Window lengths must be distinct increasing positive bp, <= max distance.')
    if not thresholds or sorted(set(thresholds)) != thresholds or any(not isinstance(x,(int,float)) or not 0 <= x <= 1 for x in thresholds):
        raise DataError('Invalid r2 thresholds.')
    if ld['distant_same_chromosome_min_bp'] <= ld['max_distance_bp']:
        raise DataError('Distant context must lie beyond the exhaustive local-pair window.')
    if ld['export_rows_per_sheet'] > 1048575 or ld['batch_pairs'] > 65536:
        raise DataError('Excel row/batch safety limits exceeded.')
    if not 0 < ld['r2_agreement_tolerance'] <= 1e-6:
        raise DataError('Invalid numerical comparison tolerance.')
    if cfg['minimum_output_free_gib'] <= 0: raise DataError('Output space requirement must be positive.')
    return cfg








