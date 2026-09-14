"""Independent, values-level comparison with external reference workbooks."""
from __future__ import annotations
from itertools import zip_longest
import json, math
from pathlib import Path
from .errors import DataError
from .io.tableio import TableBook
from .io.export import write_workbook
from .sources.common import sha256

COLUMN_ALIASES = {'baseline_coordinate_eligible':'source_coordinate_eligible',
                  'baseline_eligible_records':'source_eligible_records',
                  'baseline_eligible_outside_bounds':'source_eligible_outside_bounds'}
# Only engine provenance can vary with the reader available locally. No scientific
# value, query key, status, exclusion reason, population or denominator is ignored.
IGNORE_COLUMNS = {'engine'}
JOBS = {
 'coverage': [('01_realignment_evidence.xlsx',None),('02_fixed_marker_cohorts.xlsx',None),
              ('03_fixed_query_sets.xlsx',None),('04_paired_coverage_effects.xlsx',None),
              ('05_count_matched_sensitivity.xlsx',None),('06_matched_gene_background.xlsx',None)],
 'ld': [('01_genotype_coordinate_qc.xlsx',None),('02_ld_distance_profiles.xlsx',None),
        ('03_pairwise_dosage_ld.xlsx',None),('04_shared_collection_comparison.xlsx',None),
        ('05_observed_snp_neighbors.xlsx',None)]}


def normalize_name(key):return COLUMN_ALIASES.get(key,key)

def normalized_row(row):return {normalize_name(k):v for k,v in row.items() if k not in IGNORE_COLUMNS}


def same_value(a,b,tolerance):
    if a is None or b is None:return a is b
    if isinstance(a,bool) or isinstance(b,bool):return type(a) is type(b) and a==b
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        if not math.isfinite(a) or not math.isfinite(b):return False
        if isinstance(a,int) and isinstance(b,int):return a==b
        return abs(a-b)<=tolerance
    # Public terminology replaces only development-era nouns, not biological labels.
    if isinstance(a,str) and isinstance(b,str):
        b=b.replace('baseline','source').replace('v2.3.0','the exact-only policy')
    return a==b


def verify_hash_manifest(root):
    path=root/'OUTPUT_SHA256.json'
    obj=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj,dict) or not obj:raise DataError('Empty external reference manifest.')
    for name,digest in obj.items():
        rel=Path(name)
        if rel.is_absolute() or '..' in rel.parts or ':' in name or '\\' in name:raise DataError('Unsafe manifest path.')
        p=root/rel
        if not p.is_file() or p.is_symlink() or sha256(p)!=digest:raise DataError('External reference file missing/modified: '+name)
    return obj


def compare_workbooks(actual,reference,tolerance=1e-12):
    results=[];differences=[]
    with TableBook(actual) as a,TableBook(reference) as b:
        sa=[x for x in a.sheet_names if x!='README'];sb=[x for x in b.sheet_names if x!='README']
        if set(sa)!=set(sb):raise DataError('Worksheet inventory differs: '+actual.name)
        for sheet in sb:
            count=0;cells=0;mismatch=0;max_delta=0.0
            for count,pair in enumerate(zip_longest(a.table(sheet),b.table(sheet)),1):
                aa,bb=pair
                if aa is None or bb is None:
                    mismatch+=1
                    if len(differences)<200:differences.append(dict(file=actual.name,sheet=sheet,row=count+1,column='ROW_INVENTORY',actual='missing' if aa is None else 'present',expected='missing' if bb is None else 'present'))
                    continue
                ar,br=normalized_row(aa[1]),normalized_row(bb[1])
                if set(ar)!=set(br):raise DataError('Column inventory differs: '+actual.name+':'+sheet+' '+str(set(ar)^set(br)))
                for key in br:
                    av,bv=ar[key],br[key];cells+=1
                    if isinstance(av,(int,float)) and isinstance(bv,(int,float)) and not isinstance(av,bool) and not isinstance(bv,bool):max_delta=max(max_delta,abs(av-bv))
                    if not same_value(av,bv,tolerance):
                        mismatch+=1
                        if len(differences)<200:differences.append(dict(file=actual.name,sheet=sheet,row=count+1,column=key,actual=str(av),expected=str(bv)))
            results.append(dict(file=actual.name,sheet=sheet,rows=count,cells=cells,mismatched_cells=mismatch,maximum_numeric_difference=max_delta,status='PASS' if mismatch==0 else 'FAIL'))
        if a.errors or b.errors:raise DataError('Spreadsheet error values in a comparison input.')
    return results,differences


def verify_results(run,coverage=None,ld=None):
    run=run.resolve();verify_hash_manifest(run)
    if not coverage and not ld:raise DataError('Provide a reference coverage directory and/or a reference LD directory.')
    target=run/'validation';target.mkdir(exist_ok=False)
    rows=[];diff=[];inputs=[]
    for label,reference,sub in [('coverage',coverage,'alignment_coverage'),('ld',ld,'dosage_ld')]:
        if reference is None:continue
        reference=reference.resolve();manifest=verify_hash_manifest(reference)
        for file,_ in JOBS[label]:
            if file not in manifest:raise DataError('Comparison input not in reference manifest: '+file)
            ar=run/sub/file;br=reference/file
            results,details=compare_workbooks(ar,br)
            rows.extend(dict(analysis=label,**r) for r in results);diff.extend(dict(analysis=label,**r) for r in details)
            inputs.append(dict(analysis=label,filename=file,reference_sha256=sha256(br),new_sha256=sha256(ar)))
    failed=sum(r['mismatched_cells'] for r in rows)
    status='PASS' if not failed else 'FAIL'
    write_workbook(target/'numerical_agreement.xlsx',dict(sheet_checks=rows,differences=diff,comparison_inputs=inputs),[
        ('Status',status),('Tolerance','Integers, record keys, counts and strings must match; absolute tolerance 1e-12 for finite floating values.'),
        ('Coverage','Every non-README sheet of the declared numerical output workbooks is compared. All generated records are included; no preview is used.'),
        ('Alignment','No aligner is rerun here. Agreement establishes equivalent derived results for these supplied references.'),
        ('Reader provenance','The engine column is excluded because xlrd and the supported BIFF8 reader may differ; all other fields are checked.')])
    data=dict(status=status,worksheets=len(rows),rows=sum(r['rows'] for r in rows),cells=sum(r['cells'] for r in rows),mismatched_cells=failed,reference_inputs=inputs)
    (target/'agreement.json').write_text(json.dumps(data,indent=2)+'\n', encoding='utf-8')
    return data
