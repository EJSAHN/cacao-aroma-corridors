"""Extract source facts; never repair coordinates, genotypes or annotations on read."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import json
import sys
from .common import DataError, canonical_marker, chromosome, integer, key, number, text, stable_id
from ..io.tableio import TableBook, field

@dataclass
class Matrix:
    name: str
    samples: list[str]
    markers: list[dict]
    calls: list[tuple[str, ...]]
    sheet: str
    engine: str
    errors: list[dict]

def read_matrix(path: Path, name: str, cfg: dict) -> Matrix:
    with TableBook(path) as book:
        configured = cfg.get('matrix_sheet')
        if configured:
            if configured not in book.sheet_names:
                raise DataError(f'Configured matrix sheet absent: {path.name}:{configured}')
            sheet=configured
        elif len(book.sheet_names)==1:
            sheet=book.sheet_names[0]
        else:
            raise DataError(f'Multiple sheets in {path.name}; configure matrix_sheet explicitly.')
        it=book.rows(sheet)
        try:
            _,header=next(it)
        except StopIteration:
            raise DataError(f'Empty genotype source: {path.name}')
        header=[text(v) for v in header]
        normalized=[key(v) for v in header]
        if len(normalized)!=len(set(normalized)):
            raise DataError(f'Duplicate/ambiguous header in {path.name}')
        indexes={h:i for i,h in enumerate(normalized)}
        mid=indexes.get('rs',indexes.get('markerid'))
        if mid is None or 'chrom' not in indexes or 'pos' not in indexes:
            raise DataError(f'Marker/chrom/pos header not found in {path.name}')
        meta={'rs','markerid','alleles','chrom','pos','strand','assembly','center','protlsid','assaylsid','panellsid','qccode'}
        sample_indexes=[i for i,h in enumerate(normalized) if h not in meta]
        samples=[header[i] for i in sample_indexes]
        if not samples or any(not s for s in samples) or len(set(samples))!=len(samples):
            raise DataError(f'Invalid sample identifiers in {path.name}')
        markers=[]; calls=[]
        for rownum,vals in it:
            if len(vals)>len(header) and any(text(v) for v in vals[len(header):]):
                raise DataError(f'Unexpected extra matrix cells: {path.name}:{rownum}')
            vals=vals+[None]*max(0,len(header)-len(vals))
            marker=text(vals[mid])
            if not marker:
                if any(text(v) for v in vals):
                    raise DataError(f'Nonempty genotype row lacks marker ID: {path.name}:{rownum}')
                continue
            rowcalls=tuple(sys.intern(text(vals[i]).upper()) for i in sample_indexes)
            if any(len(v)>16 for v in rowcalls):
                raise DataError(f'Non-genotype content in {path.name}:{rownum}')
            get=lambda k: text(vals[indexes[k]]) if k in indexes else ''
            markers.append(dict(resource=name,source_file=path.name,source_sheet=sheet,source_row=rownum,
                marker_id=marker,canonical_id=canonical_marker(marker),chromosome=chromosome(get('chrom'),cfg.get('chromosome_aliases')),
                position_bp=integer(get('pos')),chromosome_raw=get('chrom'),position_raw=get('pos'),
                alleles_reported=get('alleles'),strand_reported=get('strand'),record_index=len(markers)))
            calls.append(rowcalls)
        return Matrix(name,samples,markers,calls,sheet,book.engine,list(book.errors))

def read_diversity(path: Path, cfg: dict):
    records=[]; references=defaultdict(list); metadata=[]
    with TableBook(path) as book:
        for rownum,row in book.table(cfg['diversity_marker_sheet'],required=['SNP marker name']):
            marker=text(field(row,'SNP marker name'))
            if not marker:
                continue
            records.append(dict(source_file=path.name,source_sheet=cfg['diversity_marker_sheet'],source_row=rownum,
                marker_id=marker,marker_key=canonical_marker(marker),
                reference_sequence_name=text(field(row,'reference sequence name')),
                reference_key=canonical_marker(field(row,'reference sequence name')),
                chromosome=chromosome(field(row,'chromosome'),cfg.get('chromosome_aliases')),
                position_bp=integer(field(row,'snp position')),chromosome_raw=text(field(row,'chromosome')),
                position_raw=text(field(row,'snp position')),assembly_reported=text(field(row,'version')),
                five_flank=text(field(row,'5flank_sequence')),three_flank=text(field(row,'3flank_sequence')),
                variation=text(field(row,'variation')),strand_reported=text(field(row,'strand')),
                remark=text(field(row,'remark marker'))))
        for rownum,row in book.table(cfg['diversity_reference_sheet'],required=['reference sequence name']):
            name=text(field(row,'reference sequence name'))
            if name:
                references[canonical_marker(name)].append(dict(source_row=rownum,sequence=text(field(row,'Sequence')),
                    reference_name=name, accession=text(field(row,'accession number')),
                    database=text(field(row,'sequence database name'))))
        if 'study' in book.sheet_names:
            it=book.rows('study'); h=next(it,None);r=next(it,None)
            if h and r:
                metadata=[dict(field=text(k),value=text(v)) for k,v in zip(h[1],r[1]) if text(k)]
        return records,references,metadata,list(book.errors)

def read_associations(path: Path, cfg: dict):
    records=[];headers=[]
    with TableBook(path) as book:
        if cfg['association_master_sheet'] not in book.sheet_names:
            raise DataError('Association master sheet absent.')
        from .common import HeaderNotFound
        for sheet in book.sheet_names:
            try:
                rows=list(book.table(sheet,required=['Chromosome','Traits']))
            except HeaderNotFound:
                continue
            if rows:
                headers.append(dict(source_sheet=sheet,columns=json.dumps(list(rows[0][1]),ensure_ascii=False)))
            for rownum,row in rows:
                c=chromosome(field(row,'Chromosome'),cfg.get('chromosome_aliases'))
                p=integer(field(row,'Position of the association peak (bp)','Position of the association peak'))
                trait=text(field(row,'Traits'))
                if not c and p is None and not trait:continue
                r=dict(source_file=path.name,source_sheet=sheet,source_row=rownum,chromosome=c,peak_bp=p,trait=trait,
                    method_reported=text(field(row,'GWAS method')),marker_filter_reported=text(field(row,'Sorting of marker')),
                    source_start_bp=integer(field(row,'Position of haplotypic bloc start')),
                    source_end_bp=integer(field(row,'Position of haplotypic bloc end')),
                    source_block_number=integer(field(row,'N° haplotypic bloc','N° hap. Bloc')),
                    p_value_reported=number(field(row,'p-value of the strongest association')),
                    explained_variance_reported=number(field(row,'Explanation rate of the trait of the strongest association')),
                    raw_record=json.dumps(row,ensure_ascii=False,default=str))
                r['event_key']=stable_id('a_', (c,p,trait,r['method_reported'],r['marker_filter_reported']))
                r['in_master']=sheet==cfg['association_master_sheet']
                records.append(r)
        return records,headers,list(book.errors)

def read_candidate_metadata(path: Path, cfg: dict):
    # Count distinct source IDs without modifying annotation or using them as a background universe.
    counts=Counter();errors=[]
    with TableBook(path) as book:
        for _,row in book.table(cfg['candidate_sheet'],required=['gene_id','start','end']):
            g=text(field(row,'gene_id'))
            if g:counts[g]+=1
        errors=list(book.errors)
    return [dict(gene_id=g,source_records=n) for g,n in sorted(counts.items())],errors
