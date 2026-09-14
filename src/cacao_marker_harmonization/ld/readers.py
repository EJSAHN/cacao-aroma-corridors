"""Extract source facts; never repair coordinates, genotypes or annotations on read."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import sys
from .common import DataError, canonical_marker, chromosome, integer, key, text
from ..io.tableio import TableBook

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
