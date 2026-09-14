"""Streaming, values-only XLSX reports using the Python standard library."""
from __future__ import annotations
import math
import os
from pathlib import Path
import re
import tempfile
from xml.sax.saxutils import escape, quoteattr
import zipfile

NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
XML = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'


def column_name(n: int) -> str:
    if n < 1:
        raise ValueError('Column index must be positive.')
    out=''
    while n:
        n,k=divmod(n-1,26);out=chr(65+k)+out
    return out


def clean_value(v):
    if v is None:return None
    if isinstance(v,str):
        if len(v)>32767:raise ValueError('Excel string limit exceeded; no truncation.')
        if re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]',v):raise ValueError('Unsupported XML control character.')
        return v
    if isinstance(v,(bool,int)):return v
    if isinstance(v,float):
        if not math.isfinite(v):raise ValueError('Nonfinite result must be explicitly represented as unavailable.')
        return v
    raise TypeError('Unsupported Excel cell type: '+type(v).__name__)


def cell_xml(ref: str, value, style=0) -> str:
    v=clean_value(value)
    if v is None:return ''
    base=f'<c r="{ref}" s="{style}"'
    if isinstance(v,str):return base+' t="inlineStr"><is><t xml:space="preserve">'+escape(v)+'</t></is></c>'
    if isinstance(v,bool):return base+' t="b"><v>'+str(int(v))+'</v></c>'
    return base+'><v>'+repr(v)+'</v></c>'


def style_for(key, value):
    if isinstance(value,str):return 2
    if isinstance(value,float):
        if any(x in key for x in ('fraction','matched_q','background_q','delta_fraction','observed_minus_background')):
            return 4
        return 3
    return 0


STYLES=XML+f'''<styleSheet xmlns="{NS}">
<numFmts count="2"><numFmt numFmtId="164" formatCode="0.########"/><numFmt numFmtId="165" formatCode="0.00%"/></numFmts>
<fonts count="2"><font><sz val="10"/><name val="Calibri"/></font><font><b/><sz val="10"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font></fonts>
<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF234B5C"/><bgColor indexed="64"/></patternFill></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="5">
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyAlignment="1"><alignment wrapText="1" vertical="center"/></xf>
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf>
<xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
<xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
</cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>'''


def write_workbook(path: Path, sheets: dict[str,list[dict]], notes: list[tuple[str,str]]):
    if path.exists():raise FileExistsError('Refusing to overwrite: '+str(path))
    if 'README' in sheets:raise ValueError('README is reserved.')
    tables={'README':[dict(Field=k,Description=v) for k,v in notes],**sheets}
    if len(tables)!=len(set(s.casefold() for s in tables)):raise ValueError('Duplicate sheet name.')
    cols={}
    for name,rr in tables.items():
        if not name or len(name)>31 or re.search(r'[\[\]:*?/\\]',name):raise ValueError('Invalid sheet name.')
        if len(rr)+1>1048576:raise ValueError('Excel row limit exceeded; no preview output.')
        cc=list(rr.columns) if isinstance(rr, RowStream) else (list(dict.fromkeys(k for r in rr for k in r)) or ['No_records'])
        if len(cc)>16384:raise ValueError('Excel column limit exceeded.')
        for c in cc:clean_value(c)
        if not isinstance(rr, RowStream):
            for r in rr:
                for v in r.values():clean_value(v)
        cols[name]=cc
    fd,tmp=tempfile.mkstemp(prefix='.xlsx_',suffix='.tmp',dir=path.parent);os.close(fd)
    try:
        with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            types=XML+'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/>'
            types+='<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            types+=''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(tables)+1))+'</Types>'
            z.writestr('[Content_Types].xml',types)
            z.writestr('_rels/.rels',XML+'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
            wb=XML+f'<workbook xmlns="{NS}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><bookViews><workbookView/></bookViews><sheets>'
            wb+=''.join(f'<sheet name={quoteattr(name)} sheetId="{i}" r:id="rId{i}"/>' for i,name in enumerate(tables,1))+'</sheets></workbook>'
            z.writestr('xl/workbook.xml',wb)
            rels=XML+'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            rels+=''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(tables)+1))
            rels+=f'<Relationship Id="rId{len(tables)+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'
            z.writestr('xl/_rels/workbook.xml.rels',rels);z.writestr('xl/styles.xml',STYLES)
            for j,(name,rr) in enumerate(tables.items(),1):
                cc=cols[name];end=column_name(len(cc))+str(len(rr)+1)
                with z.open(f'xl/worksheets/sheet{j}.xml','w') as f:
                    def w(t): f.write(t.encode('utf-8'))
                    w(XML+f'<worksheet xmlns="{NS}"><dimension ref="A1:{end}"/><sheetViews><sheetView workbookViewId="0" showGridLines="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><sheetFormatPr defaultRowHeight="15"/><cols>')
                    for i,k in enumerate(cc,1):
                        width=min(40,max(12,len(k)+2,max((len(str(r.get(k,'')))+2 for r in rr[:60]),default=0)))
                        w(f'<col min="{i}" max="{i}" width="{width}" customWidth="1"/>')
                    w('</cols><sheetData><row r="1" ht="48" customHeight="1">')
                    w(''.join(cell_xml(column_name(i)+'1',k,1) for i,k in enumerate(cc,1)));w('</row>')
                    for i,r in enumerate(rr,2):
                        w(f'<row r="{i}">')
                        w(''.join(cell_xml(column_name(jj)+str(i),r.get(k),style_for(k,r.get(k))) for jj,k in enumerate(cc,1)))
                        w('</row>')
                    w('</sheetData>')
                    if rr:w(f'<autoFilter ref="A1:{end}"/>')
                    w('</worksheet>')
        if path.exists():raise FileExistsError('Output appeared during save.')
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


class RowStream:
    """Replayable values-only rows with explicit length and schema; no truncation."""
    def __init__(self, length, factory, columns):
        self.length = int(length)
        self.factory = factory
        self.columns = list(columns)
        if self.length < 0 or len(self.columns) != len(set(self.columns)):
            raise ValueError("Invalid streaming row length or columns.")
    def __len__(self):
        return self.length
    def __bool__(self):
        return bool(self.length)
    def __iter__(self):
        count=0
        for row in self.factory():
            if list(row)!=self.columns:
                raise ValueError("Streaming row schema changed.")
            count+=1
            if count>self.length:
                raise ValueError("Streaming row count exceeds declaration.")
            yield row
        if count!=self.length:
            raise ValueError("Streaming row count differs from declaration.")
    def __getitem__(self, index):
        import itertools
        if not isinstance(index,slice) or index.start not in (None,0) or index.step not in (None,1):
            raise TypeError("Only a bounded prefix is supported.")
        stop=min(self.length,index.stop or self.length)
        return list(itertools.islice(self.factory(),stop))
