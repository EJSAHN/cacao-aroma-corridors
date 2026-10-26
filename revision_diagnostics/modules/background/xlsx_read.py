"""Streaming OOXML values and worksheet metadata. Never evaluates formulas."""
from __future__ import annotations
from pathlib import Path
import posixpath, re, zipfile
from xml.etree import ElementTree as ET
N='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
R='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
class DataError(ValueError): pass
class Book:
    def __init__(self,path):
        self.path=Path(path); self.z=zipfile.ZipFile(self.path)
        rels={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(self.z.read('xl/_rels/workbook.xml.rels'))}
        self.sheets={}
        for s in ET.fromstring(self.z.read('xl/workbook.xml')).find(N+'sheets'):
            t=rels[s.attrib[R]]; target=t.lstrip('/') if t.startswith('/') else posixpath.normpath(posixpath.join('xl',t))
            if '..' in target.split('/'):raise DataError('Unsafe worksheet relationship')
            self.sheets[s.attrib['name']]=target
        self.strings=[]
        if 'xl/sharedStrings.xml' in self.z.namelist():
            for e in ET.fromstring(self.z.read('xl/sharedStrings.xml')).iter(N+'si'):
                self.strings.append(''.join(t.text or '' for t in e.iter(N+'t')))
    def __enter__(self):return self
    def __exit__(self,*args):self.z.close()
    def shape(self,sheet):
        with self.z.open(self.sheets[sheet]) as f:
            for _,e in ET.iterparse(f,events=('start',)):
                if e.tag==N+'dimension':
                    end=e.attrib['ref'].split(':')[-1]; m=re.fullmatch(r'([A-Z]+)([0-9]+)',end)
                    if not m:raise DataError('Unsupported worksheet dimension')
                    c=0
                    for t in m[1]:c=c*26+ord(t)-64
                    return int(m[2]),c
                if e.tag==N+'sheetData':break
        raise DataError('Worksheet must declare a dimension')
    def rows(self,sheet):
        with self.z.open(self.sheets[sheet]) as f:
            sd=None
            for ev,e in ET.iterparse(f,events=('start','end')):
                if ev=='start' and e.tag==N+'sheetData':sd=e
                if ev!='end' or e.tag!=N+'row':continue
                vals={}
                for c in e:
                    if c.tag!=N+'c':continue
                    k=0
                    for t in re.match('[A-Z]+',c.attrib['r'])[0]:k=k*26+ord(t)-64
                    typ=c.attrib.get('t'); v=c.find(N+'v')
                    if c.find(N+'f') is not None:raise DataError('Analytical inputs must contain values only: '+sheet)
                    if typ=='inlineStr':v=''.join(t.text or '' for t in c.iter(N+'t'))
                    elif v is None:v=None
                    elif typ=='s':v=self.strings[int(v.text)]
                    elif typ=='b':v=bool(int(v.text))
                    elif typ=='e':raise DataError('Excel error cell in '+sheet+':'+c.attrib['r'])
                    elif typ in ('str','d'):v=v.text
                    else:
                        t=v.text; v=float(t) if any(x in t for x in '.eE') else int(t)
                    vals[k-1]=v
                yield int(e.attrib['r']),[vals.get(i) for i in range(max(vals,default=-1)+1)]
                e.clear()
                if sd is not None:sd.clear()
    def records(self,sheet):
        it=iter(self.rows(sheet)); _,header=next(it)
        if header==['No_records']:return
        if len(header)!=len(set(header)) or any(not isinstance(x,str) or not x for x in header):raise DataError('Invalid headers: '+sheet)
        for row,values in it:
            if len(values)>len(header):raise DataError('Row exceeds header width')
            yield row,{k:values[i] if i<len(values) else None for i,k in enumerate(header)}
