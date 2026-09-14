"""Read Excel source values without modifying files or evaluating formulas."""
from __future__ import annotations
import posixpath
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile
from ..evidence.common import DataError, HeaderNotFound, key, text
N = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'

class TableBook:
    """Read values without modifying workbooks or evaluating source formulas."""

    def __init__(self, path: Path):
        self.path = path
        self.zip = None
        self.xls = None
        self.biff = None
        self.errors = []
        if path.suffix.lower() == '.xls':
            try:
                import xlrd
            except ImportError:
                from .biff8 import BiffBook
                self.biff = BiffBook(path)
                self.sheet_names = list(self.biff.sheets)
                self.engine = 'restricted_biff8'
            else:
                self.xls = xlrd.open_workbook(str(path), on_demand=True)
                self.sheet_names = self.xls.sheet_names()
                self.engine = 'xlrd'
            return
        self.zip = ZipFile(path)
        rels = {r.attrib['Id']: r.attrib['Target'] for r in ET.fromstring(self.zip.read('xl/_rels/workbook.xml.rels'))}
        self.sheets = {}
        for s in ET.fromstring(self.zip.read('xl/workbook.xml')).find(N + 'sheets'):
            t = rels[s.attrib[REL]]
            target = t.lstrip('/') if t.startswith('/') else posixpath.normpath(posixpath.join('xl', t))
            if target.startswith('../'):
                raise DataError('Invalid worksheet relationship')
            self.sheets[s.attrib['name']] = target
        self.sheet_names = list(self.sheets)
        self.strings = []
        self.engine = 'ooxml_values'
        if 'xl/sharedStrings.xml' in self.zip.namelist():
            with self.zip.open('xl/sharedStrings.xml') as fh:
                it = ET.iterparse(fh, events=('start', 'end'))
                _, root = next(it)
                for ev, e in it:
                    if ev == 'end' and e.tag == N + 'si':
                        self.strings.append(''.join((t.text or '' for t in e.iter(N + 't'))))
                        root.clear()

    def close(self):
        if self.zip:
            self.zip.close()
        if self.xls:
            self.xls.release_resources()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def rows(self, sheet):
        if self.xls:
            s = self.xls.sheet_by_name(sheet)
            for i in range(s.nrows):
                values = s.row_values(i)
                for j, ct in enumerate(s.row_types(i)):
                    if ct == 5:
                        import xlrd
                        values[j] = xlrd.error_text_from_code.get(values[j], '#UNKNOWN!')
                        self.errors.append(dict(source_file=self.path.name,
                            source_sheet=sheet, source_cell=f'row{i + 1}_col{j + 1}', error=values[j]))
                yield (i + 1, values)
            return
        if self.biff:
            yield from self.biff.rows(sheet)
            return
        with self.zip.open(self.sheets[sheet]) as fh:
            it = ET.iterparse(fh, events=('start', 'end'))
            _, root = next(it)
            sheet_data = None
            for ev, e in it:
                if ev == 'start' and e.tag == N + 'sheetData':
                    sheet_data = e
                if ev != 'end' or e.tag != N + 'row':
                    continue
                data = {}
                for c in e:
                    if c.tag != N + 'c':
                        continue
                    letters = re.match('[A-Z]+', c.attrib['r']).group()
                    k = 0
                    for l in letters:
                        k = 26 * k + ord(l) - 64
                    typ = c.attrib.get('t')
                    v = c.find(N + 'v')
                    value = None
                    if c.find(N + 'f') is not None and (v is None or v.text is None):
                        raise DataError(f"Formula has no cached value in {self.path.name}:{sheet}:{c.attrib['r']}")
                    if typ == 'inlineStr':
                        value = ''.join((t.text or '' for t in c.iter(N + 't')))
                    elif v is not None:
                        if typ == 's':
                            value = self.strings[int(v.text)]
                        elif typ == 'b':
                            value = bool(int(v.text))
                        elif typ == 'e':
                            value = v.text
                            self.errors.append(dict(source_file=self.path.name,
                                source_sheet=sheet, source_cell=c.attrib['r'], error=value))
                        elif typ in ('str', 'd'):
                            value = v.text
                        else:
                            try:
                                value = float(v.text) if any((x in v.text for x in '.eE')) else int(v.text)
                            except (TypeError, ValueError):
                                value = v.text
                    data[k - 1] = value
                yield (int(e.attrib['r']), [data.get(i) for i in range(max(data, default=-1) + 1)])
                e.clear()
                if sheet_data is not None:
                    sheet_data.clear()

    def table(self, sheet, required=(), search_rows=12):
        rows = iter(self.rows(sheet))
        header = None
        for rownum, row in rows:
            normalized = [key(c) for c in row]
            if all((key(c) in normalized for c in required)) and any(normalized):
                header = [text(c).replace('\n', ' ') for c in row]
                break
            if rownum >= search_rows:
                break
        if header is None:
            raise HeaderNotFound(f'Expected header {list(required)} not found: {self.path.name}:{sheet}')
        for rownum, values in rows:
            if not any((text(v) for v in values)):
                continue
            yield (rownum, {h: values[i] if i < len(values) else None for i, h in enumerate(header) if h})

def field(row, *names):
    lookup = {key(k): v for k, v in row.items()}
    for name in names:
        if key(name) in lookup:
            return lookup[key(name)]
    return None
