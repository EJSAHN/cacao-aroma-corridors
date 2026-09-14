"""Restricted, read-only BIFF8 value reader.

Used only when xlrd is unavailable. Supports the OLE regular-stream layout and
SST/NUMBER/RK/MULRK/LABELSST/BOOLERR cells found in the supplied Nacional exports.
It does not evaluate formulas, decrypt files, repair files, or silently ignore
unsupported cell representations. Prefer xlrd for other historical XLS files.
"""
from __future__ import annotations
import struct
from pathlib import Path
from ..evidence.common import DataError

def u32(b, off=0):
    return struct.unpack_from('<I', b, off)[0]

def u16(b, off=0):
    return struct.unpack_from('<H', b, off)[0]

def stream(path: Path) -> bytes:
    b = path.read_bytes()
    if b[:8] != bytes.fromhex('d0cf11e0a1b11ae1'):
        raise DataError(f'Not an OLE XLS file: {path.name}')
    sector = 1 << u16(b, 30)
    if sector not in (512, 4096):
        raise DataError('Unsupported OLE sector size; use xlrd.')

    def get(i):
        result = b[(i + 1) * sector:(i + 2) * sector]
        if len(result) != sector:
            raise DataError('Truncated OLE sector')
        return result
    difat = list(struct.unpack_from('<109I', b, 76))
    nxt = u32(b, 68)
    for _ in range(u32(b, 72)):
        ints = list(struct.unpack(f'<{sector // 4}I', get(nxt)))
        difat.extend(ints[:-1])
        nxt = ints[-1]
    fat = []
    for i in difat[:u32(b, 44)]:
        fat.extend(struct.unpack(f'<{sector // 4}I', get(i)))

    def chain(i):
        chunks, seen = ([], set())
        while i < 4294967290:
            if i in seen or i >= len(fat):
                raise DataError('Invalid OLE allocation chain')
            seen.add(i)
            chunks.append(get(i))
            i = fat[i]
        return b''.join(chunks)
    directory = chain(u32(b, 48))
    for j in range(0, len(directory), 128):
        d = directory[j:j + 128]
        name = d[:max(u16(d, 64) - 2, 0)].decode('utf-16le')
        if name in ('Workbook', 'Book'):
            size = struct.unpack_from('<Q', d, 120)[0]
            if size < 4096:
                raise DataError('Mini-stream XLS requires xlrd.')
            return chain(u32(d, 116))[:size]
    raise DataError('No Workbook stream in XLS')

def records(b, offset=0):
    while offset + 4 <= len(b):
        typ, n = struct.unpack_from('<HH', b, offset)
        if offset + 4 + n > len(b):
            raise DataError('Truncated BIFF record')
        yield (typ, b[offset + 4:offset + 4 + n])
        offset += 4 + n

def decode_sst(chunks):

    class Cursor:
        k = 0
        i = 0

        def raw(self, n):
            out = bytearray()
            while n:
                if self.k >= len(chunks):
                    raise DataError('Truncated SST')
                if self.i == len(chunks[self.k]):
                    self.k += 1
                    self.i = 0
                if self.k >= len(chunks):
                    raise DataError('Truncated SST')
                take = min(n, len(chunks[self.k]) - self.i)
                out.extend(chunks[self.k][self.i:self.i + take])
                self.i += take
                n -= take
            return bytes(out)

        def string(self, n, wide):
            out = []
            while n:
                if self.i == len(chunks[self.k]):
                    self.k += 1
                    self.i = 0
                    if self.k >= len(chunks):
                        raise DataError('Truncated SST string')
                    wide = bool(self.raw(1)[0] & 1)
                take = min(n, (len(chunks[self.k]) - self.i) // (2 if wide else 1))
                if not take:
                    raise DataError('Unsupported SST code-unit boundary; use xlrd.')
                out.append(self.raw(take * (2 if wide else 1)).decode('utf-16le' if wide else 'latin1'))
                n -= take
            return ''.join(out)
    c = Cursor()
    _, count = struct.unpack('<II', c.raw(8))
    result = []
    for _ in range(count):
        n = u16(c.raw(2))
        flags = c.raw(1)[0]
        rich = u16(c.raw(2)) if flags & 8 else 0
        ext = u32(c.raw(4)) if flags & 4 else 0
        result.append(c.string(n, bool(flags & 1)))
        c.raw(4 * rich + ext)
    return result

class BiffBook:

    def __init__(self, path: Path):
        self.data = stream(path)
        self.sheets = {}
        self.strings = []
        chunks = []
        for t, v in records(self.data):
            if chunks and t != 60:
                self.strings = decode_sst(chunks)
                chunks = []
            if t == 47:
                raise DataError('Encrypted XLS is unsupported')
            if t == 133:
                n, wide = (v[6], v[7] & 1)
                name = v[8:8 + n * (2 if wide else 1)].decode('utf-16le' if wide else 'latin1')
                self.sheets[name] = u32(v)
            elif t == 252:
                chunks = [v]
            elif t == 60 and chunks:
                chunks.append(v)
            elif t == 10:
                break
        if chunks:
            self.strings = decode_sst(chunks)

    @staticmethod
    def rk(i):
        value = struct.unpack('<i',
            struct.pack('<I', i))[0] >> 2 if i & 2 else struct.unpack('<d',
            struct.pack('<II', 0, i & 4294967292))[0]
        return value / 100 if i & 1 else value

    def rows(self, sheet):
        previous, cells = (-1, {})
        for t, v in records(self.data, self.sheets[sheet]):
            if t == 10:
                break
            entries = []
            if t in (253, 515, 638):
                r, col, _ = struct.unpack_from('<HHH', v)
                val = self.strings[u32(v,
                    6)] if t == 253 else struct.unpack_from('<d', v, 6)[0] if t == 515 else self.rk(u32(v,
                    6))
                entries = [(r, col, val)]
            elif t == 189:
                r, col = struct.unpack_from('<HH', v)
                end = u16(v, len(v) - 2)
                entries = [(r, col + k, self.rk(u32(v, 6 + 6 * k))) for k in range(end - col + 1)]
            elif t == 517:
                r, col, _ = struct.unpack_from('<HHH', v)
                if v[7]:
                    raise DataError('Excel error cell in input')
                entries = [(r, col, bool(v[6]))]
            elif t in (6, 516, 214, 2, 3, 4, 5):
                raise DataError('Unsupported BIFF cell type or formula; use xlrd.')
            for r, col, val in entries:
                if r < previous:
                    raise DataError('Non-row-ordered BIFF cells require xlrd.')
                if r != previous and previous >= 0:
                    yield (previous + 1, [cells.get(i) for i in range(max(cells) + 1)])
                    cells = {}
                previous = r
                cells[col] = val
        if cells:
            yield (previous + 1, [cells.get(i) for i in range(max(cells) + 1)])
