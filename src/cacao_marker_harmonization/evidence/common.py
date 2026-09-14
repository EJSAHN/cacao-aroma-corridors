"""Strict source values, stable identifiers and bounded read-only input discovery."""
from __future__ import annotations
from ..errors import DataError, HeaderNotFound
import hashlib
import math
import re
import unicodedata
from pathlib import Path



def text(value) -> str:
    if value is None:
        return ''
    if isinstance(value, float) and not math.isfinite(value):
        return ''
    return str(value).strip()

def key(value) -> str:
    return re.sub('[^a-z0-9]+', '', unicodedata.normalize('NFKD', text(value)).encode('ascii','ignore').decode().lower())

def integer(value) -> int | None:
    if isinstance(value, bool):
        return None
    s = re.sub(r'[\s,\u00a0]', '', text(value))
    try:
        n = float(s)
    except (ValueError, OverflowError):
        return None
    return int(n) if math.isfinite(n) and n.is_integer() else None

def number(value) -> float | None:
    try:
        n = float(text(value))
    except (ValueError, OverflowError):
        return None
    return n if math.isfinite(n) else None

def chromosome(value, aliases=None) -> str:
    s = text(value)
    if aliases and s in aliases:
        return text(aliases[s])
    m = re.search(r'(?:^|[_\s])(?:chromosome|chr)[_\s-]*(\d+)(?:\b|_)', s, re.I)
    if m:
        return str(int(m[1]))
    if re.fullmatch(r'\d+(?:\.0+)?', s):
        return str(int(float(s)))
    return s

def canonical_marker(value) -> str:
    s = text(value)
    m = re.fullmatch(r'(\d+)\|([FRfr])\|(\d+)--(\d+)', s)
    if not m:
        m = re.fullmatch(r'(\d+)\.([FRfr])\.(\d+)\.(\d+)', s)
    return f'{m[1]}|{m[2].upper()}|{m[3]}--{m[4]}' if m else s

def dart_suffix(value) -> int | None:
    m = re.fullmatch(r'\d+\|[FR]\|\d+--(\d+)', canonical_marker(value))
    return int(m[1]) if m else None

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()

def stable_id(prefix: str, values) -> str:
    payload = '\x1f'.join(text(v) for v in values)
    return prefix + hashlib.sha256(payload.encode('utf-8')).hexdigest()[:20]

def contains(parent: Path, child: Path) -> bool:
    parent, child = parent.resolve(), child.resolve()
    return parent == child or parent in child.parents

def locate(root: Path, spec: dict) -> Path:
    relative = Path(spec['path'])
    if relative.is_absolute() or '..' in relative.parts:
        raise DataError('Source configuration paths must be relative to the source root.')
    expected = root / relative
    if expected.is_file():
        if not contains(root, expected):
            raise DataError('Source path resolves outside the input root.')
        return expected
    names = {n.casefold() for n in spec.get('aliases', [relative.name])}
    found = sorted(p for p in root.rglob('*') if p.is_file() and p.name.casefold() in names and 'original_archives' not in p.parts)
    if not found:
        raise DataError(f"Required source missing: {spec['path']}. Specify its actual relative path in evidence_config.json.")
    if any(not contains(root,p) for p in found):
        raise DataError('A source alias resolves outside the input root.')
    if len(found)>1 and len({sha256(p) for p in found})>1:
        raise DataError('Ambiguous source copies: ' + '; '.join(str(p.relative_to(root)) for p in found))
    return found[0]

def fraction(n: int, d: int) -> float | None:
    return n/d if d else None

def close_number(a, b, rtol=1e-9, atol=0.0) -> bool:
    return a is not None and b is not None and math.isclose(a,b,rel_tol=rtol,abs_tol=atol)
