"""Shared parsing, hashing, deterministic random-number, and path-validation utilities."""
from __future__ import annotations
from ..errors import DataError, HeaderNotFound
import hashlib
import math
import re
import unicodedata
from pathlib import Path
from typing import Any
import numpy as np



def text(value: Any) -> str:
    if value is None:
        return ''
    if isinstance(value, (float, np.floating)) and (not math.isfinite(float(value))):
        return ''
    return str(value).strip()

def key(value: Any) -> str:
    return re.sub('[^a-z0-9]+',
        '', unicodedata.normalize('NFKD', text(value)).encode('ascii',
        'ignore').decode().lower())

def integer(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    s = re.sub('[\\s,\\u00a0]', '', text(value))
    if not s:
        return None
    try:
        n = float(s)
    except ValueError:
        return None
    if not math.isfinite(n) or not n.is_integer():
        return None
    return int(n)

def number(value: Any) -> float | None:
    try:
        n = float(text(value))
    except ValueError:
        return None
    return n if math.isfinite(n) else None

def chromosome(value: Any, aliases: dict[str, str] | None=None) -> str:
    """Use explicit labels; never extract a RefSeq version as a chromosome."""
    s = text(value)
    if not s:
        return ''
    if aliases and s in aliases:
        return text(aliases[s])
    m = re.search('(?:^|[_\\s])(?:chromosome|chr)[_\\s-]*(\\d+)(?:\\b|_)', s, re.I)
    if m:
        return str(int(m.group(1)))
    if re.fullmatch('\\d+(?:\\.0+)?', s):
        return str(int(float(s)))
    return s

def canonical_marker(value: Any) -> str:
    """Conservative normalization of known DArT spellings, not arbitrary punctuation removal."""
    s = text(value)
    m = re.fullmatch('(\\d+)\\|([FRfr])\\|(\\d+)--(\\d+)', s)
    if m:
        return f'{m[1]}|{m[2].upper()}|{m[3]}--{m[4]}'
    m = re.fullmatch('(\\d+)\\.([FRfr])\\.(\\d+)\\.(\\d+)', s)
    if m:
        return f'{m[1]}|{m[2].upper()}|{m[3]}--{m[4]}'
    return s

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def stable_id(prefix: str, values) -> str:
    s = '\x1f'.join((text(v) for v in values))
    return prefix + hashlib.sha256(s.encode('utf-8')).hexdigest()[:16]

def rng_for(seed: int, label: str) -> np.random.Generator:
    h = hashlib.sha256(f'{seed}:{label}'.encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], 'little'))

def within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False
