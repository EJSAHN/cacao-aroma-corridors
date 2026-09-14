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

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()

def contains(parent: Path, child: Path) -> bool:
    parent, child = parent.resolve(), child.resolve()
    return parent == child or parent in child.parents

