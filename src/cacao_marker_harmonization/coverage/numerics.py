"""Closed-interval coverage with paired arms and bounded, seeded sampling."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import hashlib
import numpy as np
from .common import DataError

@dataclass
class Queries:
    name: str
    rows: list[dict]

    def __post_init__(self):
        self.ids = [r['query_id'] for r in self.rows]
        if len(self.ids) != len(set(self.ids)):
            raise DataError('Duplicate query ID in ' + self.name)
        self.chrom = np.array([str(r['chromosome']) for r in self.rows], dtype=object)
        self.start = np.array([r['start_bp'] for r in self.rows], dtype=np.int64)
        self.end = np.array([r['end_bp'] for r in self.rows], dtype=np.int64)
        if np.any(self.start < 1) or np.any(self.end < self.start):
            raise DataError('Invalid closed interval in ' + self.name)
        self.groups = {c: np.flatnonzero(self.chrom == c) for c in sorted(set(self.chrom))}

    def __len__(self):
        return len(self.rows)


def positions(rows: list[dict], chrom_key='chromosome', pos_key='position_bp') -> dict:
    result = defaultdict(list)
    for r in rows:
        result[str(r[chrom_key])].append(int(r[pos_key]))
    return {c: np.unique(np.array(v, dtype=np.int64)) for c, v in sorted(result.items())}


def nearest(q: Queries, points: dict) -> np.ndarray:
    """Zero inside the closed interval; infinity if its chromosome has no marker."""
    result = np.full(len(q), np.inf)
    for c, indices in q.groups.items():
        p = points.get(c)
        if p is None or not len(p):
            continue
        s, e = q.start[indices], q.end[indices]
        k = np.searchsorted(p, s, side='left')
        left = np.full(len(k), np.inf)
        right = np.full(len(k), np.inf)
        ok = k > 0
        left[ok] = s[ok] - p[k[ok] - 1]
        ok = k < len(p)
        right[ok] = np.maximum(p[k[ok]] - e[ok], 0)
        result[indices] = np.minimum(left, right)
    return result


def counts(distances: np.ndarray, windows: list[int]) -> np.ndarray:
    # O(n log n + w log n); no query-by-marker matrix.
    return np.searchsorted(np.sort(distances), np.asarray(windows), side='right').astype(np.int64)


def paired_counts(a: np.ndarray, b: np.ndarray, c: np.ndarray, window: int) -> dict:
    if not (a.shape == b.shape == c.shape):
        raise DataError('Paired arms have different query sets.')
    ha, hb, hc = a <= window, b <= window, c <= window
    if np.any(hb & ~ha):
        raise DataError('B coverage exceeds A despite B being a subset of A.')
    na, nb, nc = int(ha.sum()), int(hb.sum()), int(hc.sum())
    gain, loss = int((~hb & hc).sum()), int((hb & ~hc).sum())
    if gain - loss != nc - nb:
        raise DataError('Paired transition accounting failed.')
    return dict(n_queries=len(a), A_recovered=na, B_recovered=nb, C_recovered=nc,
                lost_on_selection=na - nb, gained_after_coordinate_change=gain,
                lost_after_coordinate_change=loss, both_B_C=int((hb & hc).sum()),
                neither_B_C=int((~hb & ~hc).sum()), selection_delta_count=nb - na,
                coordinate_delta_count=nc - nb, total_delta_count=nc - na)


def rng_for(seed: int, label: str):
    h = hashlib.sha256((str(seed) + '\x1f' + label).encode('utf-8')).digest()
    return np.random.Generator(np.random.PCG64(int.from_bytes(h[:16], 'big')))


def quotas_by_chromosome(panel_points: dict[str, dict], chromosomes: list[str]) -> dict:
    if not panel_points:
        raise DataError('No panels for count-matched sampling.')
    return {c: min(len(p.get(c, [])) for p in panel_points.values()) for c in chromosomes}


def sample_points(points: dict, quotas: dict, rng) -> dict:
    out = {}
    for chrom in sorted(quotas):
        n = quotas[chrom]
        p = points.get(chrom, np.array([], dtype=np.int64))
        if n < 0 or n > len(p):
            raise DataError('Invalid chromosome quota.')
        if n:
            out[chrom] = p.copy() if n == len(p) else np.sort(rng.choice(p, size=n, replace=False))
    return out


def gene_neighborhood_counts(genes: list[dict], radius: int) -> np.ndarray:
    """Other annotated gene midpoints within +/- radius; NOT marker density."""
    if radius < 0:
        raise DataError('Negative neighborhood radius.')
    out = np.zeros(len(genes), dtype=np.int64)
    group = defaultdict(list)
    for i, g in enumerate(genes):
        group[g['chromosome']].append(i)
    for indices in group.values():
        ii = np.array(indices)
        mid = np.array([(genes[i]['start_bp'] + genes[i]['end_bp']) // 2 for i in ii])
        sorted_mid = np.sort(mid)
        out[ii] = np.searchsorted(sorted_mid, mid + radius, side='right') - np.searchsorted(sorted_mid, mid - radius, side='left') - 1
    return out


def bin_cuts(values, nbins: int) -> np.ndarray:
    if nbins < 1 or len(values) == 0:
        raise DataError('Empty background or invalid bin count.')
    if nbins == 1:
        return np.array([], dtype=float)
    return np.unique(np.quantile(np.asarray(values), np.arange(1, nbins) / nbins, method='linear'))


def bin_of(value, cuts) -> int:
    return int(np.searchsorted(cuts, value, side='right'))


def sample_strata(pools: dict, required: dict, rng) -> np.ndarray:
    picks = []
    for s in sorted(required):
        n = required[s]
        pool = pools.get(s, np.array([], dtype=np.int64))
        if len(pool) < n:
            raise DataError('Insufficient background stratum; no pooling or replacement allowed.')
        if n:
            picks.extend(rng.choice(pool, size=n, replace=False).tolist())
    if len(picks) != len(set(picks)):
        raise DataError('A background draw duplicated a gene.')
    return np.array(picks, dtype=np.int64)
