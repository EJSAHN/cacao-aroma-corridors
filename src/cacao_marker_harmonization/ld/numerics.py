"""Pairwise-complete unphased dosage correlations and fixed, bounded pair inventories."""
from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
import numpy as np
from .common import DataError


def rng_for(seed: int, label: str):
    digest = hashlib.sha256((str(seed) + ':' + label).encode('utf-8')).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], 'little'))


def array_digest(a: np.ndarray) -> str:
    a = np.asarray(a, dtype='<i8')
    h = hashlib.sha256(json.dumps(list(a.shape), separators=(',', ':')).encode())
    h.update(a.tobytes(order='C'))
    return h.hexdigest()


@dataclass
class LDValues:
    n: np.ndarray
    r2: np.ndarray
    status: np.ndarray  # 0 valid; 1 too few; 2 zero variance on pairwise complete calls


def pairwise_ld(dosages: np.ndarray, pairs: np.ndarray, min_n: int, batch: int = 2048) -> LDValues:
    """Squared Pearson correlation on the same nonmissing samples in each pair.

    The input is marker x sample, values 0/1/2/NaN. No imputation, phasing,
    linkage model, Hardy-Weinberg assumption, or population adjustment is fitted.
    """
    d = np.asarray(dosages)
    p = np.asarray(pairs, dtype=np.int64)
    if d.ndim != 2 or p.ndim != 2 or p.shape[1] != 2:
        raise DataError('Invalid dosage/pair array dimensions.')
    if min_n < 3 or batch < 1:
        raise DataError('Invalid minimum pair samples or batch size.')
    if np.any(np.isinf(d)) or not np.all(np.isin(d[np.isfinite(d)], [0, 1, 2])):
        raise DataError('Dosages must be 0, 1, 2 or missing.')
    if len(p) and (p.min() < 0 or p.max() >= len(d) or np.any(p[:, 0] >= p[:, 1])):
        raise DataError('Pair indices must be unique nonself ordered indices.')
    if len(p) and len(np.unique(p[:,0] * len(d) + p[:,1])) != len(p):
        raise DataError('Duplicate genotype pair indices.')
    n = np.zeros(len(p), dtype=np.int32)
    out = np.full(len(p), np.nan, dtype=np.float64)
    status = np.ones(len(p), dtype=np.uint8)
    for lo in range(0, len(p), batch):
        hi = min(lo + batch, len(p))
        x = d[p[lo:hi, 0]].astype(np.float64)
        y = d[p[lo:hi, 1]].astype(np.float64)
        mask = np.isfinite(x) & np.isfinite(y)
        x = np.where(mask, x, 0.0); y = np.where(mask, y, 0.0)
        nn = mask.sum(axis=1).astype(np.int64)
        sx, sy = x.sum(axis=1), y.sum(axis=1)
        # Integer sufficient statistics avoid cancellation from missing-value means.
        cov = nn * (x * y).sum(axis=1) - sx * sy
        vx = nn * (x * x).sum(axis=1) - sx * sx
        vy = nn * (y * y).sum(axis=1) - sy * sy
        eligible = nn >= min_n
        good = eligible & (vx > 0) & (vy > 0)
        rr = np.full(hi - lo, np.nan)
        rr[good] = cov[good] ** 2 / (vx[good] * vy[good])
        if np.any(rr[good] < -1e-12) or np.any(rr[good] > 1 + 1e-12):
            raise DataError('Squared correlation outside numerical bounds.')
        rr[good] = np.clip(rr[good], 0.0, 1.0)
        n[lo:hi] = nn
        out[lo:hi] = rr
        status[lo:hi] = np.where(~eligible, 1, np.where(good, 0, 2))
    return LDValues(n, out, status)


def near_pair_codes(chroms, positions, max_bp: int, max_pairs: int) -> np.ndarray:
    """Enumerate every same-chromosome pair at <= max_bp; no variant-count cap."""
    chroms = np.asarray(chroms, dtype=str); positions = np.asarray(positions, dtype=np.int64)
    n = len(positions)
    if len(chroms) != n or max_bp < 0 or np.any(positions < 1):
        raise DataError('Invalid coordinates for pair enumeration.')
    chunks = []; total = 0
    for chrom in sorted(set(chroms)):
        ids = np.flatnonzero(chroms == chrom)
        ids = ids[np.lexsort((ids, positions[ids]))]
        pos = positions[ids]
        ends = np.searchsorted(pos, pos + max_bp, side='right')
        count = int(np.maximum(0, ends - np.arange(len(ids)) - 1).sum())
        total += count
        if total > max_pairs:
            raise DataError(f'Pair inventory exceeds configured safety limit {max_pairs}; nothing is silently sampled.')
        if count == 0:
            continue
        code = np.empty(count, dtype=np.int64); cursor = 0
        for k, end in enumerate(ends):
            others = ids[k + 1:int(end)]
            m = len(others)
            if m:
                left = np.minimum(ids[k], others); right = np.maximum(ids[k], others)
                code[cursor:cursor + m] = left * n + right; cursor += m
        chunks.append(code)
    if not chunks:
        return np.empty(0, dtype=np.int64)
    return np.sort(np.concatenate(chunks))


def sample_context_codes(chroms, positions, count: int, min_far_bp: int,
                         seed: int, label: str, which: str):
    """Uniform without-replacement pair ranks, independent of genotype values.

    Context is descriptive: interchromosomal pairs or distant same-chromosome
    pairs. These are not a null distribution for formal inference.
    """
    chroms = np.asarray(chroms, dtype=str); positions = np.asarray(positions, dtype=np.int64)
    n = len(positions)
    if which not in ('interchromosomal', 'distant_same_chromosome'):
        raise DataError('Unknown context class.')
    segments = []; sizes = []
    groups = {c: np.flatnonzero(chroms == c) for c in sorted(set(chroms))}
    if which == 'interchromosomal':
        cc = list(groups)
        for i, c in enumerate(cc):
            for e in cc[i + 1:]:
                a, b = groups[c], groups[e]
                segments.append((a, b, None)); sizes.append(len(a) * len(b))
    else:
        for c, ids in groups.items():
            ids = ids[np.lexsort((ids, positions[ids]))]
            pp = positions[ids]
            for k in range(len(ids)):
                start = max(k + 1, int(np.searchsorted(pp, pp[k] + min_far_bp, side='left')))
                if start < len(ids):
                    segments.append((ids, None, (k, start))); sizes.append(len(ids) - start)
    total = sum(sizes)
    size = min(count, total)
    if not size:
        return np.empty(0, dtype=np.int64), total
    ranks = np.sort(rng_for(seed, label + ':' + which).choice(total, size=size, replace=False))
    cum = np.cumsum(sizes, dtype=np.int64)
    seg = np.searchsorted(cum, ranks, side='right')
    prev = np.where(seg > 0, cum[np.maximum(seg - 1, 0)], 0)
    codes = []
    for ss, rank in zip(seg, ranks - prev):
        a, b, info = segments[int(ss)]
        if info is None:
            i, j = int(a[int(rank) // len(b)]), int(b[int(rank) % len(b)])
        else:
            k, start = info; i, j = int(a[k]), int(a[start + int(rank)])
        codes.append(min(i, j) * n + max(i, j))
    result = np.sort(np.array(codes, dtype=np.int64))
    if len(set(result)) != len(result):
        raise DataError('Duplicate context pair rank.')
    return result, total


def pair_inventory(old_c, old_p, new_c, new_p, cfg: dict, seed: int, label: str):
    n = len(new_p)
    old = near_pair_codes(old_c, old_p, cfg['max_distance_bp'], cfg['max_union_pairs'])
    new = near_pair_codes(new_c, new_p, cfg['max_distance_bp'], cfg['max_union_pairs'])
    ic, ni = sample_context_codes(new_c, new_p, cfg['context_pairs_per_class'], cfg['distant_same_chromosome_min_bp'], seed, label, 'interchromosomal')
    far, nf = sample_context_codes(new_c, new_p, cfg['context_pairs_per_class'], cfg['distant_same_chromosome_min_bp'], seed, label, 'distant_same_chromosome')
    codes = np.unique(np.concatenate([old, new, ic, far]))
    if len(codes) > cfg['max_union_pairs']:
        raise DataError('Combined pair inventory exceeds safety limit; increase explicitly or use a smaller declared distance.')
    pairs = np.column_stack((codes // n, codes % n)).astype(np.int32)
    flags = {k: np.isin(codes, a, assume_unique=True) for k, a in [('original_local', old), ('corrected_local', new), ('interchromosomal_context', ic), ('distant_context', far)]}
    return pairs, flags, dict(original_local_pairs=len(old), corrected_local_pairs=len(new),
                             interchromosomal_population_pairs=ni, interchromosomal_sampled_pairs=len(ic),
                             distant_population_pairs=nf, distant_sampled_pairs=len(far),
                             union_pairs=len(pairs), pair_index_sha256=array_digest(pairs))


def summary_statistics(values: LDValues, mask, thresholds):
    rr = values.r2[mask]
    valid = np.isfinite(rr); v = rr[valid]
    nn = values.n[mask]
    result = dict(pair_count=int(len(rr)), valid_pairs=int(len(v)),
                  insufficient_samples=int(np.sum(values.status[mask] == 1)),
                  zero_pair_variance=int(np.sum(values.status[mask] == 2)),
                  mean_r2=float(v.mean()) if len(v) else None,
                  median_r2=float(np.median(v)) if len(v) else None,
                  q90_r2=float(np.quantile(v, 0.9)) if len(v) else None,
                  q95_r2=float(np.quantile(v, 0.95)) if len(v) else None,
                  min_complete_samples=int(nn.min()) if len(nn) else None,
                  max_complete_samples=int(nn.max()) if len(nn) else None)
    for t in thresholds:
        suffix = str(t).replace('.', '_')
        result['pairs_r2_ge_' + suffix] = int(np.sum(v >= t))
        result['fraction_r2_ge_' + suffix] = float(np.mean(v >= t)) if len(v) else None
    return result
