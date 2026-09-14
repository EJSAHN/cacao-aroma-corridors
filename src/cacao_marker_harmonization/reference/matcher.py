"""Exhaustive exact full-tag matching in both orientations using NumPy seed lookup.

One SNP-free seed selects candidates. Both source alleles are then compared over
all tag bases, without mismatches, indels, soft clipping or reference-position
priors. Every contig is searched; reported coordinates are 1-based inclusive.
Uniqueness always means uniqueness among exact full-tag placements only.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import hashlib
import math
import re
from typing import Callable, Iterable
import numpy as np
from .common import DataError, text

_TRANS = bytes.maketrans(b'ACGT', b'TGCA')
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b'ACGT'):
    _LUT[_b] = _i


def reverse_complement(seq: bytes) -> bytes:
    return seq.translate(_TRANS)[::-1]


def clean_dna(value: object, allow_empty: bool = False) -> bytes | None:
    s = re.sub(r'\s+', '', text(value)).upper()
    if (not s and not allow_empty) or re.fullmatch('[ACGT]*', s) is None:
        return None
    return s.encode('ascii')


def encode_seed(seed: bytes) -> int:
    code = 0
    for b in seed:
        v = int(_LUT[b])
        if v == 255:
            raise DataError('A seed contains a non-ACGT base.')
        code = (code << 2) | v
    return code


def select_seed(ref: bytes, snp_index: int, k: int) -> tuple[int, bytes]:
    """Choose the highest-entropy SNP-free k-mer; deterministic position tie-break."""
    choices = []
    for off in range(len(ref) - k + 1):
        if off <= snp_index < off + k:
            continue
        seq = ref[off:off + k]
        counts = [seq.count(b) for b in b'ACGT']
        entropy = -sum((c / k) * math.log2(c / k) for c in counts if c)
        pairs = len(set(seq[j:j + 2] for j in range(k - 1)))
        choices.append((entropy, pairs, -abs(off + k / 2 - len(ref) / 2), -off, seq))
    if not choices:
        raise DataError('No SNP-free seed fits within this full tag.')
    best = max(choices)
    return -best[3], best[4]


@dataclass(frozen=True)
class Query:
    query_id: str
    source_row: int
    reference_key: str
    ref_tag: bytes
    alt_tag: bytes
    snp_index: int


@dataclass(frozen=True)
class Oriented:
    query_id: str
    ref_tag: bytes
    alt_tag: bytes
    snp_index: int
    strand: str
    seed_offset: int


@dataclass(frozen=True)
class Hit:
    query_id: str
    sequence_id: str
    tag_start_bp: int
    tag_end_bp: int
    snp_position_bp: int
    strand: str
    matched_source_allele: str
    assembly_base: str
    other_source_base: str

    def as_dict(self) -> dict:
        return dict(self.__dict__)


def build_queries(rows: list[dict], min_tag_length: int, seed_k: int):
    queries = []
    audit = []
    seen_rows = set()
    for row in rows:
        rn = int(row['source_row'])
        if rn in seen_rows:
            raise DataError('Duplicate source row in sequence evidence.')
        seen_rows.add(rn)
        qid = f'diversity_row_{rn}'
        left = clean_dna(row.get('five_flank'), True)
        right = clean_dna(row.get('three_flank'), True)
        variant = re.fullmatch(r'([ACGT])\s*>\s*([ACGT])', text(row.get('variation')).upper())
        status = 'ELIGIBLE'
        q = None
        if left is None or right is None or variant is None or variant[1] == variant[2]:
            status = 'UNSUPPORTED_SEQUENCE_OR_VARIATION'
        else:
            ref = left + variant[1].encode('ascii') + right
            alt = left + variant[2].encode('ascii') + right
            original = clean_dna(row.get('reference_sequence'))
            if original is None:
                status = 'UNSUPPORTED_REFERENCE_TAG'
            elif ref != original:
                status = 'FLANK_REFERENCE_TAG_MISMATCH'
            elif row.get('ref_lookup_status') != 'exact_ref_tag':
                status = 'PRIOR_REFERENCE_KEY_NOT_CONFIRMED'
            elif not text(row.get('reference_key')):
                status = 'MISSING_REFERENCE_KEY'
            elif len(ref) < min_tag_length:
                status = 'TAG_TOO_SHORT'
            else:
                try:
                    select_seed(ref, len(left), seed_k)
                except DataError:
                    status = 'NO_SNP_FREE_SEED'
                else:
                    q = Query(qid, rn, text(row['reference_key']), ref, alt, len(left))
        if q is not None:
            queries.append(q)
        audit.append(dict(query_id=qid, source_row=rn, reference_key=text(row.get('reference_key')),
            source_marker_id=text(row.get('marker_id')), source_chromosome=text(row.get('chromosome')),
            source_position_bp=row.get('position_bp'), source_variation=text(row.get('variation')),
            eligibility=status, tag_length=len(q.ref_tag) if q else None,
            source_snp_index_0based=q.snp_index if q else None,
            sequence_sha256=hashlib.sha256(q.ref_tag).hexdigest() if q else None))
    return queries, audit


class ExactTagMatcher:
    def __init__(self, queries: list[Query], seed_k: int = 12,
                 chunk_bases: int = 1_000_000, max_hits: int = 1_000_000):
        if not 4 <= seed_k <= 12:
            raise DataError('seed_k must be in [4,12] to bound dense seed memory.')
        if not 100 <= chunk_bases <= 10_000_000:
            raise DataError('chunk_bases must be in [100,10000000].')
        if max_hits < 1:
            raise DataError('max_hits must be positive.')
        if len({q.query_id for q in queries}) != len(queries):
            raise DataError('Duplicate query IDs.')
        self.k = seed_k
        self.chunk_bases = chunk_bases
        self.max_hits = max_hits
        self.buckets = defaultdict(list)
        self.present = np.zeros(4 ** seed_k, dtype=np.bool_)
        for q in queries:
            if len(q.ref_tag) != len(q.alt_tag) or not 0 <= q.snp_index < len(q.ref_tag):
                raise DataError('Malformed query lengths/index.')
            if q.ref_tag[q.snp_index] == q.alt_tag[q.snp_index]:
                raise DataError('Source alleles must differ.')
            if q.ref_tag[:q.snp_index] + q.ref_tag[q.snp_index + 1:] != q.alt_tag[:q.snp_index] + q.alt_tag[q.snp_index + 1:]:
                raise DataError('Queries must differ only at the reported SNP.')
            if clean_dna(q.ref_tag.decode('ascii')) != q.ref_tag or clean_dna(q.alt_tag.decode('ascii')) != q.alt_tag:
                raise DataError('Queries must be uppercase ACGT.')
            for strand in ('+', '-'):
                ref = q.ref_tag if strand == '+' else reverse_complement(q.ref_tag)
                alt = q.alt_tag if strand == '+' else reverse_complement(q.alt_tag)
                ix = q.snp_index if strand == '+' else len(ref) - 1 - q.snp_index
                off, seed = select_seed(ref, ix, seed_k)
                code = encode_seed(seed)
                self.present[code] = True
                self.buckets[code].append(Oriented(q.query_id, ref, alt, ix, strand, off))
        self.total_hits = 0
        self.total_bases = 0
        self.total_seed_candidates = 0

    def scan(self, sequence_id: str, sequence: bytes,
             progress: Callable[[str, int, int], None] | None = None) -> Iterable[Hit]:
        seq = sequence.upper()
        k = self.k
        last = len(seq) - k + 1
        if last <= 0:
            self.total_bases += len(seq)
            return
        for lo in range(0, last, self.chunk_bases):
            hi = min(last, lo + self.chunk_bases)
            n = hi - lo
            data = _LUT[np.frombuffer(seq[lo:hi + k - 1], dtype=np.uint8)]
            code = np.zeros(n, dtype=np.uint32)
            valid = np.ones(n, dtype=np.bool_)
            for j in range(k):
                part = data[j:j + n]
                code <<= 2
                code |= part & 3
                valid &= part < 4
            positions = np.flatnonzero(valid & self.present[code])
            self.total_seed_candidates += len(positions)
            for i in positions:
                seed_pos = lo + int(i)
                for q in self.buckets[int(code[i])]:
                    start = seed_pos - q.seed_offset
                    if start < 0 or start + len(q.ref_tag) > len(seq):
                        continue
                    if seq.startswith(q.ref_tag, start):
                        allele = 'source_ref'
                        base, other = q.ref_tag[q.snp_index], q.alt_tag[q.snp_index]
                    elif seq.startswith(q.alt_tag, start):
                        allele = 'source_alt'
                        base, other = q.alt_tag[q.snp_index], q.ref_tag[q.snp_index]
                    else:
                        continue
                    self.total_hits += 1
                    if self.total_hits > self.max_hits:
                        raise DataError('Exact-hit limit exceeded. No truncated completed report will be produced; increase max_hits only after reviewing repeat content.')
                    yield Hit(q.query_id, sequence_id, start + 1, start + len(q.ref_tag),
                              start + q.snp_index + 1, q.strand, allele, chr(base), chr(other))
            if progress is not None:
                progress(sequence_id, hi, max(0, last))
        self.total_bases += len(seq)


def summarize_mappings(audit: list[dict], hits: list[Hit], aliases: dict):
    indexed = defaultdict(list)
    for h in hits:
        indexed[h.query_id].append(h)
    result = []
    for row in audit:
        hs = indexed.get(row['query_id'], [])
        r = dict(row, exact_placements=len(hs), unique_exact_placement=False)
        if row['eligibility'] != 'ELIGIBLE':
            r['mapping_status'] = row['eligibility']
        elif not hs:
            r['mapping_status'] = 'NO_EXACT_FULL_TAG_MATCH'
        elif len(hs) != 1:
            r['mapping_status'] = 'MULTIPLE_EXACT_FULL_TAG_PLACEMENTS'
        else:
            h = hs[0]
            chrom = aliases.get(h.sequence_id)
            r.update(unique_exact_placement=True, mapping_status='UNIQUE_EXACT_PLACEMENT',
                     mapped_sequence_id=h.sequence_id, mapped_chromosome=chrom,
                     mapped_snp_bp=h.snp_position_bp, tag_start_bp=h.tag_start_bp,
                     tag_end_bp=h.tag_end_bp, mapped_strand=h.strand,
                     matched_source_allele=h.matched_source_allele, assembly_base=h.assembly_base,
                     other_source_base=h.other_source_base, primary_chromosome=chrom is not None)
            r['source_chromosome_agrees'] = row['source_chromosome'] == chrom if chrom is not None else None
            r['source_position_delta_bp'] = row['source_position_bp'] - h.snp_position_bp if r['source_chromosome_agrees'] and row.get('source_position_bp') is not None else None
        result.append(r)
    return result
