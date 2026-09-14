"""Pinned reference acquisition and streaming FASTA/GFF3 validation.

Downloaded bytes are never trusted by filename alone. Assembly headers, all
primary sequence aliases and the expected toplevel length are checked. SHA-256
values are recorded, not represented as publisher-verified checksums unless an
expected digest is explicitly provided. No TLS validation is disabled.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote, urlsplit
import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
import time
import urllib.error
import urllib.request
from .common import DataError, sha256


@dataclass
class ReferenceAsset:
    path: Path
    role: str
    origin: str
    url: str
    sha256: str
    bytes: int

    def report(self):
        return dict(role=self.role, filename=self.path.name, origin=self.origin,
                    source_url=self.url, sha256=self.sha256, bytes=self.bytes)


def open_text(path: Path):
    return gzip.open(path, 'rt', encoding='utf-8-sig', newline=None) if path.suffix.lower() == '.gz' else path.open('r', encoding='utf-8-sig', newline=None)


def iter_fasta(path: Path, max_contig_bases: int):
    """One contig in memory. A truncated gzip raises; duplicate sequence IDs fail."""
    seen = set()
    header = None
    sequence = bytearray()
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    if not sequence:
                        raise DataError('Empty FASTA sequence: ' + header)
                    yield header.split()[0], header, bytes(sequence)
                header = line[1:].strip()
                if not header or header.split()[0] in seen:
                    raise DataError(f'Empty or duplicate FASTA sequence name at line {lineno}.')
                seen.add(header.split()[0])
                sequence = bytearray()
            else:
                if header is None:
                    raise DataError('FASTA sequence appears before a header.')
                if not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', line):
                    raise DataError(f'Invalid FASTA sequence characters at line {lineno}.')
                sequence.extend(line.upper().encode('ascii'))
                if len(sequence) > max_contig_bases:
                    raise DataError('Reference contig exceeds the configured memory bound.')
    if header is None or not sequence:
        raise DataError('Missing or empty final FASTA sequence.')
    yield header.split()[0], header, bytes(sequence)


def validate_sequence_header(seqid: str, header: str, spec: dict):
    token = spec['assembly_name']
    if token not in header:
        raise DataError(f'Reference header lacks the required assembly name {token}: {seqid}. Supply the pinned Ensembl toplevel FASTA; a filename is not assembly verification.')
    # Validate length separately after sequence loading; no numeric accession parsing.


def validate_fasta_inventory(inventory: list[dict], spec: dict):
    if not inventory:
        raise DataError('No reference sequences scanned.')
    total = sum(r['length_bp'] for r in inventory)
    expected = int(spec['expected_total_bases'])
    if total != expected:
        raise DataError(f'Reference total length {total:,} differs from the pinned assembly expectation {expected:,}. Whole-genome exact-match uniqueness cannot be certified from a partial or different FASTA.')
    counts = {}
    for row in inventory:
        chrom = row['chromosome']
        if chrom is not None:
            counts[chrom] = counts.get(chrom, 0) + 1
    required = set(spec['chromosomes'])
    if set(counts) != required or any(v != 1 for v in counts.values()):
        raise DataError('FASTA must contain exactly one named sequence for each configured primary chromosome.')
    return total


def _prefix_validate(path: Path, role: str):
    try:
        with open_text(path) as f:
            lines = [f.readline() for _ in range(12)]
    except (OSError, EOFError, UnicodeError) as exc:
        raise DataError(f'Unreadable or truncated {role} reference: {path.name}: {exc}') from exc
    first = next((s.strip() for s in lines if s.strip()), '')
    if role == 'fasta' and not first.startswith('>'):
        raise DataError('Reference download is not FASTA (possibly an HTML error page).')
    if role == 'gff' and not any(s.startswith('##gff-version 3') for s in lines):
        raise DataError('Annotation download is not a declared GFF3 file.')


def download_asset(urls: list[str], dest: Path, role: str, max_bytes: int, log,
                   timeout: int = 60, attempts: int = 2) -> tuple[str, str]:
    """Download to the configured cache atomically; clean only newly created partials."""
    errors = []
    for url in urls:
        if urlsplit(url).scheme != 'https':
            raise DataError('Only HTTPS reference downloads are supported.')
        for attempt in range(attempts):
            tmp = None
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'CacaoMarkerHarmonization/3.0.0 (research reference download)'})
                log.info('Downloading %s from %s', role, url)
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    if urlsplit(response.geturl()).scheme != 'https':
                        raise DataError('Reference download redirected to a non-HTTPS URL.')
                    length = response.headers.get('Content-Length')
                    expected = int(length) if length and length.isdigit() else None
                    if expected is not None and expected > max_bytes:
                        raise DataError('Reference download exceeds configured byte limit.')
                    fd, name = tempfile.mkstemp(prefix='.download_', suffix=''.join(dest.suffixes), dir=dest.parent)
                    tmp = Path(name)
                    received = 0
                    digest = hashlib.sha256()
                    last_log = 0
                    with os.fdopen(fd, 'wb') as out:
                        while True:
                            data = response.read(1024 * 1024)
                            if not data:
                                break
                            received += len(data)
                            if received > max_bytes:
                                raise DataError('Reference download exceeded configured byte limit.')
                            out.write(data)
                            digest.update(data)
                            if received - last_log >= 20 * 1024 * 1024:
                                log.info('Downloaded %s: %.1f MiB', role, received / 2 ** 20)
                                last_log = received
                        out.flush()
                        os.fsync(out.fileno())
                    if expected is not None and received != expected:
                        raise DataError('Incomplete reference HTTP response.')
                _prefix_validate(tmp, role)
                if dest.exists():
                    raise DataError('Refusing to overwrite an existing reference cache file.')
                os.replace(tmp, dest)
                tmp = None
                return url, digest.hexdigest()
            except (OSError, urllib.error.URLError, ValueError) as exc:
                errors.append(f'{url}: {type(exc).__name__}: {exc}')
                log.warning('Download attempt failed: %s', errors[-1])
                if attempt + 1 < attempts:
                    time.sleep(1)
            finally:
                if tmp is not None and tmp.exists():
                    tmp.unlink()
    raise DataError('Reference retrieval failed. Existing data and prior results were not changed. Use --reference-dir for an existing pinned FASTA/GFF3 or retry when the connection is available.\n' + '\n'.join(errors))


def acquire_asset(role: str, spec: dict, cache: Path, search_dirs: list[Path],
                  explicit: Path | None, offline: bool, log) -> ReferenceAsset:
    conf = spec[role]
    filename = conf['filename']
    if Path(filename).name != filename:
        raise DataError('Reference asset filenames must be plain filenames.')
    candidates = []
    if explicit is not None:
        if not explicit.is_file():
            raise DataError('Explicit reference file is missing: ' + str(explicit))
        candidates = [explicit.resolve()]
    else:
        for root in [cache, *search_dirs]:
            # Bounded, exact-name lookup. No recursive scan of the user's drives.
            for name in (filename, filename[:-3] if filename.endswith('.gz') else filename):
                p = root / name
                if p.is_file() and p.resolve() not in candidates:
                    candidates.append(p.resolve())
    if candidates:
        p = candidates[0]
        _prefix_validate(p, role)
        digest = sha256(p)
        expected = conf.get('sha256')
        if expected and expected != digest:
            raise DataError('Reference file differs from the configured trusted SHA-256.')
        sidecar = cache / (filename + '.provenance.json')
        recorded_url = ''
        if p == (cache / filename).resolve() and sidecar.is_file():
            saved = json.loads(sidecar.read_text(encoding='utf-8'))
            if saved.get('sha256') != digest:
                raise DataError('Cached reference differs from the digest recorded when acquired; it was NOT overwritten.')
            recorded_url = saved.get('source_url', '')
        log.info('Reusing %s reference without copying: %s', role, p)
        return ReferenceAsset(p, role, 'existing_read_only_file', recorded_url, digest, p.stat().st_size)
    if offline:
        raise DataError(f'Offline mode: {filename} not found. Use --reference-dir or --fasta/--gff.')
    cache.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(cache).free < float(spec.get('minimum_cache_free_gib', 2)) * 2 ** 30:
        raise DataError('Insufficient space on the reference-cache drive.')
    p = cache / filename
    url, digest = download_asset(conf['urls'], p, role, int(conf['maximum_bytes']), log)
    expected = conf.get('sha256')
    if expected and expected != digest:
        raise DataError('Downloaded bytes differ from the configured SHA-256. File retained for inspection, not used.')
    asset = ReferenceAsset(p, role, 'downloaded_pinned_release', url, digest, p.stat().st_size)
    sidecar = cache / (filename + '.provenance.json')
    sidecar.write_text(json.dumps(asset.report(), indent=2) + '\n', encoding='utf-8')
    return asset


def parse_attributes(value: str) -> dict[str, str]:
    result = {}
    for token in value.split(';'):
        if not token:
            continue
        k, sep, v = token.partition('=')
        if sep:
            result[unquote(k)] = unquote(v)
    return result


def read_gff(path: Path, seq_lengths: dict, aliases: dict, spec: dict):
    genes = []
    headers = []
    errors = []
    seen = {}
    build_seen = False
    declared_gff3 = False
    with open_text(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line == '##FASTA':
                break
            if line.startswith('#'):
                if line.startswith('##gff-version 3'):
                    declared_gff3 = True
                if spec['assembly_name'] in line or spec['assembly_accession'] in line:
                    build_seen = True
                if len(headers) < 100:
                    headers.append(dict(line_number=n, text=line))
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                errors.append(dict(line_number=n, issue='invalid_column_count'))
                continue
            seqid, source, kind, start, end, score, strand, phase, attrs = fields
            if kind not in spec['gff_gene_features']:
                continue
            attr = parse_attributes(attrs)
            raw_id = attr.get('ID', '')
            gene_id = raw_id[5:] if raw_id.startswith('gene:') else raw_id
            try:
                s, e = int(start), int(end)
            except ValueError:
                errors.append(dict(line_number=n, issue='non_integer_gene_boundary', gene_id=gene_id))
                continue
            issue = ''
            if not gene_id:
                issue = 'missing_gene_id'
            elif seqid not in seq_lengths:
                issue = 'gene_sequence_absent_from_fasta'
            elif s < 1 or e < s or e > seq_lengths[seqid]:
                issue = 'gene_outside_reference_sequence'
            elif strand not in ('+', '-', '.', '?'):
                issue = 'invalid_gene_strand'
            if issue:
                errors.append(dict(line_number=n, issue=issue, gene_id=gene_id, sequence_id=seqid))
                continue
            gene = dict(gene_id=gene_id, original_gene_id=raw_id, gene_name=attr.get('Name', ''),
                        sequence_id=seqid, chromosome=aliases.get(seqid), start_bp=s, end_bp=e,
                        strand=strand, gene_feature=kind, biotype=attr.get('biotype', attr.get('gene_biotype', '')),
                        annotation_source=source, original_line=n,
                        description=attr.get('description', ''))
            if gene_id in seen:
                errors.append(dict(line_number=n, issue='duplicate_gene_id', gene_id=gene_id))
                # Keep evidence for both rows; downstream joins exclude ambiguous IDs.
            seen[gene_id] = True
            genes.append(gene)
    if not declared_gff3 or not genes:
        raise DataError('GFF3 has no usable gene features or version declaration.')
    return genes, headers, errors, build_seen
