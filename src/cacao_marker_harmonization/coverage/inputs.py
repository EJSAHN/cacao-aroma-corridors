"""Build fixed query sets and identical-record coordinate cohorts."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path, PurePosixPath
from .common import DataError, contains, integer, text
from ..io.tableio import TableBook
from .numerics import Queries


def read_sheet(path: Path, sheet: str, required=()) -> list[dict]:
    with TableBook(path) as book:
        if sheet not in book.sheet_names:
            raise DataError(f'Missing sheet {sheet} in {path.name}')
        rows = [r for _, r in book.table(sheet)]
        if book.errors:
            raise DataError(f'Excel error cells in {path.name}:{sheet}')
    if rows and set(rows[0]) == {'No_records'}:
        rows = []
    for r in rows:
        if not set(required).issubset(r):
            raise DataError(f'Missing required columns in {path.name}:{sheet}: {list(required)}')
    return rows


def strict_bool(v, context: str) -> bool:
    if v is True or v == 1:
        return True
    if v is False or v == 0:
        return False
    raise DataError('Expected an explicit Boolean: ' + context)


def safe_filename(name) -> str:
    name = text(name)
    if not name or name in ('.', '..') or '/' in name or '\\' in name or ':' in name:
        raise DataError('Unsafe provenance filename: ' + name)
    return name








def source_path(root: Path, rel: str) -> Path:
    parts = PurePosixPath(rel.replace('\\', '/')).parts
    if not parts or any(p in ('.', '..') or ':' in p for p in parts) or rel.startswith(('/', '\\')):
        raise DataError('Unsafe relative source path.')
    p = root.joinpath(*parts)
    if not contains(root, p) or not p.is_file():
        raise DataError('Original input missing: ' + str(p))
    return p




def valid_point(chrom, pos, lengths: dict) -> bool:
    n = integer(pos)
    return text(chrom) in lengths and n is not None and 1 <= n <= lengths[text(chrom)]


def valid_interval(row: dict, lengths: dict) -> bool:
    s, e = integer(row.get('start_bp')), integer(row.get('end_bp'))
    return s is not None and e is not None and s <= e and valid_point(row.get('chromosome'), s, lengths) and valid_point(row.get('chromosome'), e, lengths)


def marker_key(resource, rownum):
    n = integer(rownum)
    if n is None or n < 1:
        raise DataError('Invalid original marker row.')
    return text(resource), n


def build_cohorts(inventory: list[dict], all_mapped: list[dict], lengths: dict, cfg: dict):
    """B and C always contain exactly the same eligible original record keys."""
    original, mapped = {}, {}
    for r in inventory:
        k = marker_key(r['resource'], r['source_row'])
        if k in original:
            raise DataError('Duplicate source marker record key.')
        original[k] = r
    for r in all_mapped:
        k = marker_key(r['resource'], r['original_row'])
        if k in mapped or k not in original:
            raise DataError('Duplicate or unknown mapped record key.')
        b = original[k]
        if text(b['marker_id']) != text(r['original_marker_id']) or text(b.get('chromosome')) != text(r.get('original_chromosome')) or integer(b.get('position_bp')) != integer(r.get('original_position_bp')):
            raise DataError('Original identity/coordinate drift between source and mapping.')
        mapped[k] = r
    if set(original) != set(mapped):
        raise DataError('Baseline and mapping inventories cover different marker records.')
    resources = cfg['comparison_resources'] + cfg['original_only_resources']
    cohorts, ledger, resource_states = {}, [], []
    for resource in resources:
        rr = [(k, original[k], mapped[k]) for k in sorted(original) if k[0] == resource]
        if not rr:
            raise DataError('Configured resource absent: ' + resource)
        A, B, C = [], [], []
        for k, old, m in rr:
            source_eligible = strict_bool(old['eligible_coordinate'], 'source eligible_coordinate')
            original_in_bounds = valid_point(old.get('chromosome'), old.get('position_bp'), lengths)
            eligible = source_eligible and original_in_bounds
            accepted = strict_bool(m['accepted_for_coordinate_comparison'], 'mapping acceptance')
            if accepted and not valid_point(m.get('reference_chromosome'), m.get('reference_snp_bp'), lengths):
                raise DataError('Accepted reference coordinate is outside the declared chromosome.')
            paired = eligible and accepted and resource in cfg['comparison_resources']
            rowid = str(k[1])
            base = dict(record_key=resource + ':row:' + rowid, original_row=k[1], marker_id=old['marker_id'])
            if eligible:
                A.append(dict(**base, chromosome=text(old['chromosome']), position_bp=integer(old['position_bp'])))
            if paired:
                B.append(dict(**base, chromosome=text(old['chromosome']), position_bp=integer(old['position_bp'])))
                C.append(dict(**base, chromosome=text(m['reference_chromosome']), position_bp=integer(m['reference_snp_bp'])))
            reason = 'PAIRED_B_C' if paired else 'CONTEXT_ONLY_ORIGINAL' if resource in cfg['original_only_resources'] else 'UNPAIRED_ORIGINAL_COORDINATE_INVALID' if accepted else m['coordinate_status']
            ledger.append(dict(resource=resource, **base, in_A=eligible, in_B=paired, in_C=paired,
                               mapping_accepted=accepted, cohort_reason=reason,
                               source_coordinate_eligible=source_eligible, original_within_primary_bounds=original_in_bounds,
                               original_chromosome=text(old.get('chromosome')), original_position_bp=integer(old.get('position_bp')),
                               reference_chromosome=text(m.get('reference_chromosome')), reference_snp_bp=integer(m.get('reference_snp_bp')),
                               identity_basis=m.get('identity_basis'), source_tag_key=m.get('source_tag_key')))
        cohorts[resource] = {'A': A}
        if resource in cfg['comparison_resources']:
            if not B:
                raise DataError('No paired eligible markers: ' + resource)
            assert [x['record_key'] for x in B] == [x['record_key'] for x in C]
            cohorts[resource].update(B=B, C=C)
        resource_states.append(dict(resource=resource, all_source_records=len(rr), A_records=len(A),
                                    paired_records=len(B),
                                    source_eligible_records=sum(strict_bool(o['eligible_coordinate'],resource) for _,o,_ in rr),
                                    source_eligible_outside_bounds=sum(strict_bool(o['eligible_coordinate'],resource) and not valid_point(o.get('chromosome'),o.get('position_bp'),lengths) for _,o,_ in rr),
                                    mapping_accepted_but_unpaired=sum(strict_bool(m['accepted_for_coordinate_comparison'],resource) and not (strict_bool(o['eligible_coordinate'],resource) and valid_point(o.get('chromosome'),o.get('position_bp'),lengths)) for _,o,m in rr),
                                    mapping_accepted_records=sum(strict_bool(m['accepted_for_coordinate_comparison'],resource) for _,_,m in rr),
                                    status='PAIRED_ANALYSIS' if B else 'ORIGINAL_COORDINATE_CONTEXT_ONLY'))
    return cohorts, ledger, resource_states


def build_queries(candidate_rows, annotation_checks, memberships, source_queries, lengths):
    cand = {r['gene_id']: r for r in candidate_rows}
    if len(cand) != len(candidate_rows):
        raise DataError('Candidate IDs are not unique.')
    check = {r['source_gene_id']: r for r in annotation_checks}
    if len(check) != len(annotation_checks) or set(check) != set(cand):
        raise DataError('Annotation candidate inventory differs from source.')
    accepted, query_audit = [], []
    for gid in sorted(cand):
        g, a = cand[gid], check[gid]
        if text(a['source_chromosome']) != text(g['chromosome']) or integer(a['source_start_bp']) != integer(g['start_bp']) or integer(a['source_end_bp']) != integer(g['end_bp']):
            raise DataError('Candidate source coordinates changed: ' + gid)
        ok = strict_bool(a['accepted_original_interval'], gid) and a['status'] == 'EXACT_ID_AND_INTERVAL' and valid_interval(g, lengths)
        if ok:
            if text(a['annotation_chromosome']) != text(g['chromosome']) or integer(a['annotation_start_bp']) != integer(g['start_bp']) or integer(a['annotation_end_bp']) != integer(g['end_bp']):
                raise DataError('Candidate acceptance contradicts annotation coordinates.')
            accepted.append(dict(query_id=gid, gene_id=gid, chromosome=text(g['chromosome']),
                                 start_bp=integer(g['start_bp']), end_bp=integer(g['end_bp']),
                                 biotype=a['biotype'], pathways=g['pathways'], source_annotation_rows=g['original_rows']))
        query_audit.append(dict(query_set='unique_candidate_genes', query_id=gid, included=ok, reason=a['status']))
    if not accepted:
        raise DataError('No annotation-verified candidate genes.')
    queries = {'unique_candidate_genes': Queries('unique_candidate_genes', accepted)}
    by_gene = {r['query_id']: r for r in accepted}
    pathway = defaultdict(set)
    for r in memberships:
        if r['gene_id'] not in cand:
            raise DataError('Pathway membership references an unknown candidate gene.')
        if r['gene_id'] in by_gene:
            pathway[r['pathway']].add(r['gene_id'])
    for name, ids in sorted(pathway.items()):
        label = 'pathway:' + name
        queries[label] = Queries(label, [by_gene[g] for g in sorted(ids)])
    for name, rr in source_queries.items():
        valid = []
        coords = set()
        for r in rr:
            ok = valid_interval(r, lengths)
            coord = (text(r['chromosome']), integer(r['start_bp']), integer(r['end_bp']))
            if coord in coords:
                raise DataError('Duplicate positions/intervals in source query set ' + name)
            coords.add(coord)
            if ok:
                valid.append(dict(query_id=r['query_id'], chromosome=coord[0], start_bp=coord[1], end_bp=coord[2]))
            query_audit.append(dict(query_set=name, query_id=r['query_id'], included=ok,
                                    reason='SOURCE_POSITION_RETAINED_NO_IDENTITY_INFERENCE' if ok else 'INVALID_OR_OUTSIDE_PRIMARY_CHROMOSOME'))
        if valid:
            queries[name] = Queries(name, valid)
    return queries, query_audit, set(cand)
