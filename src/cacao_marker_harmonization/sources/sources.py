"""Parse original resources and preserve source-level identity and annotation evidence."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .common import DataError, HeaderNotFound, canonical_marker, chromosome, integer, key, number, sha256, stable_id, text
from ..io.tableio import TableBook, field

@dataclass
class Resource:
    spec: dict
    path: Path
    markers: pd.DataFrame
    samples: list[str]
    genotypes: np.ndarray | None
    engine: str
    metadata: list[dict]

def locate(root: Path, relative: str, aliases: list[str], optional=False) -> Path | None:
    expected = root / relative
    if expected.is_file():
        return expected
    wanted = {x.casefold() for x in aliases}
    found = sorted((p for p in root.rglob('*') if p.is_file() and p.name.casefold() in wanted and ('original_archives' not in p.parts)))
    if not found:
        if optional:
            return None
        raise DataError(f'Input not found: {relative}. Searched only under {root}. Set the path in config/study.json; do not rename unrelated files.')
    if len(found) > 1 and len({sha256(p) for p in found}) != 1:
        raise DataError(f'Different copies of {relative} found: ' + '; '.join((str(p.relative_to(root)) for p in found)))
    return found[0]

def load_associations(path: Path, cfg: dict) -> dict[str, pd.DataFrame]:
    """Separate original worksheet records, master events, unique peaks, and inherited intervals."""
    records = []
    with TableBook(path) as book:
        if cfg['association_master_sheet'] not in book.sheet_names:
            raise DataError('Configured association master worksheet is missing.')
        for sheet in book.sheet_names:
            try:
                table = list(book.table(sheet, required=['Chromosome', 'Traits']))
            except HeaderNotFound:
                continue
            for rownum, row in table:
                chrom = chromosome(field(row, 'Chromosome'), cfg['chromosome_aliases'])
                peak = integer(field(row,
                    'Position of the association peak (bp)', 'Position of the association peak'))
                trait = text(field(row, 'Traits'))
                if not chrom and peak is None and (not trait):
                    continue
                start = integer(field(row, 'Position of haplotypic bloc start'))
                end = integer(field(row, 'Position of haplotypic bloc end'))
                block = integer(field(row, 'N° haplotypic bloc', 'N° hap. Bloc'))
                method = text(field(row, 'GWAS method'))
                sorting = text(field(row, 'Sorting of marker'))
                event_key = stable_id('a_', (chrom, peak, trait, method, sorting))
                valid = chrom in cfg['chromosomes'] and peak is not None and (peak > 0) and bool(trait)
                interval_ok = valid and start is not None and (end is not None) and (0 < start <= peak <= end)
                records.append(dict(source_file=path.name,
                    source_sheet=sheet, source_row=rownum, chromosome=chrom, peak_bp=peak,
                    trait=trait, method_reported=method, marker_filter_reported=sorting,
                    source_start_bp=start, source_end_bp=end, source_block_number=block,
                    p_value_reported=number(field(row, 'p-value of the strongest association')),
                    explained_variance_reported=number(field(row, 'Explanation rate of the trait of the strongest association')),
                    event_key=event_key, valid_peak=valid, interval_valid_as_labeled=bool(interval_ok),
                    method_field_warning=method in {'G7', 'MAF5'}, raw_record=json.dumps(row,
                    ensure_ascii=False, default=str)))
    raw = pd.DataFrame(records)
    if raw.empty:
        raise DataError('No association rows parsed.')
    master_name = cfg['association_master_sheet']
    raw['in_master'] = raw.source_sheet.eq(master_name)
    master = raw[raw.in_master & raw.valid_peak].copy()
    if master.empty:
        raise DataError('The association master has no valid peak rows.')
    master_intervals = {}
    for k, g in master.groupby('event_key', sort=False):
        ints = set(((int(r.source_start_bp),
            int(r.source_end_bp), integer(r.source_block_number)) for r in g.itertuples() if r.interval_valid_as_labeled))
        master_intervals[k] = ints
    audit = []
    for r in raw.itertuples():
        matched = master_intervals.get(r.event_key, set())
        actual = (integer(r.source_start_bp), integer(r.source_end_bp), integer(r.source_block_number))
        rotated = (integer(r.source_end_bp), integer(r.source_block_number), integer(r.source_start_bp))
        state = 'master' if r.in_master else 'not_in_master' if r.event_key not in master_intervals else 'matches_master_interval' if actual in matched else 'rotation_matches_master' if rotated in matched else 'disagrees_or_ambiguous'
        audit.append(dict(source_sheet=r.source_sheet,
            source_row=r.source_row, event_key=r.event_key, master_relation=state,
            source_start_bp=actual[0], source_end_bp=actual[1], source_block_number=actual[2],
            suggested_start_bp=rotated[0] if state == 'rotation_matches_master' else None,
            suggested_end_bp=rotated[1] if state == 'rotation_matches_master' else None))
    meta = []
    for k, g in master.groupby('event_key', sort=False):
        payloads = g[['source_start_bp',
            'source_end_bp', 'source_block_number', 'p_value_reported',
            'explained_variance_reported']].fillna('').astype(str).drop_duplicates()
        first = g.iloc[0].to_dict()
        first['master_occurrences'] = len(g)
        first['distinct_payloads'] = len(payloads)
        first['master_conflict'] = len(payloads) > 1
        first['all_source_rows'] = ';'.join(map(str, g.source_row))
        meta.append(first)
    events = pd.DataFrame(meta)
    fruity_keys = set(raw.loc[raw.source_sheet.eq(cfg['association_fruity_subset_sheet']) & raw.valid_peak,
        'event_key'])
    events['in_fruity_subset'] = events.event_key.isin(fruity_keys)

    def peaks(frame):
        out = []
        for (c, p), g in frame.groupby(['chromosome', 'peak_bp'], sort=True):
            out.append(dict(query_id=stable_id('peak_',
                (c, int(p))), chromosome=c, start_bp=int(p), end_bp=int(p),
                n_events=len(g), n_traits=g.trait.nunique(), traits='; '.join(sorted(set(g.trait)))))
        return pd.DataFrame(out,
            columns=['query_id', 'chromosome', 'start_bp', 'end_bp', 'n_events',
            'n_traits', 'traits'])
    ig = master[master.interval_valid_as_labeled].drop_duplicates(['chromosome',
        'source_start_bp', 'source_end_bp'])
    intervals = pd.DataFrame([dict(query_id=stable_id('interval_',
        (r.chromosome, r.source_start_bp, r.source_end_bp)), chromosome=r.chromosome,
        start_bp=int(r.source_start_bp), end_bp=int(r.source_end_bp)) for r in ig.itertuples()],
        columns=['query_id', 'chromosome', 'start_bp', 'end_bp'])
    extra = raw[raw.valid_peak & ~raw.event_key.isin(master.event_key)].copy()
    union_keys = raw[raw.valid_peak].drop_duplicates('event_key')
    counts = []
    for sheet, g in raw.groupby('source_sheet', sort=False):
        counts.append(dict(source_sheet=sheet,
            n_rows=len(g), valid_peaks=int(g.valid_peak.sum()), unique_event_keys=g.loc[g.valid_peak,
            'event_key'].nunique(), keys_absent_from_master=g.loc[g.valid_peak & ~g.event_key.isin(master.event_key),
            'event_key'].nunique(), invalid_intervals_as_labeled=int((~g.interval_valid_as_labeled).sum()),
            method_field_warnings=int(g.method_field_warning.sum())))
    return dict(raw_records=raw,
        worksheet_audit=pd.DataFrame(counts), interval_audit=pd.DataFrame(audit),
        master_events=events, extra_records=extra, unique_peaks=peaks(events),
        fruity_peaks=peaks(events[events.in_fruity_subset]), source_intervals=intervals,
        union_peaks_sensitivity=peaks(union_keys), master_raw=master)

def load_candidates(path: Path, cfg: dict) -> dict[str, pd.DataFrame]:
    """Consolidate consistent gene IDs and union only explicit positive pathway flags."""
    records = []
    qualifiers = []
    with TableBook(path) as book:
        for rownum, row in book.table(cfg['candidate_sheet'], required=['gene_id', 'start', 'end']):
            gene = text(field(row, 'gene_id'))
            if not gene:
                continue
            c = chromosome(field(row, 'Chromosome'), cfg['chromosome_aliases'])
            start = integer(field(row, 'start'))
            end = integer(field(row, 'end'))
            r = dict(gene_id=gene,
                chromosome=c, start_bp=start, end_bp=end, gene_function=text(field(row,
                'gene_function')), source_file=path.name, source_sheet=cfg['candidate_sheet'],
                source_row=rownum, valid_interval=c in cfg['chromosomes'] and start is not None and (end is not None) and (0 < start <= end))
            for label, col in cfg['pathway_columns'].items():
                val = text(field(row, col)).lower()
                r[label + '_source_flag'] = val
                if val not in {'', 'x', '1', '1.0', 'yes', 'true', '0', 'no', 'false'}:
                    qualifiers.append(dict(gene_id=gene,
                        pathway=label, source_row=rownum, source_flag=val, primary_membership=False,
                        interpretation='Qualitative annotation retained; not counted as an explicit positive pathway flag.'))
                r[label] = val in {'x', '1', '1.0', 'yes', 'true'}
            records.append(r)
    raw = pd.DataFrame(records)
    unique = []
    members = []
    conflicts = []
    for gene, g in raw.groupby('gene_id', sort=True):
        coords = g[['chromosome', 'start_bp', 'end_bp']].drop_duplicates()
        if len(coords) != 1 or not g.valid_interval.all():
            conflicts.extend(g.to_dict('records'))
            continue
        r = g.iloc[0].to_dict()
        r['query_id'] = gene
        r['record_count'] = len(g)
        r['original_rows'] = ';'.join(map(str, g.source_row))
        r['gene_function'] = '; '.join(sorted(set(g.gene_function) - {''}))
        flags = []
        for label in cfg['pathway_columns']:
            r[label] = bool(g[label].any())
            if r[label]:
                flags.append(label)
                members.append(dict(gene_id=gene, pathway=label))
        r['pathways'] = '; '.join(flags)
        r['has_pathway_annotation'] = bool(flags)
        unique.append(r)
    return dict(raw_records=raw,
        unique_genes=pd.DataFrame(unique), memberships=pd.DataFrame(members,
        columns=['gene_id', 'pathway']), coordinate_conflicts=pd.DataFrame(conflicts),
        qualifier_annotations=pd.DataFrame(qualifiers), duplicates=raw[raw.gene_id.duplicated(False)].copy())

def load_resource(path: Path, spec: dict, cfg: dict) -> Resource:
    """Read one resource and retain raw identifiers, coordinates, and encoding metadata."""
    markers = []
    samples = []
    matrix = None
    metadata = []
    with TableBook(path) as book:
        engine = book.engine
        if spec['format'] == 'matrix':
            rows = iter(book.rows(book.sheet_names[0]))
            _, header = next(rows)
            header = [text(x) for x in header]
            hn = {key(h): i for i, h in enumerate(header)}
            mid = hn.get('rs', hn.get('markerid'))
            if mid is None or 'chrom' not in hn or 'pos' not in hn:
                raise DataError(f'Unexpected genotype header: {path.name}')
            meta_keys = {'rs',
                'markerid', 'alleles', 'chrom', 'pos', 'strand', 'assembly',
                'center', 'protlsid', 'assaylsid', 'panellsid', 'qccode'}
            sample_idx = [i for i, h in enumerate(header) if key(h) not in meta_keys]
            samples = [header[i] for i in sample_idx]
            if len(set(samples)) != len(samples) or any((not s for s in samples)):
                raise DataError(f'Duplicate/empty sample labels: {path.name}')
            genotypes = []
            for rownum, values in rows:
                values = values + [None] * (len(header) - len(values))
                marker = text(values[mid])
                if not marker:
                    continue
                markers.append(dict(marker_id=marker,
                    chromosome_raw=text(values[hn['chrom']]), position_raw=text(values[hn['pos']]),
                    alleles_reported=text(values[hn['alleles']]) if 'alleles' in hn else '',
                    strand_reported=text(values[hn['strand']]) if 'strand' in hn else '',
                    reference_sequence_name='', flank_5='', flank_3='', assembly_reported='',
                    remark='', source_row=rownum))
                calls = [text(values[i]).upper() for i in sample_idx]
                if any((len(v) > 16 for v in calls)):
                    raise DataError(f'Non-genotype content in sample columns: {path.name}:{rownum}')
                genotypes.append(calls)
            matrix = np.array(genotypes)
        else:
            for rownum, row in book.table('marker',
                required=['marker name'] if spec['format'] == 'tropgene' else ['SNP marker name']):
                marker = text(field(row, 'SNP marker name', 'marker name'))
                if not marker:
                    continue
                markers.append(dict(marker_id=marker,
                    chromosome_raw=text(field(row, 'chromosome')), position_raw=text(field(row,
                    'snp position')), alleles_reported=text(field(row, 'variation')),
                    strand_reported=text(field(row, 'strand')), reference_sequence_name=text(field(row,
                    'reference sequence name')), flank_5=text(field(row, '5flank_sequence')),
                    flank_3=text(field(row, '3flank_sequence')), assembly_reported=text(field(row,
                    'version')), remark=text(field(row, 'remark marker')), source_row=rownum))
            if 'dna_sample' in book.sheet_names:
                samples = [text(field(r,
                    'DNA sample ID')) for _, r in book.table('dna_sample', required=['DNA sample ID']) if text(field(r,
                    'DNA sample ID'))]
                sample_counts = Counter(samples)
                metadata.extend([dict(resource=spec['id'],
                    field='sample_metadata_rows', value=str(len(samples))), dict(resource=spec['id'],
                    field='sample_metadata_unique_ids', value=str(len(sample_counts)))])
                metadata.extend((dict(resource=spec['id'],
                    field='repeated_sample_metadata_id', value=json.dumps({'sample_id': sid,
                    'occurrences': n}, ensure_ascii=False)) for sid, n in sample_counts.items() if n > 1))
                samples = list(sample_counts)
            if 'study' in book.sheet_names:
                rr = list(book.rows('study'))[:2]
                if len(rr) == 2:
                    metadata.extend((dict(resource=spec['id'],
                        field=text(h), value=text(v)) for h, v in zip(rr[0][1], rr[1][1]) if text(h)))
    metadata.extend((dict(resource=spec['id'],
        field='source_excel_error', value=json.dumps(e, ensure_ascii=False)) for e in book.errors))
    df = pd.DataFrame(markers)
    if df.empty:
        raise DataError(f'No markers found: {path.name}')
    df.insert(0, 'resource', spec['id'])
    df['chromosome'] = df.chromosome_raw.map(lambda v: chromosome(v, cfg['chromosome_aliases']))
    df['position_bp'] = df.position_raw.map(integer)
    df['canonical_id'] = df.marker_id.map(canonical_marker)
    df['reference_key'] = df.reference_sequence_name.map(canonical_marker)
    df['has_position'] = df.position_bp.map(lambda p: p is not None and (not pd.isna(p)) and (p > 0))
    df['placed_chromosome'] = df.chromosome.isin(cfg['chromosomes'])
    df['eligible_coordinate'] = df.has_position & df.placed_chromosome & ~df.marker_id.str.startswith('#')
    ambiguous = set(df.groupby('marker_id').filter(lambda g: len(g[['chromosome',
        'position_bp']].drop_duplicates()) > 1).marker_id)
    df['conflicting_marker_id'] = df.marker_id.isin(ambiguous)
    df['eligible_coordinate'] &= ~df.conflicting_marker_id
    df['record_index'] = np.arange(len(df))
    return Resource(spec, path, df, samples, matrix, engine, metadata)
