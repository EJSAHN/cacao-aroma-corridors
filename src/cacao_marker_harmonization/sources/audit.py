"""Summarize source provenance, sample labels, and alternative marker-identity matches."""
from __future__ import annotations
from itertools import combinations
from collections import Counter
import numpy as np
import pandas as pd
from .common import integer

def resource_summary(resources):
    return pd.DataFrame([dict(resource=r.spec['id'],
        label=r.spec['label'], role=r.spec['role'], source_file=r.path.name,
        reader=r.engine, marker_records=len(r.markers), unique_marker_ids=r.markers.marker_id.nunique(),
        unique_canonical_ids=r.markers.canonical_id.nunique(), sample_ids=len(r.samples),
        coordinate_rows=int(r.markers.eligible_coordinate.sum()), unique_placed_coordinates=len(r.markers[r.markers.eligible_coordinate].drop_duplicates(['chromosome',
        'position_bp'])), unplaced_rows=int((r.markers.has_position & ~r.markers.placed_chromosome).sum()),
        missing_or_invalid_positions=int((~r.markers.has_position).sum()),
        marker_reference_name_disagreements=int(((r.markers.reference_key != '') & (r.markers.reference_key != r.markers.canonical_id)).sum()),
        coordinate_review=r.spec['coordinate_status'], genotype_encoding=r.spec['genotype_encoding']) for r in resources])

def sample_overlap(resources):
    out = []
    for a, b in combinations(resources, 2):
        aa, bb = (set(a.samples), set(b.samples))
        out.append(dict(resource_a=a.spec['id'],
            resource_b=b.spec['id'], n_a=len(aa), n_b=len(bb), shared_ids=len(aa & bb),
            a_contained_in_b=bool(aa) and aa.issubset(bb), b_contained_in_a=bool(bb) and bb.issubset(aa),
            only_a='; '.join(sorted(aa - bb)), only_b='; '.join(sorted(bb - aa)),
            note='Shared labels are not a genotype identity test.'))
    return pd.DataFrame(out)

def identity_audit(resources):
    """Compare unambiguous identity keys without accepting alternative joins as remapping."""
    summaries = []
    details = []
    collisions = []
    id_summary = []
    for r in resources:
        for field in ['canonical_id', 'reference_key']:
            for val, g in r.markers[r.markers[field] != ''].groupby(field, sort=False):
                if len(g) > 1:
                    collisions.append(dict(resource=r.spec['id'],
                        key_type=field, key=val, n_rows=len(g), n_marker_ids=g.marker_id.nunique(),
                        n_coordinates=len(g[['chromosome', 'position_bp']].drop_duplicates())))
    for a, b in combinations(resources, 2):
        aa = set(a.markers.marker_id)
        bb = set(b.markers.marker_id)
        ca = set(a.markers.canonical_id)
        cb = set(b.markers.canonical_id)
        id_summary.append(dict(resource_a=a.spec['id'],
            resource_b=b.spec['id'], shared_raw_ids=len(aa & bb), raw_jaccard=len(aa & bb) / len(aa | bb) if aa | bb else None,
            shared_canonical_ids=len(ca & cb), canonical_jaccard=len(ca & cb) / len(ca | cb) if ca | cb else None))
        for af, bf in [('canonical_id',
            'canonical_id'), ('canonical_id', 'reference_key'), ('reference_key',
            'canonical_id')]:
            a0 = a.markers[a.markers[af] != '']
            b0 = b.markers[b.markers[bf] != '']
            if a0.empty or b0.empty:
                continue
            a1 = a0[~a0[af].duplicated(False)]
            b1 = b0[~b0[bf].duplicated(False)]
            am = a1.set_index(af, drop=False)
            bm = b1.set_index(bf, drop=False)
            shared = sorted(set(am.index) & set(bm.index))
            tally = Counter()
            for ident in shared:
                ra, rb = (am.loc[ident], bm.loc[ident])
                valid = bool(ra.eligible_coordinate and rb.eligible_coordinate)
                samechr = valid and ra.chromosome == rb.chromosome
                delta = integer(rb.position_bp - ra.position_bp) if samechr else None
                kind = 'unavailable_coordinate' if not valid else 'different_chromosome' if not samechr else 'exact_coordinate' if delta == 0 else 'same_chromosome_offset'
                tally[kind] += 1
                details.append(dict(resource_a=a.spec['id'],
                    resource_b=b.spec['id'], key_type_a=af, key_type_b=bf, matching_key=ident,
                    marker_a=ra.marker_id, marker_b=rb.marker_id, chromosome_a=ra.chromosome,
                    position_a_bp=integer(ra.position_bp), chromosome_b=rb.chromosome,
                    position_b_bp=integer(rb.position_bp), position_b_minus_a_bp=delta,
                    classification=kind, alleles_a=ra.alleles_reported, alleles_b=rb.alleles_reported,
                    note='Alternative reference-key joins are diagnostic hypotheses, not accepted remapping.' if af != bf else 'ID equivalence does not establish sequence or allele equivalence.'))
            summaries.append(dict(resource_a=a.spec['id'],
                resource_b=b.spec['id'], key_type_a=af, key_type_b=bf, shared_unique_keys=len(shared),
                ambiguous_a=int(a0[af].duplicated(False).sum()), ambiguous_b=int(b0[bf].duplicated(False).sum()),
                exact_coordinates=tally['exact_coordinate'], same_chromosome_offsets=tally['same_chromosome_offset'],
                different_chromosomes=tally['different_chromosome'], unavailable_coordinates=tally['unavailable_coordinate']))
    return (pd.DataFrame(id_summary), pd.DataFrame(summaries), pd.DataFrame(details), pd.DataFrame(collisions))

def genotype_symbols(resources):
    """Count literal symbols and per-marker symbol patterns without imputing genotypes."""
    totals = []
    patterns = []
    for r in resources:
        if r.genotypes is None:
            continue
        vals, ct = np.unique(r.genotypes, return_counts=True)
        totals.extend((dict(resource=r.spec['id'], symbol=str(v), count=int(n)) for v, n in zip(vals, ct)))
        count = Counter((tuple(sorted(set(row))) for row in r.genotypes))
        patterns.extend((dict(resource=r.spec['id'],
            symbols=';'.join(v), marker_rows=n) for v, n in count.most_common()))
    return (pd.DataFrame(totals), pd.DataFrame(patterns))

def gene_interval_links(genes, intervals):
    """Report candidate/source-interval positional relationships, not candidate discovery."""
    rows = []
    for r in genes.itertuples():
        s = intervals[intervals.chromosome.eq(r.chromosome)]
        overlapping = s[(s.start_bp <= r.end_bp) & (s.end_bp >= r.start_bp)]
        contained = overlapping[(overlapping.start_bp <= r.start_bp) & (overlapping.end_bp >= r.end_bp)]
        if s.empty:
            nearest = None
        else:
            distances = np.maximum(np.maximum(s.start_bp.to_numpy() - r.end_bp,
                r.start_bp - s.end_bp.to_numpy()), 0)
            nearest = int(distances.min())
        rows.append(dict(gene_id=r.gene_id,
            chromosome=r.chromosome, start_bp=r.start_bp, end_bp=r.end_bp,
            n_overlapping_source_intervals=len(overlapping), n_containing_source_intervals=len(contained),
            distance_to_nearest_source_interval_bp=nearest, overlapping_interval_ids=';'.join(overlapping.query_id),
            interpretation="Positional overlap only; not a reconstruction of the authors' candidate-selection rule."))
    return pd.DataFrame(rows)
