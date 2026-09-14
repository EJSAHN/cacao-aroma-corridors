"""Project only unique, sequence-backed marker identities into a separate table.

A linked Nacional marker is not an independent sequence validation. No genotype
allele polarity, causal effect, LD proxy, or source association boundary is inferred.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from .common import DataError, text


def group(rows, field):
    result = defaultdict(list)
    for r in rows:
        v = text(r.get(field))
        if v:
            result[v].append(r)
    return result


def project_markers(inventory: list[dict], mapping: list[dict], direct_resource: str,
                    linked_resources: list[str]):
    byrow = {int(r['source_row']): r for r in mapping}
    bykey = group(mapping, 'reference_key')
    if len(byrow) != len(mapping):
        raise DataError('Duplicate mapping source rows.')
    marker_counts = Counter((r['resource'], text(r.get('canonical_id'))) for r in inventory)
    out = []
    for old in inventory:
        resource = text(old['resource'])
        key = text(old.get('canonical_id'))
        r = dict(resource=resource, original_row=old['source_row'], original_marker_id=old['marker_id'],
                 original_marker_key=key, original_chromosome=old.get('chromosome'),
                 original_position_bp=old.get('position_bp'), identity_basis='', coordinate_status='',
                 accepted_for_coordinate_comparison=False, genotype_allele_polarity='NOT_ESTABLISHED')
        candidates = []
        if resource == direct_resource:
            m = byrow.get(int(old['source_row']))
            if m is not None and m['source_marker_id'] != old['marker_id']:
                raise DataError('Evidence/source Diversity source-row identity mismatch.')
            candidates = [m] if m is not None else []
            r['identity_basis'] = 'own_flank_variation_reference_tag'
        elif resource in linked_resources:
            candidates = bykey.get(key, [])
            r['identity_basis'] = 'unique_reference_sequence_key_link_not_independent_tag'
        else:
            r['coordinate_status'] = 'NO_SEQUENCE_IDENTITY_BRIDGE'
            out.append(r)
            continue
        if not candidates:
            r['coordinate_status'] = 'NO_REFERENCE_TAG_KEY'
        elif len(candidates) != 1:
            r['coordinate_status'] = 'AMBIGUOUS_REFERENCE_TAG_KEY'
        else:
            m = candidates[0]
            r['mapped_query_id'] = m['query_id']
            r['source_tag_key'] = m['reference_key']
            r['tag_mapping_status'] = m['mapping_status']
            if m['mapping_status'] != 'UNIQUE_EXACT_PLACEMENT':
                r['coordinate_status'] = m['mapping_status']
            else:
                r.update(reference_sequence_id=m['mapped_sequence_id'], reference_chromosome=m.get('mapped_chromosome'),
                         reference_snp_bp=m['mapped_snp_bp'], reference_tag_start_bp=m['tag_start_bp'],
                         reference_tag_end_bp=m['tag_end_bp'], source_tag_strand_on_reference=m['mapped_strand'],
                         assembly_base=m['assembly_base'], other_source_base=m['other_source_base'])
                r['coordinate_status'] = 'DIRECT_UNIQUE_EXACT_TAG' if resource == direct_resource else 'PROJECTED_VIA_UNIQUE_TAG_KEY'
                if marker_counts[(resource, key)] > 1:
                    r['coordinate_status'] = 'DUPLICATE_PANEL_MARKER_KEY_REVIEW'
                r['accepted_for_coordinate_comparison'] = bool(m.get('primary_chromosome')) and marker_counts[(resource, key)] == 1
                r['source_chromosome_agrees'] = old.get('chromosome') == m.get('mapped_chromosome') if m.get('primary_chromosome') else None
                r['source_minus_reference_snp_bp'] = old['position_bp'] - m['mapped_snp_bp'] if r['source_chromosome_agrees'] and old.get('position_bp') is not None else None
                r['source_equals_tag_start_1based'] = old.get('position_bp') == m['tag_start_bp'] if r['source_chromosome_agrees'] else None
                r['source_equals_tag_start_0based'] = old.get('position_bp') == m['tag_start_bp'] - 1 if r['source_chromosome_agrees'] else None
        out.append(r)
    counts = Counter((r['resource'], r['coordinate_status']) for r in out)
    summary = [dict(resource=n, coordinate_status=s, records=count,
                    eligible_primary_coordinates=sum(x['accepted_for_coordinate_comparison'] for x in out if x['resource'] == n and x['coordinate_status'] == s))
               for (n, s), count in sorted(counts.items())]
    return out, summary


def candidate_checks(candidates: list[dict], genes: list[dict], aliases: dict):
    # Match exact source gene IDs, after removal of GFF's structural 'gene:' prefix.
    # Name is a separately labeled fallback; ambiguous matches are never selected.
    ids = group(genes, 'gene_id')
    names = group(genes, 'gene_name')
    out = []
    for old in candidates:
        gid = text(old['gene_id'])
        gs = ids.get(gid, [])
        basis = 'GFF_ID'
        if not gs:
            gs = names.get(gid, [])
            basis = 'GFF_Name_exact'
        r = dict(source_gene_id=gid, source_chromosome=old.get('chromosome'),
                 source_start_bp=old.get('start_bp'), source_end_bp=old.get('end_bp'),
                 source_annotation_records=old.get('record_count'), source_pathways=old.get('pathways', ''),
                 matching_annotation_records=len(gs), identity_field=basis,
                 accepted_original_interval=False)
        if not gs:
            r['status'] = 'SOURCE_GENE_ID_NOT_FOUND_IN_THIS_ANNOTATION'
        elif len(gs) != 1:
            r['status'] = 'AMBIGUOUS_ANNOTATION_ID'
        else:
            g = gs[0]
            r.update(annotation_gene_id=g['gene_id'], annotation_sequence_id=g['sequence_id'],
                     annotation_chromosome=g['chromosome'], annotation_start_bp=g['start_bp'],
                     annotation_end_bp=g['end_bp'], annotation_strand=g['strand'], biotype=g['biotype'])
            same_chr = g['chromosome'] is not None and old.get('chromosome') == g['chromosome']
            same = same_chr and old.get('start_bp') == g['start_bp'] and old.get('end_bp') == g['end_bp']
            r['status'] = 'EXACT_ID_AND_INTERVAL' if same else 'ID_PRESENT_COORDINATES_DIFFER'
            r['accepted_original_interval'] = same
        out.append(r)
    return out


def peak_links(peaks: list[dict], marker_rows: list[dict], mapped_rows: list[dict]):
    """Audit coordinate co-location; do not infer source causal SNP identity.

A same-coordinate match is only supporting coordinate context. It does NOT
establish that an association was generated by that SNP/allele.
"""
    index = defaultdict(list)
    m_by_row = {(r['resource'], r['original_row']): r for r in mapped_rows}
    for old in marker_rows:
        key = (old.get('chromosome'), old.get('position_bp'))
        index[key].append(m_by_row[(old['resource'], old['source_row'])])
    out = []
    for p in peaks:
        rs = index.get((p['chromosome'], p['start_bp']), [])
        r = dict(query_id=p['query_id'], source_chromosome=p['chromosome'], source_peak_bp=p['start_bp'],
                 source_events=p.get('n_events'), source_traits=p.get('n_traits'),
                 colocated_marker_rows=len(rs), unique_exact_tag_linked_rows=0,
                 coordinate_rewrite_applied=False, association_marker_identity='NOT_ESTABLISHED')
        validated = [x for x in rs if x.get('accepted_for_coordinate_comparison')]
        unchanged = [x for x in validated if x['reference_chromosome'] == p['chromosome'] and x['reference_snp_bp'] == p['start_bp']]
        r['unique_exact_tag_linked_rows'] = len(validated)
        r['same_original_and_reference_coordinate_rows'] = len(unchanged)
        r['coordinate_context'] = 'COLOCATED_REFERENCE_SNP' if unchanged else 'ONLY_DIFFERENT_REFERENCE_POSITION' if validated else 'NO_SEQUENCE_BACKED_COORDINATE_CONTEXT'
        r['interpretation'] = 'No association SNP identity, source boundary, effect direction or independent replication is inferred from coordinate co-location.'
        out.append(r)
    return out
