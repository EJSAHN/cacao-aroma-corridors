from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from cacao_aroma_pipeline.models import ParsedDataset, SupplementaryPayload


def _prepare_coords(df: pd.DataFrame, chr_col: str, pos_col: str, keep_cols: list[str] | None = None) -> pd.DataFrame:
    if df.empty or chr_col not in df.columns or pos_col not in df.columns:
        return pd.DataFrame(columns=[chr_col, pos_col])
    keep = [chr_col, pos_col]
    if keep_cols:
        keep.extend([col for col in keep_cols if col in df.columns])
    frame = df[keep].copy()
    frame = frame.dropna(subset=[chr_col, pos_col]).copy()
    frame[chr_col] = frame[chr_col].astype(str)
    frame[pos_col] = pd.to_numeric(frame[pos_col], errors="coerce")
    frame = frame.dropna(subset=[pos_col]).copy()
    frame[pos_col] = frame[pos_col].astype(int)
    return frame.sort_values([chr_col, pos_col]).reset_index(drop=True)


def _nearest_marker_details(
    queries: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    query_chr_col: str,
    query_pos_col: str,
    target_chr_col: str = "chromosome",
    target_pos_col: str = "position",
    query_id_col: str,
    target_id_col: str = "marker_id",
) -> pd.DataFrame:
    queries_p = _prepare_coords(queries, query_chr_col, query_pos_col, keep_cols=[query_id_col])
    targets_p = _prepare_coords(targets, target_chr_col, target_pos_col, keep_cols=[target_id_col])

    results = []
    if queries_p.empty or targets_p.empty:
        return pd.DataFrame(
            columns=[
                query_id_col,
                "chromosome",
                "query_position",
                "nearest_marker_id",
                "nearest_marker_position",
                "distance_bp",
            ]
        )

    target_groups = {chrom: grp.reset_index(drop=True) for chrom, grp in targets_p.groupby(target_chr_col, sort=False)}
    for row in queries_p.itertuples(index=False):
        chrom = getattr(row, query_chr_col)
        position = int(getattr(row, query_pos_col))
        qid = getattr(row, query_id_col)
        group = target_groups.get(chrom)
        if group is None or group.empty:
            results.append(
                {
                    query_id_col: qid,
                    "chromosome": chrom,
                    "query_position": position,
                    "nearest_marker_id": None,
                    "nearest_marker_position": None,
                    "distance_bp": None,
                }
            )
            continue
        positions = group[target_pos_col].to_numpy()
        idx = np.searchsorted(positions, position)
        candidates = []
        if idx < len(group):
            candidates.append(group.iloc[idx])
        if idx > 0:
            candidates.append(group.iloc[idx - 1])
        best = None
        best_distance = None
        for candidate in candidates:
            distance = abs(int(candidate[target_pos_col]) - position)
            if best is None or distance < best_distance:
                best = candidate
                best_distance = distance
        results.append(
            {
                query_id_col: qid,
                "chromosome": chrom,
                "query_position": position,
                "nearest_marker_id": None if best is None else best[target_id_col],
                "nearest_marker_position": None if best is None else int(best[target_pos_col]),
                "distance_bp": best_distance,
            }
        )
    return pd.DataFrame(results)


def _nearest_interval_marker_details(
    genes: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    gene_id_col: str = "candidate_gene_id",
) -> pd.DataFrame:
    genes_p = genes.copy()
    genes_p = genes_p.dropna(subset=["chromosome", "start_bp", "end_bp"]).copy()
    genes_p["chromosome"] = genes_p["chromosome"].astype(str)
    genes_p["start_bp"] = pd.to_numeric(genes_p["start_bp"], errors="coerce")
    genes_p["end_bp"] = pd.to_numeric(genes_p["end_bp"], errors="coerce")
    genes_p = genes_p.dropna(subset=["start_bp", "end_bp"]).copy()
    genes_p["start_bp"] = genes_p["start_bp"].astype(int)
    genes_p["end_bp"] = genes_p["end_bp"].astype(int)

    targets_p = _prepare_coords(targets, "chromosome", "position", keep_cols=["marker_id"])
    target_groups = {chrom: grp.reset_index(drop=True) for chrom, grp in targets_p.groupby("chromosome", sort=False)}
    rows = []
    for gene in genes_p.itertuples(index=False):
        chrom = str(gene.chromosome)
        start = int(gene.start_bp)
        end = int(gene.end_bp)
        gid = getattr(gene, gene_id_col)
        group = target_groups.get(chrom)
        if group is None or group.empty:
            rows.append(
                {
                    gene_id_col: gid,
                    "chromosome": chrom,
                    "start_bp": start,
                    "end_bp": end,
                    "nearest_marker_id": None,
                    "nearest_marker_position": None,
                    "distance_bp": None,
                    "inside_gene_body": False,
                }
            )
            continue
        positions = group["position"].to_numpy()
        left_idx = np.searchsorted(positions, start)
        right_idx = np.searchsorted(positions, end, side="right")
        if left_idx < right_idx:
            inside = group.iloc[left_idx]
            rows.append(
                {
                    gene_id_col: gid,
                    "chromosome": chrom,
                    "start_bp": start,
                    "end_bp": end,
                    "nearest_marker_id": inside["marker_id"],
                    "nearest_marker_position": int(inside["position"]),
                    "distance_bp": 0,
                    "inside_gene_body": True,
                }
            )
            continue
        candidates = []
        if left_idx < len(group):
            candidates.append(group.iloc[left_idx])
        if left_idx > 0:
            candidates.append(group.iloc[left_idx - 1])
        best = None
        best_distance = None
        for candidate in candidates:
            pos = int(candidate["position"])
            distance = min(abs(pos - start), abs(pos - end))
            if best is None or distance < best_distance:
                best = candidate
                best_distance = distance
        rows.append(
            {
                gene_id_col: gid,
                "chromosome": chrom,
                "start_bp": start,
                "end_bp": end,
                "nearest_marker_id": None if best is None else best["marker_id"],
                "nearest_marker_position": None if best is None else int(best["position"]),
                "distance_bp": best_distance,
                "inside_gene_body": False,
            }
        )
    return pd.DataFrame(rows)


def build_pairwise_overlap_report(datasets: list[ParsedDataset], windows_bp: list[int]) -> dict[str, pd.DataFrame]:
    summary_rows = []
    detail_frames: list[pd.DataFrame] = []

    for left, right in itertools.combinations(sorted(datasets, key=lambda d: d.dataset_id), 2):
        left_markers = left.markers.copy()
        right_markers = right.markers.copy()

        left_ids = set(left_markers["marker_id"].dropna().astype(str))
        right_ids = set(right_markers["marker_id"].dropna().astype(str))
        shared_ids = left_ids & right_ids
        left_canonical = set(left_markers.get("canonical_marker_id", pd.Series(dtype=object)).dropna().astype(str))
        right_canonical = set(right_markers.get("canonical_marker_id", pd.Series(dtype=object)).dropna().astype(str))
        shared_canonical = left_canonical & right_canonical

        row = {
            "left_dataset_id": left.dataset_id,
            "right_dataset_id": right.dataset_id,
            "left_marker_count": left.marker_count,
            "right_marker_count": right.marker_count,
            "shared_marker_ids": len(shared_ids),
            "shared_marker_id_jaccard": (len(shared_ids) / len(left_ids | right_ids)) if (left_ids or right_ids) else None,
            "shared_canonical_marker_ids": len(shared_canonical),
            "shared_canonical_marker_id_jaccard": (len(shared_canonical) / len(left_canonical | right_canonical)) if (left_canonical or right_canonical) else None,
            "left_has_coordinates": left.has_coordinates,
            "right_has_coordinates": right.has_coordinates,
        }

        if left.has_coordinates and right.has_coordinates:
            exact = _nearest_marker_details(
                left_markers.rename(columns={"marker_id": "left_marker_id"}),
                right_markers,
                query_chr_col="chromosome",
                query_pos_col="position",
                query_id_col="left_marker_id",
            )
            detail = exact.copy()
            detail["left_dataset_id"] = left.dataset_id
            detail["right_dataset_id"] = right.dataset_id
            detail_frames.append(detail)
            for window in windows_bp:
                hit_count = int((exact["distance_bp"].fillna(window + 1) <= window).sum())
                row[f"left_markers_within_{window}_bp_of_right"] = hit_count
            exact_coords = int((exact["distance_bp"] == 0).sum())
            row["exact_coordinate_matches"] = exact_coords
        else:
            for window in windows_bp:
                row[f"left_markers_within_{window}_bp_of_right"] = None
            row["exact_coordinate_matches"] = None

        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    detail = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    return {
        "pairwise_overlap_summary": summary,
        "pairwise_nearest_marker_detail": detail,
    }


def build_association_overlap_report(
    supplementary: SupplementaryPayload,
    datasets: list[ParsedDataset],
    windows_bp: list[int],
) -> dict[str, pd.DataFrame]:
    associations = supplementary.associations.copy()
    if associations.empty:
        return {"association_overlap_summary": pd.DataFrame(), "association_nearest_marker_detail": pd.DataFrame()}

    associations["chromosome"] = associations["chromosome"].astype("Int64")
    summary_rows = []
    detail_frames: list[pd.DataFrame] = []

    for dataset in datasets:
        if not dataset.has_coordinates:
            continue
        nearest = _nearest_marker_details(
            associations,
            dataset.markers,
            query_chr_col="chromosome",
            query_pos_col="peak_position_bp",
            query_id_col="association_id",
        )
        merged = associations.merge(nearest, on="association_id", how="left")
        merged["dataset_id"] = dataset.dataset_id
        merged["dataset_name"] = dataset.dataset_name
        detail_frames.append(merged)
        total = len(merged)
        for window in windows_bp:
            hit_count = int((merged["distance_bp"].fillna(window + 1) <= window).sum())
            summary_rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "dataset_name": dataset.dataset_name,
                    "window_bp": window,
                    "total_associations": total,
                    "associations_with_hit": hit_count,
                    "fraction_with_hit": (hit_count / total) if total else None,
                    "unique_traits_with_hit": merged.loc[merged["distance_bp"].fillna(window + 1) <= window, "trait"].nunique(),
                }
            )
    summary = pd.DataFrame(summary_rows)
    detail = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    return {
        "association_overlap_summary": summary,
        "association_nearest_marker_detail": detail,
    }


def build_candidate_gene_overlap_report(
    supplementary: SupplementaryPayload,
    datasets: list[ParsedDataset],
    windows_bp: list[int],
) -> dict[str, pd.DataFrame]:
    genes = supplementary.candidate_genes.copy()
    if genes.empty:
        return {"candidate_overlap_summary": pd.DataFrame(), "candidate_nearest_marker_detail": pd.DataFrame(), "candidate_pathway_summary": pd.DataFrame()}

    detail_frames = []
    summary_rows = []
    pathway_rows = []

    pathway_cols = [
        col
        for col in [
            "pathway_monoterpene",
            "pathway_l_phe",
            "pathway_fatty_acid",
            "pathway_sugar",
            "pathway_maillard",
        ]
        if col in genes.columns
    ]

    for dataset in datasets:
        if not dataset.has_coordinates:
            continue
        nearest = _nearest_interval_marker_details(genes, dataset.markers)
        merged = genes.merge(nearest, on="candidate_gene_id", how="left", suffixes=("", "_nearest"))
        merged["dataset_id"] = dataset.dataset_id
        merged["dataset_name"] = dataset.dataset_name
        detail_frames.append(merged)

        total = len(merged)
        for window in windows_bp:
            hit_mask = merged["distance_bp"].fillna(window + 1) <= window
            summary_rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "dataset_name": dataset.dataset_name,
                    "window_bp": window,
                    "total_candidate_genes": total,
                    "candidate_genes_with_hit": int(hit_mask.sum()),
                    "fraction_with_hit": (float(hit_mask.mean()) if total else None),
                    "gene_body_hits": int((merged["distance_bp"] == 0).sum()),
                }
            )
            for pathway_col in pathway_cols:
                pathway_mask = merged[pathway_col].fillna(False).astype(bool)
                if not pathway_mask.any():
                    continue
                pathway_rows.append(
                    {
                        "dataset_id": dataset.dataset_id,
                        "dataset_name": dataset.dataset_name,
                        "window_bp": window,
                        "pathway": pathway_col,
                        "pathway_gene_count": int(pathway_mask.sum()),
                        "pathway_gene_hits": int((hit_mask & pathway_mask).sum()),
                        "fraction_pathway_genes_hit": float(((hit_mask & pathway_mask).sum()) / pathway_mask.sum()),
                    }
                )

    return {
        "candidate_overlap_summary": pd.DataFrame(summary_rows),
        "candidate_nearest_marker_detail": pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame(),
        "candidate_pathway_summary": pd.DataFrame(pathway_rows),
    }


def build_pipeline_summary_report(
    datasets: list[ParsedDataset],
    supplementary: SupplementaryPayload,
    inventory: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    dataset_rows = []
    for dataset in datasets:
        dataset_rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "dataset_name": dataset.dataset_name,
                "parser_type": dataset.parser_type,
                "source_path": dataset.source_path,
                "source_archive": dataset.source_archive,
                "marker_count": dataset.marker_count,
                "sample_count": dataset.sample_count,
                "has_genotypes": dataset.has_genotypes,
                "has_coordinates": dataset.has_coordinates,
                "note_count": len(dataset.notes),
                "notes": " | ".join(dataset.notes),
            }
        )

    supplementary_rows = [
        {
            "source_files": " | ".join(supplementary.source_files),
            "sensory_trait_count": len(supplementary.sensory_traits),
            "association_count": len(supplementary.associations),
            "association_summary_count": len(supplementary.association_summary),
            "candidate_gene_count": len(supplementary.candidate_genes),
            "notes": " | ".join(supplementary.notes),
        }
    ]
    return {
        "dataset_summary": pd.DataFrame(dataset_rows),
        "supplementary_summary": pd.DataFrame(supplementary_rows),
        "raw_inventory_snapshot": inventory.copy(),
    }
