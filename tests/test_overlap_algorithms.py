from __future__ import annotations

import pandas as pd

from cacao_aroma_pipeline.analysis.overlaps import _nearest_interval_marker_details, _nearest_marker_details


def test_nearest_marker_details_basic():
    queries = pd.DataFrame(
        {
            "association_id": [1, 2, 3],
            "chromosome": [1, 1, 2],
            "peak_position_bp": [100, 210, 500],
        }
    )
    targets = pd.DataFrame(
        {
            "marker_id": ["m1", "m2", "m3"],
            "chromosome": [1, 1, 2],
            "position": [90, 300, 480],
        }
    )
    result = _nearest_marker_details(
        queries,
        targets,
        query_chr_col="chromosome",
        query_pos_col="peak_position_bp",
        query_id_col="association_id",
    )
    assert list(result["nearest_marker_id"]) == ["m1", "m2", "m3"]
    assert list(result["distance_bp"]) == [10, 90, 20]


def test_nearest_interval_marker_details_gene_body_and_flanking():
    genes = pd.DataFrame(
        {
            "candidate_gene_id": [1, 2],
            "chromosome": [1, 1],
            "start_bp": [100, 400],
            "end_bp": [200, 450],
        }
    )
    targets = pd.DataFrame(
        {
            "marker_id": ["m1", "m2"],
            "chromosome": [1, 1],
            "position": [150, 500],
        }
    )
    result = _nearest_interval_marker_details(genes, targets)
    assert result.loc[result["candidate_gene_id"] == 1, "inside_gene_body"].item() is True
    assert result.loc[result["candidate_gene_id"] == 1, "distance_bp"].item() == 0
    assert result.loc[result["candidate_gene_id"] == 2, "inside_gene_body"].item() is False
    assert result.loc[result["candidate_gene_id"] == 2, "distance_bp"].item() == 50
