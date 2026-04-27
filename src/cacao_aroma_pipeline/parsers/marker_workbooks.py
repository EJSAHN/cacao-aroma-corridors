from __future__ import annotations

from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.models import ParsedDataset
from cacao_aroma_pipeline.utils import dataframe_from_key_value_row, rename_with_fallback, slugify, standardize_marker_frame


MARKER_MAPPING = {
    "snp marker name": "marker_id",
    "marker name": "marker_id",
    "reference sequence name": "reference_sequence_name",
    "chromosome": "chromosome",
    "snp position": "position",
    "variation": "alleles",
    "strand": "strand",
}


def parse_marker_metadata_workbook(
    path: Path,
    *,
    dataset_name: str | None = None,
    source_archive: str | None = None,
) -> ParsedDataset:
    excel = pd.ExcelFile(path)
    dataset_name = dataset_name or path.stem
    dataset_id = slugify(dataset_name)
    parser_type = "marker_metadata_workbook"

    metadata_rows = []
    if "study" in excel.sheet_names:
        study_raw = excel.parse("study", header=None)
        if len(study_raw) >= 2:
            metadata_rows.append(dataframe_from_key_value_row(study_raw.iloc[0].tolist(), study_raw.iloc[1].tolist()))
    marker_raw = excel.parse("marker", header=0)
    marker_std = rename_with_fallback(marker_raw.copy(), MARKER_MAPPING)
    if "marker_id" not in marker_std.columns:
        first_col = marker_std.columns[0]
        marker_std = marker_std.rename(columns={first_col: "marker_id"})
    for col in ["chromosome", "position", "alleles", "strand"]:
        if col not in marker_std.columns:
            marker_std[col] = pd.NA
    markers = standardize_marker_frame(
        marker_std,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
    )
    metadata = pd.concat(metadata_rows, ignore_index=True) if metadata_rows else pd.DataFrame(columns=["field", "value"])
    metadata = pd.concat(
        [
            metadata,
            pd.DataFrame(
                [
                    ("dataset_name", dataset_name),
                    ("source_path", str(path)),
                    ("source_archive", source_archive or ""),
                    ("marker_count", int(markers["marker_id"].nunique(dropna=True))),
                ],
                columns=["field", "value"],
            ),
        ],
        ignore_index=True,
    )
    return ParsedDataset(
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
        file_hash="",
        markers=markers,
        samples=pd.DataFrame(columns=["sample_id"]),
        genotype_matrix=pd.DataFrame(),
        metadata=metadata,
        raw_marker_table=marker_raw,
        notes=["Marker metadata workbook does not include a genotype matrix."],
    )
