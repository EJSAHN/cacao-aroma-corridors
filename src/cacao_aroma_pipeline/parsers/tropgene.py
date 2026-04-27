from __future__ import annotations

from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.models import ParsedDataset
from cacao_aroma_pipeline.utils import (
    dataframe_from_key_value_row,
    rename_with_fallback,
    slugify,
    standardize_marker_frame,
    standardize_sample_frame,
)


MARKER_NAME_CANDIDATES = {
    "marker name": "marker_id",
    "snp marker name": "marker_id",
    "reference sequence name": "reference_sequence_name",
    "chromosome": "chromosome",
    "snp position": "position",
    "variation": "alleles",
    "strand": "strand",
}


SAMPLE_CANDIDATES = {
    "dna sample id": "sample_id",
    "germplasm name": "germplasm_name",
    "dna sample origin": "origin",
    "accession number": "accession_number",
}


def _load_sheet(excel: pd.ExcelFile, sheet_name: str) -> pd.DataFrame:
    return excel.parse(sheet_name, header=0)


def _standardize_sample_table(
    sample_df: pd.DataFrame,
    *,
    dataset_id: str,
    dataset_name: str,
    parser_type: str,
    source_path: str,
    source_archive: str | None,
) -> pd.DataFrame:
    standardized = rename_with_fallback(sample_df.copy(), SAMPLE_CANDIDATES)
    if "sample_id" not in standardized.columns:
        if "dna_sample_id" in standardized.columns:
            standardized = standardized.rename(columns={"dna_sample_id": "sample_id"})
        else:
            first_col = standardized.columns[0]
            standardized = standardized.rename(columns={first_col: "sample_id"})
    if "germplasm_name" not in standardized.columns:
        standardized["germplasm_name"] = standardized["sample_id"]
    if "dna_sample_id" not in standardized.columns:
        standardized["dna_sample_id"] = standardized["sample_id"]
    columns = [col for col in ["sample_id", "germplasm_name", "dna_sample_id", "origin", "accession_number"] if col in standardized.columns]
    standardized = standardized[columns].drop_duplicates()
    return standardize_sample_frame(
        standardized,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=source_path,
        source_archive=source_archive,
    )


def _standardize_marker_table(
    marker_df: pd.DataFrame,
    *,
    dataset_id: str,
    dataset_name: str,
    parser_type: str,
    source_path: str,
    source_archive: str | None,
) -> pd.DataFrame:
    standardized = rename_with_fallback(marker_df.copy(), MARKER_NAME_CANDIDATES)
    if "marker_id" not in standardized.columns:
        first_col = standardized.columns[0]
        standardized = standardized.rename(columns={first_col: "marker_id"})
    for col in ["chromosome", "position", "alleles", "strand"]:
        if col not in standardized.columns:
            standardized[col] = pd.NA
    keep_cols = ["marker_id", "chromosome", "position", "alleles", "strand"] + [
        col for col in standardized.columns if col not in {"marker_id", "chromosome", "position", "alleles", "strand"}
    ]
    standardized = standardized[keep_cols].drop_duplicates()
    return standardize_marker_frame(
        standardized,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=source_path,
        source_archive=source_archive,
    )


def _normalize_matrix_orientation(matrix: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    top_left = str(matrix.columns[0]).strip().lower()
    if top_left == "marker/dnasample":
        markers_as_rows = matrix.copy()
        markers_as_rows = markers_as_rows.rename(columns={markers_as_rows.columns[0]: "marker_id"})
        markers_as_rows["marker_id"] = markers_as_rows["marker_id"].astype(str)
        markers_as_rows = markers_as_rows.set_index("marker_id")
        return markers_as_rows, "markers_as_rows"
    if top_left == "dnasample/marker":
        samples_as_rows = matrix.copy()
        samples_as_rows = samples_as_rows.rename(columns={samples_as_rows.columns[0]: "sample_id"})
        samples_as_rows["sample_id"] = samples_as_rows["sample_id"].astype(str)
        markers_as_rows = samples_as_rows.set_index("sample_id").T
        markers_as_rows.index.name = "marker_id"
        return markers_as_rows, "samples_as_rows"
    raise ValueError("Unsupported TropGene data_matrix orientation.")


def parse_tropgene_study_workbook(
    path: Path,
    *,
    dataset_name: str | None = None,
    source_archive: str | None = None,
) -> ParsedDataset:
    excel = pd.ExcelFile(path)
    dataset_name = dataset_name or path.stem
    dataset_id = slugify(dataset_name)
    parser_type = "tropgene_study_workbook"

    study_raw = excel.parse("study", header=None)
    metadata = dataframe_from_key_value_row(study_raw.iloc[0].tolist(), study_raw.iloc[1].tolist())

    sample_raw = _load_sheet(excel, "dna_sample")
    marker_raw = _load_sheet(excel, "marker")
    matrix_raw = _load_sheet(excel, "data_matrix")

    samples = _standardize_sample_table(
        sample_raw,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
    )
    markers = _standardize_marker_table(
        marker_raw,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
    )

    genotype_matrix, orientation = _normalize_matrix_orientation(matrix_raw)
    notes = [f"Original TropGene data_matrix orientation: {orientation}"]

    metadata = pd.concat(
        [
            metadata,
            pd.DataFrame(
                [
                    ("dataset_name", dataset_name),
                    ("source_path", str(path)),
                    ("source_archive", source_archive or ""),
                    ("matrix_orientation", orientation),
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
        samples=samples,
        genotype_matrix=genotype_matrix,
        metadata=metadata,
        raw_marker_table=marker_raw,
        raw_sample_table=sample_raw,
        notes=notes,
    )
