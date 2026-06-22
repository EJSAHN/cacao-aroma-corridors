from __future__ import annotations

from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.models import ParsedDataset
from cacao_aroma_pipeline.utils import slugify, standardize_marker_frame, standardize_sample_frame


def parse_tassel_like_matrix(
    path: Path,
    *,
    dataset_name: str | None = None,
    source_archive: str | None = None,
) -> ParsedDataset:
    excel = pd.ExcelFile(path)
    sheet_name = excel.sheet_names[0]
    matrix = excel.parse(sheet_name, header=0)
    matrix.columns = [str(col).strip() for col in matrix.columns]

    dataset_name = dataset_name or path.stem
    dataset_id = slugify(dataset_name)
    parser_type = "tassel_like_matrix"

    if {"rs#", "chrom", "pos"}.issubset(matrix.columns):
        marker_id_col = "rs#"
        metadata_cols = [c for c in ["rs#", "chrom", "pos"] if c in matrix.columns]
        if "alleles" in matrix.columns:
            metadata_cols.append("alleles")
        if "strand" in matrix.columns:
            metadata_cols.append("strand")
    elif {"marker_id", "alleles", "chrom", "pos", "strand"}.issubset(matrix.columns):
        marker_id_col = "marker_id"
        metadata_cols = [c for c in ["marker_id", "alleles", "chrom", "pos", "strand"] if c in matrix.columns]
    else:
        raise ValueError(f"Unsupported TASSEL-like header structure: {path}")

    sample_cols = [col for col in matrix.columns if col not in metadata_cols]
    marker_frame = pd.DataFrame(
        {
            "marker_id": matrix[marker_id_col],
            "chromosome": matrix.get("chrom"),
            "position": matrix.get("pos"),
            "alleles": matrix.get("alleles"),
            "strand": matrix.get("strand"),
        }
    )
    marker_frame = standardize_marker_frame(
        marker_frame,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
    )

    samples = pd.DataFrame({"sample_id": sample_cols})
    samples["germplasm_name"] = samples["sample_id"]
    samples["dna_sample_id"] = samples["sample_id"]
    samples = standardize_sample_frame(
        samples,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
    )

    genotype_matrix = matrix[sample_cols].copy()
    genotype_matrix.index = matrix[marker_id_col].astype(str)

    metadata_rows = [
        ("dataset_name", dataset_name),
        ("sheet_name", sheet_name),
        ("marker_count", int(marker_frame["marker_id"].nunique(dropna=True))),
        ("sample_count", int(len(sample_cols))),
        ("source_path", str(path)),
        ("source_archive", source_archive or ""),
    ]
    metadata = pd.DataFrame(metadata_rows, columns=["field", "value"])

    notes = []
    if genotype_matrix.astype(str).eq("M").any().any():
        notes.append("Observed 'M' allele code in genotype matrix; downstream interpretation should be dataset-aware.")

    return ParsedDataset(
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        parser_type=parser_type,
        source_path=str(path),
        source_archive=source_archive,
        file_hash="",
        markers=marker_frame,
        samples=samples,
        genotype_matrix=genotype_matrix,
        metadata=metadata,
        raw_marker_table=matrix[metadata_cols].copy(),
        notes=notes,
    )
