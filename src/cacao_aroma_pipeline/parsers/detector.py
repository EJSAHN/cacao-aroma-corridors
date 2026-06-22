from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.exceptions import DetectionError


def _safe_excel(path: Path) -> pd.ExcelFile | None:
    try:
        return pd.ExcelFile(path)
    except Exception:
        return None


def detect_parser_type(path: Path) -> str:
    stem = path.stem.lower()
    match = re.search(r"mmc([1-7])", stem)
    if match:
        return f"supplementary_{match.group(1)}"

    excel = _safe_excel(path)
    if excel is None:
        raise DetectionError(f"Unable to open spreadsheet for detection: {path}")

    sheets_lower = {sheet.lower(): sheet for sheet in excel.sheet_names}
    if {"study", "dna_sample", "marker", "data_matrix"}.issubset(sheets_lower):
        return "tropgene_study_workbook"
    if {"study", "marker"}.issubset(sheets_lower) and "data_matrix" not in sheets_lower:
        return "marker_metadata_workbook"

    first_sheet = excel.sheet_names[0]
    probe = excel.parse(first_sheet, nrows=3, header=None)
    row0 = [str(x).strip().lower() for x in probe.iloc[0].tolist()]
    if len(row0) >= 4 and row0[:3] == ["rs#", "chrom", "pos"]:
        return "tassel_like_matrix"
    if len(row0) >= 5 and row0[:5] == ["marker_id", "alleles", "chrom", "pos", "strand"]:
        return "tassel_like_matrix"

    probe_a1 = str(probe.iloc[0, 0]).strip().lower() if probe.shape[0] and probe.shape[1] else ""
    if probe_a1.startswith("table s1"):
        return "supplementary_1"
    if probe_a1.startswith("table s2"):
        return "supplementary_2"
    if probe_a1.startswith("table s3"):
        return "supplementary_3"
    if probe_a1.startswith("table s4"):
        return "supplementary_4"

    raise DetectionError(f"Unsupported spreadsheet structure: {path}")
