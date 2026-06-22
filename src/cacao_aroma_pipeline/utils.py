from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from cacao_aroma_pipeline.constants import MAX_EXCEL_SHEETNAME_LEN


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def file_md5(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify(value: str) -> str:
    cleaned = re.sub(r"[^0-9a-zA-Z]+", "_", value.strip())
    cleaned = cleaned.strip("_").lower()
    return cleaned or "dataset"


def safe_sheet_name(name: str) -> str:
    cleaned = re.sub(r"[:\\/?*\[\]]+", "_", name).strip()
    return cleaned[:MAX_EXCEL_SHEETNAME_LEN] or "Sheet1"


def normalize_chromosome(value: object) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip()
    if not text:
        return None

    # Prefer explicit chromosome labels when present, e.g. "chromosome_4", "chr4".
    label_match = re.search(r"(?:chromosome|chr)[_\s-]*(\d+)\b", text, flags=re.IGNORECASE)
    if label_match:
        return str(int(label_match.group(1)))

    if re.fullmatch(r"\d+(?:\.0)?", text):
        return str(int(float(text)))

    return text




def canonicalize_marker_id(value: object) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip().lower()
    if not text or text == "nan":
        return None
    text = re.sub(r"[^0-9a-z]+", "", text)
    return text or None

def coerce_integer(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, float):
        if np.isnan(value):
            return None
        return int(round(value))
    text = str(value).strip().replace(",", "").replace(" ", "")
    if not text or text.lower() == "nan":
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def first_non_null(values: Sequence[object]) -> object | None:
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and np.isnan(value):
            continue
        if str(value).strip() == "":
            continue
        return value
    return None


def truthy_flag(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"x", "1", "true", "yes", "y"}


def dataframe_from_key_value_row(header: Sequence[object], row: Sequence[object]) -> pd.DataFrame:
    cols = [str(x).strip() for x in header]
    values = list(row)
    return pd.DataFrame({"field": cols, "value": values})


def rename_with_fallback(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    reverse = {str(k).strip().lower(): v for k, v in mapping.items()}
    renamed = {}
    for col in df.columns:
        key = str(col).strip().lower()
        if key in reverse:
            renamed[col] = reverse[key]
    return df.rename(columns=renamed)


def standardize_marker_frame(
    df: pd.DataFrame,
    *,
    dataset_id: str,
    dataset_name: str,
    parser_type: str,
    source_path: str,
    source_archive: str | None,
) -> pd.DataFrame:
    frame = df.copy()
    if "marker_id" not in frame.columns:
        raise ValueError("marker_id column is required in standardized marker frame.")
    if "chromosome" in frame.columns:
        frame["chromosome"] = frame["chromosome"].map(normalize_chromosome)
    if "position" in frame.columns:
        frame["position"] = frame["position"].map(coerce_integer)
    frame["marker_id"] = frame["marker_id"].astype(str).str.strip()
    frame["canonical_marker_id"] = frame["marker_id"].map(canonicalize_marker_id)
    frame.insert(0, "source_archive", source_archive)
    frame.insert(0, "source_path", source_path)
    frame.insert(0, "parser_type", parser_type)
    frame.insert(0, "dataset_name", dataset_name)
    frame.insert(0, "dataset_id", dataset_id)
    return frame


def standardize_sample_frame(
    df: pd.DataFrame,
    *,
    dataset_id: str,
    dataset_name: str,
    parser_type: str,
    source_path: str,
    source_archive: str | None,
) -> pd.DataFrame:
    frame = df.copy()
    if "sample_id" not in frame.columns:
        raise ValueError("sample_id column is required in standardized sample frame.")
    frame.insert(0, "source_archive", source_archive)
    frame.insert(0, "source_path", source_path)
    frame.insert(0, "parser_type", parser_type)
    frame.insert(0, "dataset_name", dataset_name)
    frame.insert(0, "dataset_id", dataset_id)
    return frame


def ordered_unique(values: Iterable[object]) -> list[object]:
    seen = set()
    output = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        output.append(value)
    return output
