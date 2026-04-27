from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.models import SupplementaryPayload
from cacao_aroma_pipeline.utils import coerce_integer, truthy_flag


def _read_table(path: Path, header_row: int = 1, sheet_name: str | None = None) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet_name, header=header_row)


def _first_cell_text(path: Path) -> str:
    try:
        probe = pd.read_excel(path, sheet_name=0, header=None, nrows=1)
    except Exception:
        return ""
    if probe.empty:
        return ""
    return str(probe.iloc[0, 0]).strip().lower()


def _supplementary_slot(path: Path) -> int | None:
    stem = path.stem.lower()
    match = re.search(r"mmc([1-7])", stem)
    if match:
        return int(match.group(1))

    a1 = _first_cell_text(path)
    if a1.startswith("table s1"):
        return 1
    if a1.startswith("table s2"):
        return 2
    if a1.startswith("table s3"):
        return 3
    if a1.startswith("table s4"):
        return 4
    return None


def _normalize_numeric_string(value: object) -> object:
    if pd.isna(value):
        return value
    text = str(value).strip()
    if not text:
        return pd.NA
    if re.fullmatch(r"[\d\s]+", text):
        return text.replace(" ", "")
    return value


def _standardize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(col).strip().replace("\n", " ") for col in df.columns]
    return df


def _normalize_associations(path: Path) -> pd.DataFrame:
    excel = pd.ExcelFile(path)
    frames = []
    for sheet_name in excel.sheet_names:
        # mmc2 sheets vary: some need header=1, others header=0.
        header_candidates = [1, 0]
        frame = None
        for header_row in header_candidates:
            candidate = _read_table(path, header_row=header_row, sheet_name=sheet_name)
            candidate = _standardize_column_names(candidate)
            if {"Chromosome", "Traits"}.issubset(set(candidate.columns)) or {"Chromosome", "Position of the association peak"}.issubset(set(candidate.columns)):
                frame = candidate
                break
        if frame is None:
            continue

        frame = frame.rename(
            columns={
                "Chromosome": "chromosome",
                "Position of the association peak (bp)": "peak_position_bp",
                "Position of the association peak": "peak_position_bp",
                "Position of haplotypic bloc start": "haploblock_start_bp",
                "Position of haplotypic bloc end": "haploblock_end_bp",
                "N° haplotypic bloc": "haploblock_number",
                "N° hap. Bloc": "haploblock_number",
                "GWAS method": "gwas_method",
                "GWAS method ": "gwas_method",
                "Sorting of marker": "marker_sorting",
                "Sorting of marker ": "marker_sorting",
                "Traits": "trait",
                "p-value of the strongest association": "p_value",
                "Explanation rate of the trait of the strongest association": "trait_explained_rate",
                "Associations detected": "associations_detected",
            }
        )

        if "chromosome" not in frame.columns:
            continue
        if "trait" not in frame.columns:
            continue

        for col in ["peak_position_bp", "haploblock_start_bp", "haploblock_end_bp", "haploblock_number", "associations_detected"]:
            if col in frame.columns:
                frame[col] = frame[col].map(_normalize_numeric_string).map(coerce_integer)
        frame["chromosome"] = frame["chromosome"].map(_normalize_numeric_string).map(coerce_integer)
        if "p_value" in frame.columns:
            frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
        if "trait_explained_rate" in frame.columns:
            frame["trait_explained_rate"] = pd.to_numeric(frame["trait_explained_rate"], errors="coerce")

        frame = frame.dropna(subset=["chromosome", "trait"], how="all").copy()
        frame["source_sheet"] = sheet_name
        frame["source_file"] = path.name
        frames.append(frame)

    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    combined.insert(0, "association_id", range(1, len(combined) + 1))
    return combined


def _normalize_candidate_genes(path: Path) -> pd.DataFrame:
    frame = _read_table(path, header_row=1, sheet_name=None)
    if isinstance(frame, dict):
        frame = next(iter(frame.values()))
    frame = _standardize_column_names(frame)
    frame = frame.rename(
        columns={
            "Chromosome": "chromosome",
            "start": "start_bp",
            "end": "end_bp",
            "gene_id": "gene_id",
            "gene_function": "gene_function",
            "Monoterpene pathway": "pathway_monoterpene",
            "L-phenylalanine degradation pathway": "pathway_l_phe",
            "Fatty acid pathway": "pathway_fatty_acid",
            "Sugar pathway": "pathway_sugar",
            "Maillard precursor": "pathway_maillard",
            "Defense (L-phe)": "defense_l_phe",
            "Defense (FA/SS)": "defense_fa_ss",
            "Defense (pyr)": "defense_pyrazine",
            "Implication in defense": "defense_note",
        }
    )
    for col in ["chromosome", "start_bp", "end_bp"]:
        if col in frame.columns:
            frame[col] = frame[col].map(_normalize_numeric_string).map(coerce_integer)
    for col in [
        "pathway_monoterpene",
        "pathway_l_phe",
        "pathway_fatty_acid",
        "pathway_sugar",
        "pathway_maillard",
        "defense_l_phe",
        "defense_fa_ss",
        "defense_pyrazine",
    ]:
        if col in frame.columns:
            frame[col] = frame[col].map(truthy_flag)
    frame = frame.dropna(subset=["gene_id"], how="all").copy()
    frame.insert(0, "candidate_gene_id", range(1, len(frame) + 1))
    return frame


def _normalize_sensory_traits(path: Path) -> pd.DataFrame:
    frame = _read_table(path, header_row=1, sheet_name=None)
    if isinstance(frame, dict):
        frame = next(iter(frame.values()))
    frame = _standardize_column_names(frame)
    if "Trait" not in frame.columns:
        first_col = frame.columns[0]
        frame = frame.rename(columns={first_col: "Trait"})
    frame = frame.rename(columns={"Trait": "trait"})
    frame = frame.dropna(subset=["trait"]).copy()
    frame.insert(0, "trait_id", range(1, len(frame) + 1))
    return frame


def parse_supplementary_files(files: list[Path]) -> SupplementaryPayload:
    payload = SupplementaryPayload(source_files=[str(path) for path in files])
    for path in files:
        slot = _supplementary_slot(path)
        if slot == 1:
            payload.sensory_traits = _normalize_sensory_traits(path)
        elif slot == 2:
            payload.associations = _normalize_associations(path)
        elif slot == 3:
            payload.association_summary = _normalize_associations(path)
        elif slot == 4:
            payload.candidate_genes = _normalize_candidate_genes(path)
        else:
            payload.notes.append(f"Skipped supplementary file in core parser: {path.name}")

    payload.notes.append(
        "Supplementary row counts | sensory_traits="
        f"{len(payload.sensory_traits)} | associations={len(payload.associations)} | "
        f"association_summary={len(payload.association_summary)} | candidate_genes={len(payload.candidate_genes)}"
    )
    return payload
