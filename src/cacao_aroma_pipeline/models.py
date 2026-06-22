from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(slots=True)
class RawFileRecord:
    path: Path
    relative_path: str
    extension: str
    size_bytes: int
    md5: str
    source_archive: str | None = None
    extracted_member: str | None = None


@dataclass(slots=True)
class ParsedDataset:
    dataset_id: str
    dataset_name: str
    parser_type: str
    source_path: str
    source_archive: str | None
    file_hash: str
    markers: pd.DataFrame = field(default_factory=pd.DataFrame)
    samples: pd.DataFrame = field(default_factory=pd.DataFrame)
    genotype_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: pd.DataFrame = field(default_factory=pd.DataFrame)
    raw_marker_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    raw_sample_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    notes: list[str] = field(default_factory=list)

    @property
    def has_genotypes(self) -> bool:
        return not self.genotype_matrix.empty

    @property
    def has_coordinates(self) -> bool:
        if self.markers.empty:
            return False
        required = {"chromosome", "position"}
        if not required.issubset(set(self.markers.columns)):
            return False
        return bool(self.markers[["chromosome", "position"]].dropna(how="any").shape[0] > 0)

    @property
    def marker_count(self) -> int:
        return 0 if self.markers.empty else int(self.markers["marker_id"].nunique(dropna=True))

    @property
    def sample_count(self) -> int:
        if not self.samples.empty and "sample_id" in self.samples.columns:
            return int(self.samples["sample_id"].nunique(dropna=True))
        if not self.genotype_matrix.empty:
            return int(len(self.genotype_matrix.columns))
        return 0


@dataclass(slots=True)
class SupplementaryPayload:
    source_files: list[str]
    sensory_traits: pd.DataFrame = field(default_factory=pd.DataFrame)
    associations: pd.DataFrame = field(default_factory=pd.DataFrame)
    association_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    candidate_genes: pd.DataFrame = field(default_factory=pd.DataFrame)
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class RunContext:
    project_root: Path
    raw_dir: Path
    processed_dir: Path
    results_dir: Path
    staging_dir: Path
    report_dir: Path
    normalized_dir: Path
    log_dir: Path
    config: dict[str, Any]
