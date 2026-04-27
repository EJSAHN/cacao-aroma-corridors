from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from cacao_aroma_pipeline.models import RunContext
from cacao_aroma_pipeline.utils import ensure_dir


DEFAULT_CONFIG: dict[str, Any] = {
    "project": {
        "raw_dir": "data_raw",
        "processed_dir": "data_processed",
        "results_dir": "results",
        "log_subdir": "logs",
        "normalized_subdir": "normalized",
        "report_subdir": "reports",
        "staging_subdir": "staging",
    },
    "analysis": {
        "overlap_windows_bp": [0, 10_000, 50_000, 250_000],
        "association_windows_bp": [0, 10_000, 50_000, 250_000],
        "candidate_windows_bp": [0, 50_000, 250_000],
        "emit_detail_rows": True,
        "max_detail_rows_per_sheet": 200_000,
        "max_cells_per_normalized_matrix_sheet": 250_000,
        "normalized_matrix_preview_rows": 1000,
    },
    "discovery": {
        "include_archives": True,
        "archive_member_extensions": [".xlsx", ".xls"],
        "file_extensions": [".xlsx", ".xls", ".zip", ".pdf"],
    },
    "excel": {
        "freeze_header": True,
        "autofilter": True,
        "max_column_width": 60,
        "min_column_width": 10,
    },
}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    output = dict(base)
    for key, value in override.items():
        if key in output and isinstance(output[key], dict) and isinstance(value, dict):
            output[key] = _deep_merge(output[key], value)
        else:
            output[key] = value
    return output


def load_config(project_root: Path, config_path: Path | None = None) -> dict[str, Any]:
    resolved = DEFAULT_CONFIG
    candidate_paths = []
    if config_path is not None:
        candidate_paths.append(config_path)
    else:
        candidate_paths.extend(
            [
                project_root / "pipeline_config.yaml",
                project_root / "pipeline_config.yml",
            ]
        )
    for candidate in candidate_paths:
        if candidate.exists():
            user_config = yaml.safe_load(candidate.read_text(encoding="utf-8")) or {}
            resolved = _deep_merge(resolved, user_config)
            break
    return resolved


def build_run_context(project_root: Path, config: dict[str, Any]) -> RunContext:
    project_settings = config["project"]
    raw_dir = project_root / project_settings["raw_dir"]
    processed_dir = ensure_dir(project_root / project_settings["processed_dir"])
    results_dir = ensure_dir(project_root / project_settings["results_dir"])
    log_dir = ensure_dir(results_dir / project_settings["log_subdir"])
    normalized_dir = ensure_dir(processed_dir / project_settings["normalized_subdir"])
    report_dir = ensure_dir(results_dir / project_settings["report_subdir"])
    staging_dir = ensure_dir(processed_dir / project_settings["staging_subdir"])
    return RunContext(
        project_root=project_root,
        raw_dir=raw_dir,
        processed_dir=processed_dir,
        results_dir=results_dir,
        staging_dir=staging_dir,
        report_dir=report_dir,
        normalized_dir=normalized_dir,
        log_dir=log_dir,
        config=config,
    )
