from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.analysis.overlaps import (
    build_association_overlap_report,
    build_candidate_gene_overlap_report,
    build_pairwise_overlap_report,
    build_pipeline_summary_report,
)
from cacao_aroma_pipeline.discovery import build_inventory_frame, collect_candidate_spreadsheets, scan_raw_files
from cacao_aroma_pipeline.excel import write_workbook
from cacao_aroma_pipeline.exceptions import DetectionError, ParsingError
from cacao_aroma_pipeline.models import ParsedDataset, RawFileRecord, RunContext, SupplementaryPayload
from cacao_aroma_pipeline.parsers import (
    detect_parser_type,
    parse_marker_metadata_workbook,
    parse_supplementary_files,
    parse_tassel_like_matrix,
    parse_tropgene_study_workbook,
)
from cacao_aroma_pipeline.utils import ensure_dir, slugify


def _load_parsed_datasets(candidates: list[RawFileRecord], logger: logging.Logger) -> tuple[list[ParsedDataset], SupplementaryPayload]:
    datasets: list[ParsedDataset] = []
    supplementary_candidates: list[Path] = []

    for candidate in candidates:
        try:
            parser_type = detect_parser_type(candidate.path)
        except DetectionError as exc:
            logger.warning("Skipping unsupported spreadsheet %s (%s)", candidate.relative_path, exc)
            continue

        logger.info("Parsing %s with parser %s", candidate.relative_path, parser_type)

        try:
            if parser_type.startswith("supplementary_"):
                supplementary_candidates.append(candidate.path)
                continue

            if parser_type == "tassel_like_matrix":
                dataset = parse_tassel_like_matrix(
                    candidate.path,
                    dataset_name=candidate.path.stem,
                    source_archive=candidate.source_archive,
                )
            elif parser_type == "tropgene_study_workbook":
                dataset = parse_tropgene_study_workbook(
                    candidate.path,
                    dataset_name=candidate.path.stem,
                    source_archive=candidate.source_archive,
                )
            elif parser_type == "marker_metadata_workbook":
                dataset = parse_marker_metadata_workbook(
                    candidate.path,
                    dataset_name=candidate.path.stem,
                    source_archive=candidate.source_archive,
                )
            else:
                logger.warning("Skipping unsupported parser output for %s", candidate.path)
                continue
        except Exception as exc:
            logger.exception("Failed to parse %s: %s", candidate.relative_path, exc)
            continue

        dataset.file_hash = candidate.md5
        datasets.append(dataset)

    supplementary = parse_supplementary_files(sorted(supplementary_candidates)) if supplementary_candidates else SupplementaryPayload(source_files=[])
    return datasets, supplementary


def _write_normalized_datasets(context: RunContext, datasets: list[ParsedDataset]) -> list[Path]:
    outputs = []
    max_cells = int(context.config["analysis"].get("max_cells_per_normalized_matrix_sheet", 250_000))
    preview_rows = int(context.config["analysis"].get("normalized_matrix_preview_rows", 1000))

    for dataset in datasets:
        workbook_path = context.normalized_dir / f"{dataset.dataset_id}.xlsx"
        sheets: dict[str, pd.DataFrame] = {
            "metadata": dataset.metadata,
            "markers": dataset.markers,
            "samples": dataset.samples,
        }
        if not dataset.genotype_matrix.empty:
            matrix = dataset.genotype_matrix.copy()
            matrix.insert(0, "marker_id", matrix.index)
            cell_count = int(matrix.shape[0] * max(matrix.shape[1], 1))
            matrix_meta = pd.DataFrame(
                [
                    {"metric": "rows", "value": int(matrix.shape[0])},
                    {"metric": "columns", "value": int(matrix.shape[1])},
                    {"metric": "cell_count", "value": cell_count},
                    {"metric": "full_matrix_written", "value": cell_count <= max_cells},
                    {"metric": "preview_rows", "value": min(preview_rows, int(matrix.shape[0]))},
                ]
            )
            sheets["genotype_matrix_metadata"] = matrix_meta
            if cell_count <= max_cells:
                sheets["genotype_matrix"] = matrix.reset_index(drop=True)
            else:
                sheets["genotype_matrix_preview"] = matrix.head(preview_rows).reset_index(drop=True)
                dataset.notes.append(
                    f"Full genotype matrix omitted from normalized workbook because estimated cell count ({cell_count}) exceeded configured limit ({max_cells})."
                )
        if not dataset.raw_marker_table.empty:
            sheets["raw_marker_table"] = dataset.raw_marker_table
        if not dataset.raw_sample_table.empty:
            sheets["raw_sample_table"] = dataset.raw_sample_table
        if dataset.notes:
            sheets["notes"] = pd.DataFrame({"note": dataset.notes})
        write_workbook(
            workbook_path,
            sheets,
            freeze_header=context.config["excel"]["freeze_header"],
            autofilter=context.config["excel"]["autofilter"],
            min_width=context.config["excel"]["min_column_width"],
            max_width=context.config["excel"]["max_column_width"],
        )
        outputs.append(workbook_path)
    return outputs


def _write_supplementary(context: RunContext, supplementary: SupplementaryPayload) -> Path | None:
    if not supplementary.source_files:
        return None
    output_dir = ensure_dir(context.processed_dir / "supplementary")
    workbook_path = output_dir / "fruity_aroma_supplementary.xlsx"
    sheets = {
        "sensory_traits": supplementary.sensory_traits,
        "associations_all": supplementary.associations,
        "association_summary": supplementary.association_summary,
        "candidate_genes": supplementary.candidate_genes,
        "notes": pd.DataFrame({"note": supplementary.notes or ["No parser notes."]}),
    }
    write_workbook(
        workbook_path,
        sheets,
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )
    return workbook_path




def _prepare_report_sheets(raw_sheets: dict[str, pd.DataFrame], context: RunContext) -> dict[str, pd.DataFrame]:
    if not context.config["analysis"].get("emit_detail_rows", True):
        return {name: df for name, df in raw_sheets.items() if "detail" not in name.lower()}

    max_rows = int(context.config["analysis"].get("max_detail_rows_per_sheet", 200_000))
    prepared: dict[str, pd.DataFrame] = {}
    truncation_rows = []
    for name, df in raw_sheets.items():
        if "detail" in name.lower() and len(df) > max_rows:
            prepared[name] = df.head(max_rows).copy()
            truncation_rows.append(
                {
                    "sheet_name": name,
                    "original_rows": len(df),
                    "written_rows": len(prepared[name]),
                    "truncated": True,
                }
            )
        else:
            prepared[name] = df
    if truncation_rows:
        prepared["detail_truncation"] = pd.DataFrame(truncation_rows)
    return prepared


def _write_master_marker_catalog(context: RunContext, datasets: list[ParsedDataset]) -> Path:
    output_dir = ensure_dir(context.processed_dir / "combined")
    workbook_path = output_dir / "master_marker_catalog.xlsx"
    marker_catalog = pd.concat([dataset.markers for dataset in datasets if not dataset.markers.empty], ignore_index=True)
    dataset_summary = pd.DataFrame(
        [
            {
                "dataset_id": dataset.dataset_id,
                "dataset_name": dataset.dataset_name,
                "marker_count": dataset.marker_count,
                "sample_count": dataset.sample_count,
                "parser_type": dataset.parser_type,
                "has_coordinates": dataset.has_coordinates,
                "has_genotypes": dataset.has_genotypes,
            }
            for dataset in datasets
        ]
    )
    write_workbook(
        workbook_path,
        {"marker_catalog": marker_catalog, "dataset_summary": dataset_summary},
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )
    return workbook_path


def run_pipeline(context: RunContext, logger: logging.Logger) -> dict[str, Path]:
    raw_records = scan_raw_files(
        context.raw_dir,
        allowed_extensions={ext.lower() for ext in context.config["discovery"]["file_extensions"]},
    )
    inventory = build_inventory_frame(raw_records)

    inventory_dir = ensure_dir(context.processed_dir / "inventory")
    inventory_path = inventory_dir / "raw_file_inventory.xlsx"
    write_workbook(
        inventory_path,
        {"raw_file_inventory": inventory},
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )

    candidates = collect_candidate_spreadsheets(
        raw_records,
        raw_dir=context.raw_dir,
        staging_dir=context.staging_dir,
        include_archives=bool(context.config["discovery"]["include_archives"]),
        archive_member_extensions={ext.lower() for ext in context.config["discovery"]["archive_member_extensions"]},
    )

    logger.info("Discovered %s raw files and %s candidate spreadsheets.", len(raw_records), len(candidates))
    datasets, supplementary = _load_parsed_datasets(candidates, logger)
    logger.info(
        "Loaded %s datasets | supplementary counts: sensory=%s, associations=%s, association_summary=%s, candidate_genes=%s",
        len(datasets),
        len(supplementary.sensory_traits),
        len(supplementary.associations),
        len(supplementary.association_summary),
        len(supplementary.candidate_genes),
    )

    normalized_outputs = _write_normalized_datasets(context, datasets)
    supplementary_output = _write_supplementary(context, supplementary)
    marker_catalog_output = _write_master_marker_catalog(context, datasets)

    reports = {}
    overlap_windows = list(context.config["analysis"]["overlap_windows_bp"])
    association_windows = list(context.config["analysis"]["association_windows_bp"])
    candidate_windows = list(context.config["analysis"]["candidate_windows_bp"])

    pairwise_report = _prepare_report_sheets(build_pairwise_overlap_report(datasets, overlap_windows), context)
    pairwise_report_path = context.report_dir / "pairwise_marker_overlap.xlsx"
    write_workbook(
        pairwise_report_path,
        pairwise_report,
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )
    reports["pairwise_marker_overlap"] = pairwise_report_path

    association_report = _prepare_report_sheets(build_association_overlap_report(supplementary, datasets, association_windows), context)
    association_report_path = context.report_dir / "association_marker_overlap.xlsx"
    write_workbook(
        association_report_path,
        association_report,
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )
    reports["association_marker_overlap"] = association_report_path

    candidate_report = _prepare_report_sheets(build_candidate_gene_overlap_report(supplementary, datasets, candidate_windows), context)
    candidate_report_path = context.report_dir / "candidate_gene_marker_overlap.xlsx"
    write_workbook(
        candidate_report_path,
        candidate_report,
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )
    reports["candidate_gene_marker_overlap"] = candidate_report_path

    summary_report = build_pipeline_summary_report(datasets, supplementary, inventory)
    summary_report_path = context.report_dir / "pipeline_summary.xlsx"
    write_workbook(
        summary_report_path,
        summary_report,
        freeze_header=context.config["excel"]["freeze_header"],
        autofilter=context.config["excel"]["autofilter"],
        min_width=context.config["excel"]["min_column_width"],
        max_width=context.config["excel"]["max_column_width"],
    )
    reports["pipeline_summary"] = summary_report_path

    outputs = {
        "inventory": inventory_path,
        "marker_catalog": marker_catalog_output,
        **reports,
    }
    if supplementary_output is not None:
        outputs["supplementary"] = supplementary_output
    for idx, path in enumerate(normalized_outputs, start=1):
        outputs[f"normalized_{idx}"] = path
    return outputs


def validate_project(context: RunContext) -> pd.DataFrame:
    checks = []
    checks.append({"check": "project_root_exists", "status": context.project_root.exists(), "detail": str(context.project_root)})
    checks.append({"check": "raw_dir_exists", "status": context.raw_dir.exists(), "detail": str(context.raw_dir)})
    checks.append({"check": "processed_dir_exists", "status": context.processed_dir.exists(), "detail": str(context.processed_dir)})
    checks.append({"check": "results_dir_exists", "status": context.results_dir.exists(), "detail": str(context.results_dir)})
    if context.raw_dir.exists():
        spreadsheet_count = sum(1 for path in context.raw_dir.rglob("*") if path.is_file() and path.suffix.lower() in {".xlsx", ".xls", ".zip"})
        checks.append({"check": "raw_spreadsheet_or_archive_count_gt_zero", "status": spreadsheet_count > 0, "detail": str(spreadsheet_count)})
    return pd.DataFrame(checks)
