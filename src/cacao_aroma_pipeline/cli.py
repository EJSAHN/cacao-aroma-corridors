from __future__ import annotations

import argparse
from pathlib import Path

from cacao_aroma_pipeline.config import build_run_context, load_config
from cacao_aroma_pipeline.excel import write_workbook
from cacao_aroma_pipeline.logging_utils import setup_logging
from cacao_aroma_pipeline.pipeline import run_pipeline, validate_project


def _base_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cacao-pipeline", description="Spreadsheet-first cacao aroma pipeline.")
    parser.add_argument("--project-root", type=Path, default=Path("."), help="Path to the project root.")
    parser.add_argument("--config", type=Path, default=None, help="Optional YAML config path.")
    return parser


def _build_parser() -> argparse.ArgumentParser:
    parser = _base_parser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("run", help="Run the full pipeline.")
    subparsers.add_parser("inventory", help="Run the full pipeline inventory and normalization workflow.")
    subparsers.add_parser("validate", help="Validate project structure and write a validation workbook.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    project_root = args.project_root.resolve()
    config = load_config(project_root, args.config.resolve() if args.config else None)
    context = build_run_context(project_root, config)
    logger = setup_logging(context.log_dir, args.command)

    logger.info("Project root: %s", project_root)
    logger.info("Command: %s", args.command)

    if args.command in {"run", "inventory"}:
        outputs = run_pipeline(context, logger)
        for name, path in outputs.items():
            logger.info("Wrote %s -> %s", name, path)
        return 0

    if args.command == "validate":
        validation = validate_project(context)
        output_path = context.report_dir / "validation.xlsx"
        write_workbook(
            output_path,
            {"validation": validation},
            freeze_header=context.config["excel"]["freeze_header"],
            autofilter=context.config["excel"]["autofilter"],
            min_width=context.config["excel"]["min_column_width"],
            max_width=context.config["excel"]["max_column_width"],
        )
        logger.info("Wrote validation report -> %s", output_path)
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
