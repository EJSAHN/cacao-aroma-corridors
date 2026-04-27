# Reproducibility checklist

1. Create or activate a Python environment with the dependencies in `environment.yml`.
2. Install the repository with `pip install -e .`.
3. Obtain the raw input spreadsheets from the primary sources listed in `docs/DATA_SOURCES.md`.
4. Place input files under `data_raw/` using the recommended layout.
5. Run:

```bash
python scripts/run_pipeline.py --project-root . run
```

6. Confirm that the following files are generated:

```text
results/reports/pairwise_marker_overlap.xlsx
results/reports/association_marker_overlap.xlsx
results/reports/candidate_gene_marker_overlap.xlsx
results/reports/pipeline_summary.xlsx
```

7. Use `results/logs/run.log` and `data_processed/inventory/raw_file_inventory.xlsx` to audit parser decisions.
