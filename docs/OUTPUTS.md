# Output workbooks

The pipeline writes Excel workbooks under `data_processed/` and `results/`.

## `data_processed/inventory/raw_file_inventory.xlsx`
Inventory of discovered raw files, supported formats, parser decisions, and notes.

## `data_processed/normalized/*.xlsx`
Per-dataset normalized marker or genotype summaries. Large genotype matrices are summarized to keep workbooks practical.

## `data_processed/combined/master_marker_catalog.xlsx`
Combined marker catalog with dataset labels, marker identifiers, canonicalized marker identifiers, chromosome, position, and coordinate availability.

## `data_processed/supplementary/fruity_aroma_supplementary.xlsx`
Parsed reference supplementary layers, including sensory traits, association records, association summaries, and candidate genes.

## `results/reports/pairwise_marker_overlap.xlsx`
Pairwise marker continuity across datasets using raw/canonical marker IDs and coordinate-based windows.

## `results/reports/association_marker_overlap.xlsx`
Recovery of published fruity-aroma association zones across target marker panels under exact and flanking-window criteria.

## `results/reports/candidate_gene_marker_overlap.xlsx`
Candidate-gene and pathway-level corridor recovery across marker panels.

## `results/reports/pipeline_summary.xlsx`
Run-level summaries for datasets, supplementary parsing, and raw input inventory.
