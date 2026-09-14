# Output schema

Each run is isolated under `analysis_<timestamp>`. Directory names identify analysis
functions, not separate versioned dependencies. All numerical output is XLSX. Logs,
configuration snapshots, hashes and completion files accompany the workbooks.

| Directory | Contents |
|---|---|
| `source_data` | Resources, original association/candidate rows and distinct units, sample overlap, marker identity/coordinate audit |
| `source_evidence` | Tag/identifier correspondence, coordinate hypotheses, annotation-record crosswalk, source genotype concordance |
| `reference_mapping` | Exact full-tag mapping, projected coordinate records, coordinate conventions, reference gene catalog, unchanged association context |
| `exact_policy_coverage` | Fixed-cohort coverage under the exact-only mapping policy |
| `alignment_coverage` | Seven-function evidence in six workbooks: bounded alignment decisions; marker cohorts; query sets; paired coverage and exact-policy sensitivity; count-matched sensitivity; matched-gene background |
| `dosage_ld` | Genotype/coordinate QC; distance profiles; full pair tables; fixed shared-collection comparisons; observed-SNP neighbor summaries |
| `provenance` | Inputs, reference assets, software versions, statuses and limitations |
| `validation` | Optional externally referenced numerical agreement tables and summary |

## Main interpretation files

`alignment_coverage/04_paired_coverage_effects.xlsx` separates A/B/C changes.
`alignment_coverage/05_count_matched_sensitivity.xlsx` reports quotas and all draws.
`alignment_coverage/06_matched_gene_background.xlsx` records backgrounds, balance,
deficient strata and all descriptive draw summaries.
`dosage_ld/02_ld_distance_profiles.xlsx` and
`dosage_ld/04_shared_collection_comparison.xlsx` summarize dosage evidence.

`RUN_SUMMARY.txt` is a reading guide, not an approval of genetic-effect transfer.
`analysis_summary.json` and `OUTPUT_SHA256.json` allow completion and integrity
checks. `source_allele_tags.fa` contains the generated biallelic query tags.
`alignment_provenance.json` identifies the SAM and aligner configuration. Complete
SAM/index/FASTA files reside in the explicitly chosen cache and are not duplicated
into output ZIPs. Local execution paths, when needed to reuse local software, are
kept in `local_execution_paths.json` beside run directories rather than in tables.

Counts of marker records, unique positions, genes, source annotation rows, shared
samples, queried pairs and background replicates have different denominators. They
must not be relabeled as interchangeable measures. Missing numerical values are
not converted to zero. Source strings resembling formulas are stored as text.
The output preserves original `synth_tot` and other worksheet names as provenance.

A `BLOCKED_INSUFFICIENT_STRATA` background model is not a failed execution. It
means no complete comparable draw was possible under the stated conditions.
The full candidate result and coding-only sensitivity remain separate. Neither
receives a fabricated significance label or a replacement target denominator.
