# Cacao marker harmonization

Version 3.0.0. Reproducible analysis of marker identity, reference-derived SNP
coordinates, positional representation and conditional genotype-dosage correlations
in public cacao resources. Numerical deliverables are values-only Excel workbooks.
No figures, phenotype predictions, new GWAS, imputed genotypes or reconstructed
haplotypes are generated.

## Scope

The workflow starts with the original spreadsheets and the declared Criollo v2
FASTA/GFF3. It does **not** require existing result workbooks or previous run folders.
Original files are read-only. Source annotation records, analysis units and reference
positions are distinct objects; source peaks and gene boundaries are never shifted
merely to agree with target markers.

1. Parse source records and identify overlapping summary worksheets, unique peak
   positions, source-defined intervals and repeated candidate annotations.
2. Inspect marker/reference-sequence names, source tag sequences, coordinate
   conventions, genotype symbols and shared sample calls.
3. Search complete biallelic tags exactly in both orientations across the entire
   pinned reference and compare source candidate intervals with the gene catalog.
4. Reassess every eligible tag with pinned Bowtie 2 end-to-end alignment. Exclude
   ambiguous, censored, gap-proximal and allele-discordant focal positions.
5. Compare A (all eligible original markers), B (retained marker records with
   original coordinates) and C (the **same B records**, reference-derived SNP
   coordinates). Report mapping-policy sensitivity separately.
6. Perform chromosome-matched marker-count resampling and descriptive matched-gene
   background comparisons for the complete candidate set and the coding-only set.
7. Calculate conditional unphased dosage r2 on fixed Nacional marker cohorts,
   including prespecified sample sensitivity and exactly matched shared collections.

G7 and MAF5 are overlapping filter exports, not independent populations. Amazonia
is retained as original-coordinate context; its coding and sequence bridge are not
inferred. High positional coverage is not evidence of aroma-effect transferability.

## Requirements

Python 3.10 or newer; NumPy >=1.26, pandas >=2.1 and openpyxl >=3.1. An installed
xlrd >=2.0.1 is used for legacy XLS when available. A restricted, read-only BIFF8
fallback is included for supported value-only layouts and rejects unsupported
records instead of silently interpreting them. No Excel formulas are evaluated.

Bowtie 2 **2.5.5** is the only external executable. Supply an existing directory
with `bowtie2-align-s` and `bowtie2-build-s` (with `.exe` on Windows), or allow
verified official binaries to be downloaded to the specified tools directory.
Windows x86-64 and Linux x86-64 native assets are configured. The Python wrapper
invokes the native binaries directly; no Perl, WSL or Conda modification is needed.

The reference is `Criollo_cocoa_genome_V2`, `GCA_000208745.2`, Ensembl Plants 62.
Supply the complete pinned toplevel FASTA and GFF3 or allow download to the explicit
reference cache. Whole-genome length, sequence inventory, assembly declaration and
candidate annotation are checked. See [data sources](docs/data_sources.md).

## Run from the source tree

No editable installation is necessary. The entry point resolves its own `src`
directory and accepts explicit data and output paths. Replace angle-bracket values
below with local paths, retaining quotes around paths containing spaces.

```text
python -B scripts/run_analysis.py selftest
python -B scripts/run_analysis.py preflight --source-root "<source-project>"
python -B scripts/run_analysis.py run --source-root "<source-project>" --output-root "<output-directory>" --scratch "<temporary-directory>" --cache-root "<alignment-cache>" --tools-root "<native-tools>" --reference-cache "<reference-cache>" --reference-fasta "<toplevel.fa.gz>" --reference-gff "<annotation.gff3.gz>" --bowtie2-dir "<bowtie2-directory>" --offline
```

`--source-root` may be a project containing `data_raw` or that directory itself.
Output, scratch, native-tool and cache locations must be separate, nonnested paths
outside source inputs and code. The scientific configuration is `config/study.json`;
`--config` can select a different file. Defaults are the published analysis policy,
not optimized parameters or inferred biological assumptions. Changing a policy
requires interpreting and documenting a new analysis, not overriding integrity checks.

`--offline` prohibits downloads. Without it, only the configured reference and
pinned native-tool assets may be downloaded. No Python packages are installed.
Use `--threads 1` to `--threads 8` to bound native work; four threads is the reference
setting. Parallel aligner settings can change candidate reporting and are recorded.
The complete local pair inventory is processed in bounded batches, not replaced by
a sparse list of high-r2 pairs.

## Outputs and preservation

Every run creates a new `analysis_<timestamp>` directory with scientific workbooks,
configuration, summary, logs and a SHA-256 manifest, followed by a result ZIP.
Inputs, FASTA/GFF3, executables, index files and complete SAM are **not** copied into
the result archive. Verified reference/index/alignment caches are reusable; an
incomplete or changed cache is not silently accepted. Exact-search results are
recomputed from the original reference rather than loaded from an old workbook.
Full enumerated pair tables are retained and split at an indexed worksheet row
limit when necessary. There are no preview-only replacements for full tables.

`COMPLETED_WITH_DECLARED_LIMITATIONS` means execution completed; it does not certify
biological portability. A background model can be `BLOCKED_INSUFFICIENT_STRATA`;
its criteria are never silently weakened. See [output schema](docs/output_schema.md),
[analysis definitions](docs/methods.md) and [validation](docs/validation.md).

The separately distributed local convenience launcher deploys **this same source
archive** and may discover existing local caches. It is not a second scientific
implementation and is not needed by other users of this repository.

## Numerical reference comparison

A new raw-input run can optionally be compared against externally preserved
coverage and dosage workbooks. These directories are validation references only,
never required scientific inputs to the raw-input workflow.

```text
python -B scripts/run_analysis.py verify-results --analysis-run "<completed-analysis>" --reference-coverage "<coverage-reference-directory>" --reference-ld "<dosage-reference-directory>"
```

A failed comparison returns a nonzero exit code and lists differing fields.
Do not accept a release by deleting those discrepancies or changing reference data.
No old result directories or local comparison fixtures are distributed here.

## License and citation

The code is MIT licensed; source-data and external-tool licenses remain separate.
Third-party input spreadsheets are not redistributed. Retain source identifiers,
original worksheet names and acquisition details when reusing results. Software
citation metadata is provided in `CITATION.cff`.
