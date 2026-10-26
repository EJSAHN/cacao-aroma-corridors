# Cacao revision diagnostics

Version 1.0.0 reproduces the additional analyses of marker provenance, coordinate
harmonization, genotype concordance, adjusted genotype correlations, peak-level
coordinate anchors, and matched annotated-gene backgrounds.

These analyses support the revised manuscript **Coordinate conventions and marker
density shape coverage of cacao aroma loci**. The original submitted manuscript was
**Sequence-based harmonization of cacao marker resources clarifies coverage of
aroma-associated loci**. The earlier dataset title uses **Coordinate harmonization
and marker density shape the representation of cacao aroma loci**. These titles
identify successive descriptions of the same study, not separate datasets.

The original **cacao-marker-harmonization 3.0.0** pipeline remains the source-to-
reference analysis. This companion uses its published fixed tables and the original
workbooks to reproduce the additional analyses. The scientific calculation modules
are unchanged from the four verified revision stages; their SHA-256 identities are
listed in `SCIENTIFIC_SOURCE_MANIFEST.json`. The portable entry point accepts explicit
paths and does not require a particular drive, username, installation or old run folder.

## Requirements and inputs

Use Python 3.10 or newer and NumPy 1.26 or newer. The included read-only XLSX and
restricted BIFF8 readers handle the source workbook layouts. No packages are installed
and no data are downloaded by this companion. Genome alignment is not repeated here.

Obtain the complete **Supplementary Data S1, version 1.0.0** workbook identified by
the original study at https://doi.org/10.5281/zenodo.22755806. Its SHA-256 is
`f80b33927090d11b804a7dffd09333945a0fecf58243f44aa8bbbdc9dac16fbb`.
The `Contents` worksheet maps original output workbooks and worksheets to the public
catalog. The displayed workbook filename may retain an earlier journal abbreviation;
content identity is determined by the digest and catalog, not that abbreviation.

Place unedited source workbooks at the relative paths in `DATA_SOURCES.md`. Original
provider rights remain separate from this code's MIT license. Source workbooks and
the large reference table workbook are not redistributed in this code archive.

## Run

Run all four stages from the extracted companion directory:

```text
python -B run_revision.py --stage all --supplement-s1 "<Supplementary_Data_S1.xlsx>" --source-root "<data_raw>" --output-root "<new_results_directory>"
```

Alternatively choose one stage: `harmonization`, `genotypes`, `adjusted_ld`, or
`background`. The `background` stage requires only Supplementary Data S1 and may
omit `--source-root`. Keep outputs separate from source data and code. Each stage
creates a new directory and ZIP with tables, checks, summaries and manifests.
All source data and scientific modules are read-only. A failed check returns a
nonzero exit code and retains diagnostic information.

For machines where BLAS otherwise starts many threads, setting the process variables
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `OMP_NUM_THREADS` to `1` can reduce
overhead in the small adjusted-correlation matrices. This is a runtime choice, not
a scientific parameter.

The original version 3.0.0 pipeline can reproduce source parsing and reference
alignment from the provider workbooks and pinned Criollo v2 reference. Its separate
README gives that workflow. A newly generated workbook archive can differ in file
bytes even when scientific tables agree. This companion intentionally checks the
published S1 reference for exact study reproduction; it does not silently substitute
a newly generated or altered reference table.

## Stages and outputs

| Stage | Scope | Principal outputs |
|---|---|---|
| Harmonization | Source-filter denominators, own-filter recovery, sequence mapping versus identifier links, mapped strand, quotas and tag lengths | Filter and recovery summaries, per-peak coverage, identifier-link tables, mapping-rejection evidence |
| Genotypes | Literal source calls, between-export differences, fixed marker cohorts and same-tag sensitivity | Genotype and sample QC, discordant calls, unadjusted pair and neighbor summaries |
| Adjusted LD | Leave-one-chromosome-out PC and kinship sensitivity; shared-sample cross-export coordinate anchors | Adjusted distance profiles, neighbor summaries, peak assessment and model diagnostics |
| Background | Original-model replay and length/density/radius sensitivity | Model status, draw counts, balance, percentiles and matching strata |

The module directories preserve their original file names to retain traceability.
`CODE1` through `CODE4` labels in the source manifest identify the verified author
run packages; they are not additional software dependencies.

## Interpretation

G7 and MAF5 are overlapping exports, not independent populations. Dosages use the
conditional A/M/C interpretation documented in `modules/genotypes/CODING_EVIDENCE.txt`.
The exact exporter codebook was not confirmed. A positional anchor does not prove
the identity of the original discovery variant, and a missing focal genotype is
not evidence of low correlation. See the adjusted-math notes for the covariance
and nuisance-projection definitions.

Background percentiles describe the specified sampled annotation distributions;
they are not enrichment P-values. The full grid includes insufficient-stratum
models, which remain blocked without dropping candidates or relaxing matching.
The original tested-GWAS universe is not reconstructed by annotation matching.

## Testing and validation scope

Run the companion interface tests:

```text
python -B -m unittest discover -s tests -p "test_*.py"
```

The scientific tests are available in their respective module directories. Native
Windows results from the four original stages were verified before packaging.
The portable entry point is then checked against those result tables. Validation
describes the exact run scope in `VALIDATION.json`; it is not a claim that genome
alignment was rerun by this companion.

## Versions and citation

The original analysis software is version **3.0.0**; this companion is version
**1.0.0**; the original figure-generation package is a separate version **1.0.0**.
Those numbers refer to different software components. The repository name
`cacao-aroma-corridors` is a stable project identifier, not a second analysis.
The earlier dataset title is retained in the provenance of the original data
record and refers to the same study and fixed table export.

Code is MIT licensed. Cite the original data sources as well as the software.
Archival DOI metadata must be taken from an actual published archival record;
this code does not mint or predict a DOI.
