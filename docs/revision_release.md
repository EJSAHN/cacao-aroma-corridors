# Software versions and reproducibility

This repository supports **Coordinate conventions and marker density shape coverage
of cacao aroma loci**. It contains two separately versioned analysis components:

| Component | Version | Entry point | Inputs |
|---|---|---|---|
| Original marker harmonization | 3.0.0 | `scripts/run_analysis.py` | Original provider workbooks and pinned Criollo v2 reference |
| Revision diagnostics | 1.0.0 | `revision_diagnostics/run_revision.py` | Complete published Supplementary Data S1 and original provider workbooks |

The original scientific source, configuration and tests are preserved at commit
`67c8235362982b0ba1b6612eea24f914666e5290`. The revision diagnostics provide additional
analyses using the original fixed tables.

## Reproducing the study

Follow the [original pipeline instructions](../README.md) to start from source
workbooks and the reference genome. Follow the
[diagnostics instructions](../revision_diagnostics/README.md) to reproduce the
additional analyses from the published fixed tables. These are different input
routes: the diagnostics do not rerun genome alignment or reconstruct missing
discovery-variant alleles.

The exact original source archive and the diagnostics archive are preserved in
Online Resource 2 (ESM 2). Their SHA-256 digests are recorded in
[`code_archives.sha256`](code_archives.sha256). Third-party source workbooks are
obtained from their providers and are not included in this repository.

Figure and publication production is distributed separately in Online Resource 6
(ESM 6).

## Verification

The original pipeline passed 349 tests and the diagnostics passed 45 tests.
The diagnostics source manifest records the relationship between the public
files and the original analysis packages.

The four accepted native Windows analysis stages passed 1,586 checks. The separate
portable replay compared 110 tables containing 1,274,097 rows with those accepted
results: 98 tables were byte-identical; the remaining adjusted-correlation tables
agreed within absolute tolerance 1e-8 and relative tolerance 1e-7, with no cells
outside tolerance. Details are in
[`revision_diagnostics/VALIDATION.json`](../revision_diagnostics/VALIDATION.json).
The portable validation covers downstream analyses; it does not repeat genome
alignment.

For source integrity, the companion's `PACKAGE_SHA256.json` records every other
companion file, and `SCIENTIFIC_SOURCE_MANIFEST.json` links scientific files to the
original analysis stages. Record the full Git commit when using this repository;
`main` may advance as documentation or software is updated.

## Citation and archive identifiers

The root `CITATION.cff` describes the original 3.0.0 software; the diagnostics have
their own `CITATION.cff` and version. The dataset identifier
`10.5281/zenodo.22755806` refers to Supplementary Data S1, not to a software archive.
Use the software version and full Git commit to identify the code used.
