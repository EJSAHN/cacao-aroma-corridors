# Revision code snapshot

This repository supports **Coordinate conventions and marker density shape coverage
of cacao aroma loci**. It contains two separately versioned analysis components:

| Component | Version | Entry point | Inputs |
|---|---|---|---|
| Original marker harmonization | 3.0.0 | `scripts/run_analysis.py` | Original provider workbooks and pinned Criollo v2 reference |
| Revision diagnostics | 1.0.0 | `revision_diagnostics/run_revision.py` | Complete published Supplementary Data S1 and original provider workbooks |

The revision adds the diagnostics directory and documentation. The original
scientific source, configuration and tests retain their submitted bytes. The
original source is also preserved at commit
`67c8235362982b0ba1b6612eea24f914666e5290`; that historical tree contains all
submitted source files except the archive's `.gitignore`. The current tree restores
that file and extends the root README. No original calculation is replaced.

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
(ESM 6). This repository does not contain the manuscript, response letter, rendered
figures or large analysis outputs.

## Verification

Before this update, all 72 files in the original GitHub tree matched the submitted
archive by canonical Git blob identity. The original pipeline's 349 tests and the
diagnostics' 45 tests were rerun successfully. All 38 scientific source records in
the diagnostics manifest match the verified author-run packages.

The four accepted native Windows analysis stages passed 1,586 checks. The separate
portable replay compared 110 tables containing 1,274,097 rows with those accepted
results: 98 tables were byte-identical; the remaining adjusted-correlation tables
agreed within absolute tolerance 1e-8 and relative tolerance 1e-7, with no cells
outside tolerance. Details are in
[`revision_diagnostics/VALIDATION.json`](../revision_diagnostics/VALIDATION.json).
These checks do not claim a new full genome-alignment run.

For source integrity, the companion's `PACKAGE_SHA256.json` records every other
companion file, and `SCIENTIFIC_SOURCE_MANIFEST.json` links scientific files to the
verified revision stages. Record the full Git commit when using this repository;
`main` may advance after this snapshot.

## Citation and archive identifiers

The root `CITATION.cff` describes the original 3.0.0 software; the diagnostics have
their own `CITATION.cff` and version. The dataset identifier
`10.5281/zenodo.22755806` refers to Supplementary Data S1, not to a software archive.
An archival software DOI must be cited only after its record is published. This
code snapshot does not assign a DOI or claim that a GitHub release has been created.
