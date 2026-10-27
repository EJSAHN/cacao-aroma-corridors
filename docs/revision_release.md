# Software versions and reproducibility

This repository supports **Coordinate conventions and marker density shape coverage
of cacao aroma loci**. It contains two analysis components:

| Component | Version | Entry point | Inputs |
|---|---|---|---|
| Marker harmonization | 3.0.0 | `scripts/run_analysis.py` | Provider workbooks and pinned Criollo v2 reference |
| Revision diagnostics | 1.0.0 | `revision_diagnostics/run_revision.py` | Complete Supplementary Data S1 and provider workbooks |

## Reproducing the study

Select the release cited by the study from
[GitHub releases](https://github.com/EJSAHN/cacao-aroma-corridors/releases) or
[Zenodo](https://doi.org/10.5281/zenodo.23677321).
Record its version DOI and full commit identifier when reporting an analysis.

Follow the [pipeline instructions](../README.md) to start from provider workbooks
and the reference genome. Follow the
[diagnostics instructions](../revision_diagnostics/README.md) for calculations
using fixed analysis tables. The diagnostics do not repeat genome alignment or
reconstruct missing discovery-variant alleles.

Supplementary Data S1, version 1.0.0, is available at
https://doi.org/10.5281/zenodo.22755806. Obtain third-party source workbooks from
their providers as described in [data sources](data_sources.md).

## Verification

The archived pipeline passed 349 tests and the diagnostics passed 45 tests.
The reference analyses passed 1,586 checks. A comparison of 110 tables containing
1,274,097 rows found 98 byte-identical tables. The remaining adjusted-correlation
tables agreed within absolute tolerance 1e-8 and relative tolerance 1e-7, with
no cells outside tolerance. These comparisons cover downstream calculations.
Detailed records are in
[`revision_diagnostics/VALIDATION.json`](../revision_diagnostics/VALIDATION.json).

`PACKAGE_SHA256.json` records file digests for the diagnostics component.
`SCIENTIFIC_SOURCE_MANIFEST.json` links calculation files to their source
packages and retains source digests. The checksums of the archived component
ZIPs are listed in [`code_archives.sha256`](code_archives.sha256).

## Citation

The concept DOI https://doi.org/10.5281/zenodo.23677321 groups the combined
repository releases. Each archived release has its own version DOI.
The component versions identify marker harmonization 3.0.0
and revision diagnostics 1.0.0. Each component includes a `CITATION.cff` file.
The dataset DOI https://doi.org/10.5281/zenodo.22755806 identifies Supplementary
Data S1. Cite the software, data and original provider studies used in an analysis.
