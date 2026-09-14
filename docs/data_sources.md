# Data sources

Obtain third-party data from their original providers and retain their licenses and
study metadata. The listed names are the original resource labels, not inferred
publication-to-workbook equivalences. No phenotype measurements or sequencing reads
are fabricated or downloaded by this workflow.

## Required input workbooks

Paths below are relative to the selected source root (`data_raw` when present).
The exact expected path is preferred; filename aliases are tried only when the
expected path is absent, and ambiguous copies cause an error. The configuration
lists permitted aliases, worksheet names and expected panel roles.

| Relative path | Resource and role |
|---|---|
| `nacional_paper/mmc2_all_associations.xlsx` | Colonges et al. supplementary association worksheets; master `synth_tot` and its trait subsets |
| `nacional_paper/mmc4_candidate_genes.xlsx` | Colonges et al. candidate annotations; `Synth` |
| `amazonia/Cocoa Amazonia Ecuador aroma.xlsx` | TropGeneDB Amazonia genotype resource; original-coordinate context only |
| `tropgene/genotypes/cocoa_diversity/Cocoa diversity.xlsx` | TropGeneDB diversity marker/flank/reference-sequence resource |
| `tropgene/genotypes/cocoa_nacional_aroma/GBS_Nac_G7.xls` | TropGeneDB Nacional G7 filter export |
| `tropgene/genotypes/cocoa_nacional_aroma/GBS_Nac_MAF5.xls` | TropGeneDB Nacional MAF5 filter export |

Optional ID/resource context files are
`tropgene/genotypes/incompatibility_snp_mapping_brazil/Incompatibility SNP mapping Brazil.xlsx`
and
`tropgene/genotypes/reference_collection_for_aDNA_studies/Reference collection for aDNA studies.xlsx`.
Their unavailable physical coordinates are not inferred. Omitting optional resources
changes the ID/resource audit but not the configured principal coverage/LD analyses.

Download the original supplementary files from the publisher record below. Original
publisher names `1-s2.0-S0981942821005647-mmc2.xlsx` and
`1-s2.0-S0981942821005647-mmc4.xlsx` are supported filename aliases. From the TropGeneDB
COCOA studies-download interface, select the named Genotypes studies and retain the
original workbook contents. Extract the Nacional archive to expose both XLS files.
Place each unedited workbook under the relative path above; copying inputs into a
second analysis directory is unnecessary. Record an actual acquisition date locally,
not a date inferred from the database study year.

- Colonges et al. (2022), Plant Physiology and Biochemistry 171:213–225.
  DOI: https://doi.org/10.1016/j.plaphy.2021.11.006
- TropGeneDB COCOA module:
  https://tropgenedb.cirad.fr/tropgene/JSP/interface.jsp?module=COCOA
- Hamelin et al. (2013), Nucleic Acids Research 41:D1172–D1175.
  DOI: https://doi.org/10.1093/nar/gks1105

The workflow records relative input paths, file sizes and SHA-256 digests. The
absence of an accession-level matrix in the diversity export is not a claim that
the originating project had no biological samples. G7/MAF5 sample overlap is
measured directly instead of assuming independent populations.

## Reference assembly and annotation

Use the complete Criollo v2 toplevel reference, not primary chromosomes alone.
Unplaced contigs must participate in mapping ambiguity checks. Their placements
are retained in evidence tables but not included in primary-chromosome coverage.

- Assembly: `Criollo_cocoa_genome_V2`, accession `GCA_000208745.2`.
- Ensembl Plants release: **62**.
- FASTA: `Theobroma_cacao_criollo.Criollo_cocoa_genome_V2.dna.toplevel.fa.gz`.
- GFF3: `Theobroma_cacao_criollo.Criollo_cocoa_genome_V2.62.gff3.gz`.
- Argout et al. (2017), BMC Genomics 18:730.
  DOI: https://doi.org/10.1186/s12864-017-4120-9

Exact HTTPS download locations are declared in `config/study.json`. Local files
are validated and reused without copying when provided explicitly. SHA-256 is
recorded as observed provenance, not misrepresented as a publisher-authenticated
checksum. An explicit expected digest can be supplied in the configuration.

## Native aligner

Bowtie 2 2.5.5 binaries are available at
https://github.com/BenLangmead/bowtie2/releases/tag/v2.5.5 . Pinned asset digests and
sizes are recorded in the configuration. An existing executable directory is
version checked and its binary hashes recorded. Complete reference/index/SAM cache
identity includes reference data, query FASTA, binary hashes and alignment settings.
