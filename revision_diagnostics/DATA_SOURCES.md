# Input sources and roles

The input root contains the paths below. Exact bytes are checked against the
published-study manifests bundled with each scientific stage.

| Relative path | Role |
|---|---|
| `nacional_paper/mmc2_all_associations.xlsx` | Published association records and originating-filter fields |
| `nacional_paper/mmc4_candidate_genes.xlsx` | Published candidate annotations |
| `amazonia/Cocoa Amazonia Ecuador aroma.xlsx` | Original-coordinate resource context |
| `tropgene/genotypes/cocoa_diversity/Cocoa diversity.xlsx` | Diversity marker sequences and identifiers |
| `tropgene/genotypes/cocoa_nacional_aroma/GBS_Nac_G7.xls` | Nacional G7 genotype export |
| `tropgene/genotypes/cocoa_nacional_aroma/GBS_Nac_MAF5.xls` | Nacional MAF5 genotype export |
| `tropgene/genotypes/incompatibility_snp_mapping_brazil/Incompatibility SNP mapping Brazil.xlsx` | Auxiliary source identity audit |
| `tropgene/genotypes/reference_collection_for_aDNA_studies/Reference collection for aDNA studies.xlsx` | Auxiliary source identity audit |

Obtain the association and candidate workbooks from Colonges et al. (2022),
Plant Physiology and Biochemistry 171:213–225,
https://doi.org/10.1016/j.plaphy.2021.11.006 . Original publisher filenames
`1-s2.0-S0981942821005647-mmc2.xlsx` and `1-s2.0-S0981942821005647-mmc4.xlsx`
correspond to the first two files; use the paths above for the companion.

Obtain the genotype/marker exports from the TropGeneDB COCOA studies-download
interface, https://tropgenedb.cirad.fr/tropgene/JSP/interface.jsp?module=COCOA .
Select the named studies and retain the original workbook contents. The source
filenames use `GBS`; the relevant study describes DArTseq genotyping. The code
preserves those immutable filenames rather than renaming provider exports.

The harmonization, genotype and adjusted-LD stages verify the eight-workbook
study source set. The background stage uses only the complete published S1
table workbook. The analysis does not require phenotype data or sequencing
reads, and these are not inferred from genotype calls.

Large public tables are supplied separately in Supplementary Data S1, version
1.0.0, https://doi.org/10.5281/zenodo.22755806 . It contains the fixed reference
annotation, marker mapping decisions, coordinate cohorts and original summaries
needed by the additional calculations. Its input role is explicit: the companion
recomputes downstream analyses from these fixed tables and source workbooks.
The separate version 3.0.0 pipeline starts from source workbooks and the genome
reference when a new reference-alignment run is required.
