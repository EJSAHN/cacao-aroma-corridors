# Data sources and expected input layout

This repository does not redistribute raw or compiled third-party input spreadsheets. Users should obtain the relevant files from the original public sources and place them under `data_raw/`.

## Primary sources

- Published Nacional fruity-aroma GWAS and supplementary tables: Colonges et al. (2022), *Plant Physiology and Biochemistry*.
- TropGeneDB COCOA module: public cacao marker and study exports.

## Recommended local layout

The pipeline uses structure and spreadsheet content to detect supported files, but the layout below is recommended for reproducibility.

```text
data_raw/
  amazonia/
    Cocoa Amazonia Ecuador aroma.xlsx
  nacional_paper/
    mmc1_sensory_traits.xlsx
    mmc2_all_associations.xlsx
    mmc3_pathway_summary.xlsx
    mmc4_candidate_genes.xlsx
  tropgene/
    genotypes/
      cocoa_diversity/
        Cocoa diversity.xlsx
      cocoa_nacional_aroma/
        GBS_Nac_G7.xls
        GBS_Nac_MAF5.xls
      incompatibility_snp_mapping_brazil/
        Incompatibility SNP mapping Brazil.xlsx
      reference_collection_for_aDNA_studies/
        Reference collection for aDNA studies.xlsx
```

Archive files such as `Cocoa Nacional aroma.zip` may also be placed in `data_raw/`; archive scanning is enabled by default in `pipeline_config.yaml`.

## Notes for users

- Input names do not need to be hardcoded in the pipeline, but recognizable source names make inventory review easier.
- The Cocoa diversity resource is treated as a marker scaffold when no accession-level genotype matrix is detected.
- Optional auxiliary files may appear in the inventory but are not required for the core association/candidate-gene portability outputs.
