# Developmental Enviromics for Wheat G×E Dissection

R scripts for processing weather data, APSIM phenology simulation, enviromic marker development, heritability analysis, GWAS, QTL/haplotype analysis, and genomic prediction for wheat genotype-by-environment interaction analysis.

## Repository Structure

```
scripts/
├── 00_data_download/              # CMIP6/NASA weather data download & GDD calculation
├── 01_met_file_creation/          # APSIM .met file creation & weather segmentation
├── 02_phenology_apsim/            # Phenology data processing & APSIM combined outputs
├── 03_growthstage_alignment/      # Growth stage alignment reports
├── 04_enviromics_merge/           # Weather-phenology merge & enviromic marker development
├── 05_visualization/              # Weather data visualization
├── 06_heritability_variance/      # Heritability: FW independence, variance partitioning
├── 07_genotype_prep/              # VCF to 012 genotype matrix conversion
├── 08_variance_decomposition/     # GW & GW-G×W model variance partitioning
├── 09_pca_analysis/               # Genotype PCA, envirotype PCA, combined analysis
├── 10_ewas_association/           # EWAS individual/population, GAPIT GWAS
├── 11_phenotypic_plasticity/      # FW model cross-validation, envirotype factor, GWAS
├── 12_qtl_haplotype/              # QTL aggregation, haplotype classification, effect trend
├── 13_epistasis_haplotype/        # Epistasis, combined haplotype, QTL interaction
├── 14_genomic_prediction/         # Genomic selection validation (3 strategies)
├── 15_final_figures/              # Combined manuscript figures
└── 16_fw_validation/              # FW independent validation, 3-fold cross-validation
```

## Workflow

Run scripts in numerical order (00 → 16) for full analysis pipeline.

**Phase 1 — Data Preparation (00–05):** Weather data → APSIM phenology → enviromics alignment

**Phase 2 — Heritability & Variance (06–08):** FW independence test → variance decomposition → GW/G×W partitioning

**Phase 3 — Association & PCA (09–10):** Genotype/envirotype PCA → EWAS → GAPIT GWAS

**Phase 4 — Plasticity & QTL (11–13):** FW phenotypic plasticity → QTL detection → haplotype analysis → epistasis

**Phase 5 — Prediction & Figures (14–16):** Genomic prediction → manuscript figures → FW validation

## Data Sources

- NASA POWER (historical weather)
- CMIP6 (future climate projections)
- APSIM Next Generation (phenology simulation)
- Genotypic data (VCF)
- Phenotypic data (field trials)

## Requirements

- R >= 4.2.0
- APSIM Next Generation
- Key packages: GAPIT, sommer, BGLR, rrBLUP, tidyverse
