# Developmental Enviromics for Wheat G×E Dissection

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

R scripts for processing weather data, APSIM phenology simulation, enviromic marker development, heritability analysis, GWAS, QTL/haplotype analysis, and genomic prediction for wheat genotype-by-environment interaction analysis.

> **Zhang C., He J., Yu R., et al.** A developmental enviromics framework dissects wheat thousand-kernel weight genotype-by-environment interactions into modular baseline potential and phenotypic plasticity. *(under review)*

---

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
└── 16_fw_validation/              # FW independent validation, cross-validation
```

## Workflow

Run scripts in numerical order (00 → 16) for full analysis pipeline.

| Phase | Scripts | Description |
|-------|---------|-------------|
| **1 — Data Preparation** | 00–05 | Weather data → APSIM phenology → enviromics alignment |
| **2 — Heritability & Variance** | 06–08 | FW independence test → variance decomposition → GW/G×W partitioning |
| **3 — Association & PCA** | 09–10 | Genotype/envirotype PCA → EWAS → GAPIT GWAS |
| **4 — Plasticity & QTL** | 11–13 | FW phenotypic plasticity → QTL detection → haplotype analysis → epistasis |
| **5 — Prediction & Figures** | 14–16 | Genomic prediction → manuscript figures → FW validation |

## Data

```
data/raw/
├── genotype/     983_renamed.vcf.gz   (25 MB, 983 wheat lines)
├── phenotype/    TKW.txt              (96 KB, thousand-kernel weight)
└── envirotype/   EC8.csv              (712 KB, 8 environmental covariates)
```

All processed data and analysis outputs are also deposited at **Figshare**:  
[https://doi.org/10.6084/m9.figshare.30873803](https://doi.org/10.6084/m9.figshare.30873803)

## Data Sources

- **NASA POWER** — historical daily weather data (1985–2025)
- **CMIP6** — future climate projections
- **APSIM Next Generation** — crop phenology simulation
- **Genotypic data** — 16K+5K targeted genotyping array (46,325 SNPs after QC)
- **Phenotypic data** — field trials across 8 environments (2024–2025)

## Requirements

- R ≥ 4.2.0
- APSIM Next Generation
- Key R packages: `GAPIT`, `sommer`, `BGLR`, `rrBLUP`, `tidyverse`, `data.table`, `ggplot2`

For exact package versions, see session info in the Figshare deposit.

## Citation

If you use this code or data, please cite the corresponding paper (forthcoming) and the Figshare deposit:

> Zhang C., He J., Yu R., et al. A developmental enviromics framework dissects wheat thousand-kernel weight genotype-by-environment interactions into modular baseline potential and phenotypic plasticity. *(under review)*

> Zhang C. et al. (2026). Analysis code and processed data for "A developmental enviromics framework dissects wheat TKW G×E interactions." Figshare. [https://doi.org/10.6084/m9.figshare.30873803](https://doi.org/10.6084/m9.figshare.30873803)

## License

MIT License — see [LICENSE](LICENSE).
