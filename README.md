# Developmental enviromics disentangles baseline potential and phenotypic plasticity in wheat grain weight

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

R scripts for processing weather data, APSIM phenology simulation, enviromic marker development, heritability analysis, GWAS, QTL/haplotype analysis, and genomic prediction for wheat genotype-by-environment interaction analysis.

> **Zhang C., He J., Yu R., et al.** *Developmental enviromics disentangles baseline potential and phenotypic plasticity in wheat grain weight.* (under review)

---

## Repository Structure

```
scripts/
├── 00_data_download/              # NASA POWER retrieval & GDD calculation
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

## Marker sets — 46,325 vs 43,481

These two numbers refer to different sets and should not be confused:

| Number | What it is |
|---|---|
| **46,325** | Marker **backbone** after quality control (read depth ≥ 10, MAF ≥ 0.01, heterozygosity rate ≤ 10%, missing rate ≤ 10%). Used for genomic relationship matrix construction, variance-component estimation and genomic prediction. |
| **43,481** | The subset of the backbone actually **tested in the GWAS** (Figure 4). |
| **3,215** | Number of **independent tests**, obtained by partitioning the 43,481 GWAS SNPs into LD blocks (GAPIT v4 `GPART`, r² threshold = 0.7). LD-adjusted Bonferroni threshold: P = 0.05 / 3,215 = 1.56 × 10⁻⁵. |

## Environmental windows

Developmental windows are defined on anthesis-anchored APSIM developmental coordinates. The
exhaustive sliding-window scan spans 10 stages before to 30 stages after anthesis (window lengths
5–30 stages, step 1 stage), giving **637 windows** (21 pre-flowering, 239 cross-flowering,
377 post-flowering). Multiple testing was controlled by the Benjamini–Hochberg FDR across all
window × covariate combinations tested (9,872 tests).

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

- **NASA POWER** — historical daily weather data (1985–2025). This is the meteorological dataset
  underlying every analysis reported in the manuscript.
- **CMIP6** — future climate projections. Scripts for CMIP6 processing are retained in
  `00_data_download/` for completeness; **no CMIP6-derived result is reported in the manuscript**.
- **APSIM Next Generation 2024.10.7600.0** — crop phenology simulation, using the CAMP
  (Cereal Anthesis Molecular Phenology) wheat phenology model.
- **Genotypic data** — 16K+5K targeted genotyping array (46,325 SNPs after QC).
- **Phenotypic data** — field trials across 8 environments (2024–2025).

## Requirements

- R ≥ 4.2.0
- APSIM Next Generation 2024.10.7600.0 (CAMP wheat phenology model)
- Key R packages: `GAPIT`, `sommer`, `BGLR`, `rrBLUP`, `tidyverse`, `data.table`, `ggplot2`

For exact package versions, see session info in the Figshare deposit.

## Citation

If you use this code or data, please cite the corresponding paper and the Figshare deposit:

> Zhang C., He J., Yu R., et al. *Developmental enviromics disentangles baseline potential and phenotypic plasticity in wheat grain weight.* (under review)

> Zhang C. et al. (2026). Analysis code and processed data for "Developmental enviromics disentangles baseline potential and phenotypic plasticity in wheat grain weight." Figshare. [https://doi.org/10.6084/m9.figshare.30873803](https://doi.org/10.6084/m9.figshare.30873803)

## License

MIT License — see [LICENSE](LICENSE).
