# Developmental Enviromics for Wheat G×E Dissection

R scripts for processing historical and future weather data, APSIM phenology simulation, and enviromic marker development for wheat genotype-by-environment interaction analysis.

## Repository Structure

```
scripts/
├── 00_data_download/           # CMIP6/NASA weather data download & GDD calculation
├── 01_met_file_creation/        # APSIM .met file creation & weather segmentation
├── 02_phenology_apsim/          # Phenology data processing & APSIM combined outputs
├── 03_growthstage_alignment/    # Growth stage alignment reports
├── 04_enviromics_merge/         # Weather-phenology merge & enviromic marker development
└── 05_visualization/            # Weather data visualization
```

## Workflow

Run scripts in numerical order (00 → 06) within each stage.

## Data Sources

- NASA POWER (historical weather)
- CMIP6 (future climate projections)
- APSIM Next Generation (phenology simulation)

## Requirements

- R >= 4.2.0
- APSIM Next Generation
