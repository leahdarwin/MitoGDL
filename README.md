# MitoGDL Project Overview

This repository contains all data and scripts used for the analyses and figures in the accompanying manuscript "Genetic and environmental interactions outweigh mitonuclear coevolution for complex traits in Drosophila". The project is organized by trait, with each trait-specific folder containing raw data and scripts for processing. A separate `figures` folder contains the code required to reproduce all figures presented in the manuscript.  

---

## Repository Structure

```
├── development/
├── weight/
├── climb/
├── flight/
│   ├── data/ # Raw data files for this trait
│   ├── scripts/  # Data wrangling and adjustment scripts
│   └── ...
│
├── figures/
│   ├── scripts/ # Scripts to reproduce all manuscript figures
│   └── ...
│
└── README.md
```

- **Trait folders** (`development`, `weight`, `climb`, `flight`)  
  Each folder contains:  
  - **Raw trait data** (`data/`): Direct exports from experiments and adjusted data.  
  - **Wrangling scripts** (`scripts/`): Code to clean and adjust data following procedures described in the manuscript.  

- **Figures folder**  
  Contains all R scripts needed to reproduce the figures from the manuscript using the processed data.  

---

## Reproducibility

All analyses were conducted in R (version 4.4.0). The required R packages are listed and loaded within each script. To reproduce analyses or figures:  

1. Clone this repository.  
2. Open the `.Rproj` file in RStudio.  
3. Run the scripts in the trait folders to wrangle and adjust raw data.  
4. Use the scripts in `figures/` to generate the figures in the manuscript.  

---

## Citation

If you use this code or data, please cite the accompanying manuscript:  

*Author(s), Title, Year, Journal/Preprint DOI*  

