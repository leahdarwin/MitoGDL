# MitoGDL 
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
│   ├── scripts/ # Scripts to reproduce all figures
│   ├── main_figs/ # main manucript figures
│   ├── ed_figs/ # extended data figures
│   ├── supp_figs/ # supplemental data figures
│   ├── extra_data/ # Extra data related to creation of figures but not data from the primary experiment 
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

  Complex I AA changes are given in `/figures/extra_data/CI_chain_snps.tsv`.

---

## Reproducibility

Most analyses were conducted in R (version 4.4.0). The required R packages are listed and loaded within each script. To reproduce analyses or figures:  

1. Clone this repository.  
2. Open the `.Rproj` file in RStudio.  
3. Run the scripts in the trait folders to wrangle and adjust raw data.  
4. Use the scripts in `figures/` to generate the figures in the manuscript.

Protein figures were made using Python (version 3.12.12) with pymol-open-source (version 3.1.0).

---

## Citation

If you use this code or data, please cite the accompanying manuscript: 

Genetic and environmental interactions outweigh mitonuclear coevolution for complex traits in _Drosophila_
Leah J. Darwin, Faye A. Lemieux, Rebecca Z. Bachtel, Jack H. Blocker, Camille P. Brown, Jacob D. Lerman, Olivia C. Maule, Yevgeniy Raynes, David M. Rand
bioRxiv 2025.11.24.689096; doi: https://doi.org/10.1101/2025.11.24.689096 


