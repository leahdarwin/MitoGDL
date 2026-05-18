# MitoGDL 
This repository contains all data and scripts used for the analyses and figures in the accompanying manuscript "Genetic and environmental interactions outweigh mitonuclear coevolution for complex traits in Drosophila". Data and scripts are organised into top-level `data/` and `scripts/` folders, each subdivided by trait. A separate `figures/` folder holds the output figures.

---

## Repository Structure

```
├── data/
│   ├── climb/          # Raw and adjusted climbing assay data
│   ├── development/    # Raw and adjusted development time data
│   ├── flight/         # Raw and adjusted flight performance data
│   ├── survival/       # Survival assay data
│   ├── weight/         # Raw and adjusted weight data
│   └── figures/        # Extra data for figures (phylogeny, CI chain SNPs, protein structure)
│
├── scripts/
│   ├── climb/          # Wrangling and ANOVA scripts for climbing assay
│   ├── development/    # Wrangling and ANOVA scripts for development time
│   ├── flight/         # Wrangling and ANOVA scripts for flight performance
│   ├── survival/       # ANOVA scripts for survival assay
│   ├── weight/         # Wrangling and ANOVA scripts for weight
│   └── figures/        # Scripts to reproduce all manuscript figures
│
├── figures/
│   ├── main_figs/      # Main manuscript figures
│   └── supp_figs/      # Supplemental figures
│
└── README.md
```

- **`data/`**  
  Contains all raw and processed data files, organised by trait. Each subfolder mirrors the trait structure used in the analysis scripts.  
  Complex I AA changes are given in `data/figures/CI_chain_snps.tsv`.

- **`scripts/`**  
  Contains all analysis and figure scripts, organised by trait. All scripts use paths relative to the project root (open `MitoGDL.Rproj` in RStudio to set the working directory automatically).

- **`figures/`**  
  Output figures written by the scripts in `scripts/figures/`.

---

## Reproducibility

Most analyses were conducted in R (version 4.4.0). The required R packages are listed and loaded within each script. To reproduce analyses or figures:  

1. Clone this repository.  
2. Open the `.Rproj` file in RStudio.  
3. Run the scripts in `scripts/<trait>/` to wrangle and adjust raw data.  
4. Use the scripts in `scripts/figures/` to generate the figures in the manuscript.

Protein figures were made using Python (version 3.12.12) with pymol-open-source (version 3.1.0).

---

## Citation

If you use this code or data, please cite the accompanying manuscript: 

Genetic and environmental interactions outweigh mitonuclear coevolution for complex traits in _Drosophila_
Leah J. Darwin, Faye A. Lemieux, Rebecca Z. Bachtel, Jack H. Blocker, Camille P. Brown, Jacob D. Lerman, Olivia C. Maule, Yevgeniy Raynes, David M. Rand
bioRxiv 2025.11.24.689096; doi: https://doi.org/10.1101/2025.11.24.689096 

---

## ⚠️ ggtree + aplot + ggplot2 Compatibility Note

**ggplot2 4.x is currently incompatible with ggtree 3.12.0 and aplot.**

ggplot2 4.0.0 removed internal functions (including `is.waive()` and `is_ggplot()`) that ggtree and aplot depend on, causing the following errors:

```
Error in `geom_segment2()`:
! Problem while converting geom to grob.
Caused by error in `is.waive()`:
! could not find function "is.waive"
```

```
Error: package or namespace load failed for 'aplot':
 object 'is_ggplot' is not exported by 'namespace:ggplot2'
```

### Workaround

Downgrade ggplot2 to the last stable 3.x release:

```r
remove.packages("ggplot2")
install.packages("remotes")
remotes::install_version("ggplot2", version = "3.5.1")
```

Restart R after installing.

This issue is expected to be resolved in a future ggtree update. Track progress on the [ggtree GitHub](https://github.com/YuLab-SMU/ggtree/issues).


