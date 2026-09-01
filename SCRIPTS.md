# Scripts

Descriptions below are copied from each script's header comment.

## scripts/climb/

- **climb_adj.R** — Apply within-set and across-set correction terms to climbing assay data and save adjusted dataset.
- **climb_ANOVA.R** — Phenotype Analysis: Climbing Assay with Mixed Models.
  - Loads pre-processed climbing assay data
  - Fits linear mixed-effects models (LMMs) using lmerTest
  - Runs Type III ANOVAs
  - Computes R² and effect sizes
  - Generates publication-ready tables
- **wrangle_climbSet.R** — Process climbing assay data by cleaning, merging, and standardizing genotype annotations to create a reproducible dataset for downstream analyses.

## scripts/development/

- **development_adj.R** — Normalize pupation assay data by applying within-set and across-set correction terms, adjusting for larval density, and exporting a cleaned dataset for downstream analyses.
- **development_anova.R** — Fit linear mixed models to adjusted pupation (development) data, evaluate fixed effects, compute effect sizes, and format results for LaTeX-compatible tables.
- **larvalDensity_anova.R** — *(no header comment in this file — it fits `Larval_density ~ Mito * Nuc * Treatment` as a mixed model and runs `anova()` on it; consider adding a header comment.)*
- **wrangle_pupationSet.R** — Process pupation assay data by parsing raw CSVs, extracting genotype and metadata annotations, computing time-since-egg (TSE), and merging with stock genotypes to produce a clean, standardized dataset for analysis.

## scripts/figures/

- **AMMI_plot.R** — Perform AMMI analysis on adjusted climbing data and generate genotype-by-environment interaction biplots stratified by sex. Generates figure 2.
- **buildCorr_plot.R** — Calculate correlations between builds for each trait.
- **CIcontact_plot.py** — Use pymol to highlight contact sites for mtDNA and nuclear proteins given a csv of AA changes. (Dependencies: pymol, pandas)
- **CIcontact_sourceData.R** — Standalone (no PyMOL) reproduction of the mito/nuclear contact-interface computation done in CIcontact_plot.py, for exporting the Fig 5 source-data table. PyMOL is expensive to run, so this parses the mmCIF structure directly and computes the same "within 4.0 Angstrom" contact test as a nearest-neighbor distance check.
- **compile_source_data.R** — Compile every CSV under source_data/{main_figs,supp_figs}/ into one Excel workbook (source_data/SourceData.xlsx), one tab per CSV/panel, main figures first then supplementary, in manuscript figure-number order. Run this AFTER all figure scripts (and CIcontact_sourceData.R) have populated source_data/.
- **ixn_stability.R** — Tests stability of interaction terms via subsampling (without replacement) across n_mito levels c(5, 10, 15, 20) with 1000 draws each.
- **maineffects_plot.R** — Plots emmeans and adjusted averages for significant main effects.
- **mothersCurse_plot.R** — Analyze fly trait datasets (climb, flight, weight) to assess variance, coefficient of variation (CV), sex differences, and correlations. Generates bootstrapped plots, variance/correlation tables, and combined summary figures.
- **phenoCorr_plot.R** — Correlation analysis of phenotypic traits (weight, climb, flight, development) for males and females across nuclear backgrounds.
- **phenoPCA_plot.R** — Summarize fly phenotype data (weight, climb, flight, development) by sex, treatment, and genotype, merge into a combined dataset, and run PCA to visualize multivariate trait patterns.
- **phyloSignal_plot.R** — Test for phylogenetic signal (Blomberg's K) in fly traits across mitochondrial haplotypes. Data include climbing, flight, development, and weight phenotypes.
- **secondOrderInteraction_plot.R** — Plots interaction plots for significant second order effects using estimated marginal means.
- **set9Corr_plot.R** — Compare split vs. unified experiments across climbing, weight, and development phenotypes. Produces rank–rank correlation plots and stats.
- **sourceData_helpers.R** — Shared helper for exporting the tabular data behind each figure panel/table to source_data/, for the journal-required Source Data workbook. Sourced by every script in scripts/figures/.
- **surv_plot.R** — Generate boxplots for adult survival counts from egg-pick and egg-count assays, stratified by sex.
- **weightDensity_plot.R** — Plot weight vs. estimated larval density (Females/Males) with a linear fit and R²/p annotation. The mixed-model test of Larval_density ~ Mito * Nuc * Treatment lives in scripts/development/larvalDensity_anova.R; this script only builds and saves the figure.

## scripts/flight/

- **flight_adj.R** — Apply within-set and across-set corrections to raw flight assay data, generating adjusted phenotypes for downstream analysis and plotting.
- **flight_ANOVA.R** — Perform mixed-model ANOVAs on adjusted flight phenotypes, stratified by sex, and extract effect sizes and R² values for reporting.
- **wrangle_flightSet.R** — Description:
  - Reads in raw flight assay CSV files
  - Cleans and harmonizes column formats
  - Assigns Set, Build, and Treatment values
  - Joins genotype metadata from stock_genotype.csv
  - Outputs final cleaned dataset: flight.csv

## scripts/survival/

- **surv_ANOVA.R** — Perform ANOVAs on adult survival counts and produce separate LaTeX ANOVA and estimated marginal means tables via kable booktabs.

## scripts/weight/

- **weight_adj.R** — Process fly weight data by removing outliers, applying within- and across-set corrections, and generating an adjusted dataset for downstream analyses.
- **weight_ANOVA.R** — Fit linear mixed models to adjusted fly weight data, evaluate fixed effects, compute effect sizes, and format results for LaTeX-compatible tables.
- **wrangle_weightSet.R** — Process fly body weight data by merging multiple sets, correcting annotations, computing per-fly weight, and standardizing genotype labels.
