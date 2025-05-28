# Multivariable Mendelian Randomization (MVMR) Module

This folder contains all scripts used to run a multivariable Mendelian Randomization (MVMR) pipeline, extending the core Mendelian randomization analysis to include multiple exposures simultaneously.

---

## 📘 Description

This module is part of the broader project:  
**"The Effect of Body Fat Distribution on Systemic Sclerosis"**  
[DOI: 10.3390/jcm11206014](https://doi.org/10.3390/jcm11206014)

It allows for testing the *independent effects* of multiple genetically correlated traits (e.g., BMI, WHR) on a single outcome (Systemic Sclerosis), correcting for confounding between exposures.

---

## 📁 Folder structure

```bash
MVMR/
├── MRMV_skeleton.sh          # Main pipeline script
├── MRMV_fun_v2.R             # MVMR logic and visualization
├── PheWAS_clump.sh           # Clumps and formats exposures
├── PheWAS_prep.py            # Prepares NEALE exposures
├── outcome_prep.sh           # Formats outcome GWAS
```

---

## ▶️ How to use

### 1. Prepare input files:
- One outcome GWAS (formatted with `outcome_prep.sh`)
- A list of PheWAS/trait files (`trait_list.txt`)
- Optionally, primary/secondary SNP signal files per trait

### 2. Run the full pipeline:

```bash
bash MRMV_skeleton.sh \
    -o META_GWAS_SSc.tsv \
    -p trait_list.txt \
    -i none \
    -k none \
    -d MVMR_run_01
```

This will:
- Format and clump each trait
- Merge SNPs
- Prune LD with PLINK
- Format for multivariable MR
- Run `MRMV_fun_v2.R` and generate plots and tables

---

## 📦 Requirements

- `plink1.9`
- R packages:
  - `TwoSampleMR`, `ggplot2`, `dplyr`, `optparse`, `na.tools`
- Python ≥3.6 (`pandas` for `PheWAS_prep.py`)

---

## 📄 Output

- Pruned and harmonized exposures
- Scatter plots for each exposure
- `rawresults_...txt`: raw MVMR results
- `table_...txt`: summary table with OR, 95% CI, and p-values

---

## 👤 Author

**Gonzalo Villanueva-Martin, PhD**  
Instituto de Parasitologia y Biomedicina, Lopez-Neyra
Universidad de Granada
---
