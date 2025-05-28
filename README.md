# Mendelian Randomization of Body Fat Distribution and Systemic Sclerosis

This repository contains all scripts and utilities used in the publication:

**Villanueva-Martin, G. & Bossini-Castillo, L.**  
*The Effect of Body Fat Distribution on Systemic Sclerosis*  
Journal of Clinical Medicine (2022), [DOI: 10.3390/jcm11206014](https://doi.org/10.3390/jcm11206014)

---

## 🧪 Objective

This project performs a cross-tissue Mendelian Randomization (MR) analysis to evaluate the causal role of body fat distribution traits (exposures) in the susceptibility to Systemic Sclerosis (SSc) (outcome). We use public summary statistics from the **GIANT consortium** and **GWAS Catalog**, applying both classical two-sample MR and MR-PRESSO for pleiotropy assessment.

---

## 🗂️ Repository structure

```bash
.
├── scripts/                   # Main MR pipeline scripts
│   ├── skeleton.sh           # Main entrypoint for MR analysis
│   ├── MR_2SMR_fun.R
│   ├── MR_Formating.R
│   ├── outcome_prep.sh
│
├── utilities/                # Helper scripts to prepare input GWAS files
│   ├── Formating_GWAS_Catalog.sh
│   ├── Formating_GWAS_UKBB_As_Outcome.sh
│   ├── MR_Clumping_Phenos.sh
│   ├── MR_Format_FDR.sh
│   ├── Merging_GWAS_Phenos.R
│
├── data/                     # Instructions to access required public GWAS datasets
│   └── README_data_sources.txt
│
├── results/                  # Output directory created by the pipeline
│
├── .gitignore
└── README.md                 # You are here
```

---

## 📥 Input data

All summary statistics used are publicly available:

- **SSc GWAS**: GWAS Catalog [GCST90253087](https://www.ebi.ac.uk/gwas/studies/GCST90253087)
- **Body fat traits (BMI, WHR, etc.)**: GIANT consortium, [https://portals.broadinstitute.org/collaboration/giant](https://portals.broadinstitute.org/collaboration/giant)

---

## ⚙️ Requirements

- **R ≥ 4.0** with packages:
  - `TwoSampleMR`, `MRPRESSO`, `ggplot2`, `optparse`, `dplyr`, `na.tools`, `tidyr`
- **Bash** and `awk`
- [`plink1.9`](https://www.cog-genomics.org/plink/)
- Unix tools: `sed`, `cat`, `sort`

---

## ▶️ How to run

Assuming you have properly formatted input files, run:

```bash
bash scripts/skeleton.sh \
    -o META_GWAS_SSc.tsv \
    -p obesity_traits_list.txt \
    -s Overall \
    -e OR
```

This will:
1. Format the outcome and exposure GWAS files
2. Perform clumping with `plink`
3. Merge summary statistics
4. Run MR + MR-PRESSO + sensitivity analyses
5. Generate plots and result tables

---

## 🧹 Utilities

If your input GWAS files need formatting, use the scripts in `utilities/`:
- `Formating_GWAS_Catalog.sh`: for GWAS Catalog-style files
- `Formating_GWAS_UKBB_As_Outcome.sh`: for UKBB outcomes
- `MR_Clumping_Phenos.sh`: full preprocessing and clumping
- `Merging_GWAS_Phenos.R`: merges summary statistics by shared SNPs

---

## 👤 Author

**Gonzalo Villanueva-Martin PhD**  
Instituto de Parasitologia y Biomedicina, Lopez-Neyra
Universidad de Granada
---

## 📄 License

This repository is for reproducibility and reference. No license is applied by default.
This README file was generated with the help of chatGPT
