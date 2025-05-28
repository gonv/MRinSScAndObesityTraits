#!/usr/bin/env Rscript

# === Multivariable Mendelian Randomization Analysis ===
# This script harmonizes and runs MVMR between multiple exposures and a single outcome.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(devtools)
  library(TwoSampleMR)
  library(na.tools)
})

args <- commandArgs(trailingOnly = TRUE)
outcome_file <- args[1]
prune_list <- args[2]
tag <- args[3]

# === Load outcome GWAS ===
outcome <- read.table(outcome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "unknown")
outcome_exp <- format_data(
  outcome, type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", pval_col = "pval",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  samplesize_col = "samplesize"
)
outcome_exp$id.outcome <- outcome_file

# === Load exposure files ===
filenames <- read.table(prune_list, header = FALSE, stringsAsFactors = FALSE)[[1]]
exposures <- lapply(filenames, function(f) {
  read_outcome_data(f,
    sep = "\t",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    pval_col = "pval", samplesize_col = "samplesize",
    id_col = "id"
  )
})

exposure_dat <- bind_rows(exposures)
id_exposure <- unique(exposure_dat$id.outcome)

# Convert first to exposure, rest stay as outcomes
d1 <- convert_outcome_to_exposure(subset(exposure_dat, id.outcome == id_exposure[1]))
d2 <- subset(exposure_dat, id.outcome != id_exposure[1])

# Harmonize exposures against each other
d <- harmonise_data(d1, d2, action = 2)

# Reshape exposures into MVMR input format
dh1 <- subset(d, id.exposure == id.exposure[1],
              select = c(SNP, exposure, id.exposure, effect_allele.exposure,
                         other_allele.exposure, eaf.exposure, beta.exposure,
                         se.exposure, pval.exposure))
dh2 <- subset(d, id.outcome != id.exposure[1],
              select = c(SNP, outcome, id.outcome, effect_allele.outcome,
                         other_allele.outcome, eaf.outcome, beta.outcome,
                         se.outcome, pval.outcome))
colnames(dh2) <- gsub("outcome", "exposure", colnames(dh2))
dh <- rbind(dh1, dh2)
dh$exposure <- "exposure"

# Harmonize against outcome
mvdat <- mv_harmonise_data(dh, outcome_exp, harmonise_strictness = 1)

# Run MVMR analysis
resMV <- mv_multiple(mvdat, plots = TRUE, pval_threshold = 1,
                     intercept = FALSE, instrument_specific = FALSE)
ggsave("MRMV_plots_expo1.png", plot = resMV$plots[[1]])
ggsave("MRMV_plots_expo2.png", plot = resMV$plots[[2]])

# IVW results
ivwMV <- mv_ivw(mvdat, pval_threshold = 1)

# Format output
OR <- exp(resMV$result[6])
ORup <- exp(resMV$result[6] + 1.96 * resMV$result[7])
ORlp <- exp(resMV$result[6] - 1.96 * resMV$result[7])
pval <- resMV$result[8]

restable <- data.frame(
  outcome = "SSc",
  exposure = c(unique(d1$id.exposure), unique(d2$id.outcome)),
  OR = round(OR, 3),
  `95% OR-CIlp` = round(ORlp, 3),
  `95% OR-CIip` = round(ORup, 3),
  P = round(pval, 3)
)

# Save output tables
raw_name <- paste0("rawresults_", prune_list, "_", tag, ".txt")
summary_name <- paste0("table_", prune_list, "_", tag, ".txt")

write.table(resMV$result, raw_name, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(restable, summary_name, sep = "\t", quote = FALSE, row.names = FALSE)