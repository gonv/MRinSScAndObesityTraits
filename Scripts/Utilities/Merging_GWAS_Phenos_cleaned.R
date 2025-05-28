#!/usr/bin/env Rscript

# === Merge outcome and exposure GWAS files by shared SNPs (for a given chromosome) ===
# Input:
#   1. Outcome file with CHR_BP and effect columns
#   2. Exposure file with chr_pos and effect columns
#   3. Chromosome number to process
# Output:
#   Tab-separated file with merged columns for shared SNPs, no header

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
outcome_file <- args[1]
exposure_file <- args[2]
chr <- args[3]

# Load files
gwas <- read.table(outcome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pheno <- read.table(exposure_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Separate chromosomal positions
gwas <- gwas %>% separate("CHR_BP", c("CHR", "BP"), sep = "_", remove = FALSE)
pheno <- pheno %>% separate("chr_pos", c("chr", "pos"), sep = "_", remove = FALSE)

cat(sprintf("Started processing chromosome %s at %s\n", chr, format(Sys.time(), "%c")))

# Filter by chromosome
gwas_chr <- gwas %>% filter(CHR == chr)
pheno_chr <- pheno %>% filter(chr == chr)

# Merge on shared SNPs
df <- merge(gwas_chr, pheno_chr, by.x = "CHR_BP", by.y = "chr_pos")

# Select and reorder columns
merged_df <- df %>% select(CHR, CHR_BP, A1, P, BETA, SE, SNP, alt, beta_EUR, se_EUR, pval_EUR)

# Output file
output_name <- paste("Shared_Beta_SE", outcome_file, exposure_file, sep = "_")
write.table(merged_df, file = output_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("Finished processing chromosome %s at %s\n", chr, format(Sys.time(), "%c")))