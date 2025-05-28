
#!/usr/bin/env Rscript

# === Format GWAS summary statistics for MR analysis ===
# This script computes beta and standard error (SE) from GWAS odds ratios (ORs).
# Input 1: summary stats file with columns CHR_BP, A1, P, OR, SNP
# Input 2: chromosome number to process

library(optparse)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]           # GWAS summary statistics file
chromosome <- args[2]           # Chromosome to process

cat(sprintf("Started processing chromosome %s at %s\n", chromosome, format(Sys.time(), "%c")))  

# === Load and preprocess GWAS data ===
gwas <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Separate CHR_BP column into CHR and BP
gwas <- gwas %>% separate("CHR_BP", c("CHR", "BP"), sep = "_", remove = FALSE)

# Filter by chromosome
gwas_chr <- gwas %>% filter(CHR == chromosome)

# Compute SE and BETA
gwas_chr <- gwas_chr %>% 
  mutate(
    SE = round(abs(log(OR) / qnorm(P / 2)), 4),
    BETA = round(log(OR), 4)
  )

# Select and order columns
gwas_chr <- gwas_chr %>% select(CHR_BP, A1, P, BETA, SE, SNP)

# Remove rows with missing SE values
gwas_chr <- gwas_chr[!is.na(gwas_chr$SE), ]

# Write output file
output_file <- sprintf("BETA_SE_chr%s_%s_Bis", chromosome, input_file)  
write.table(gwas_chr, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Finished processing chromosome %s at %s\n", chromosome, format(Sys.time(), "%c")))  
