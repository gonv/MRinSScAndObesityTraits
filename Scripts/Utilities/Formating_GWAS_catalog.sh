#!/bin/bash
# === Format GWAS Catalog file for MR analysis ===
# Input: GWAS summary statistics with at least 9 columns
# Output: Reformatted file with columns ready for MR input

awk -F"\t" '{print $1 "\t" $3 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $11 "\t" $8 "\t" $9}' "$1" > "${1}_MR"
