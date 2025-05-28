#!/bin/bash
set -euo pipefail

# === Format outcome GWAS file for MVMR ===
# Input: GWAS file with >= 10 columns
# Output: reformatted file with standard MR columns and SNP ID

input="$1"
intermediate="${input}_i"

# Add SNP ID column (CHR_POS)
awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $1 "_" $2 "\t" $10}' "$input" > "$intermediate"

# Remove original file header
sed -i '1d' "$intermediate"

# Add standardized header and merge
echo -e "chr\tposition\teffect_allele\tother_allele\tP\tbeta\tse\tZ\tsamplesize\tSNP\trsID" > Header_META
cat Header_META "$intermediate" > "$input"

# Cleanup
rm Header_META "$intermediate"
