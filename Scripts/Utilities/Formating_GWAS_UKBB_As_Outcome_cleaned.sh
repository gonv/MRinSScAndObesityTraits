#!/bin/bash
# === Format UKBB GWAS file to be used as outcome in MR analysis ===

input=$1
output="${1}_MR"

# Remove header and format fields
sed -i '1d' "$input"
awk -F"\t" '{print $1 "\t" $2 "\t" $1 "_" $2 "_" $3 "_"$4 "\t" $5 "\t" $5 "\t" $15 "\t" $12 "\t" $13 "\t" "NA"}' "$input" > "${input}_Temp"

# Add header and combine
echo -e "chr\tCHR_BP\tSNP\tA1\tA1\tP\tBETA\tSE\tOR" > Header_outcome.txt
cat Header_outcome.txt "${input}_Temp" > "$output"

# Clean up
rm Header_outcome.txt "${input}_Temp"