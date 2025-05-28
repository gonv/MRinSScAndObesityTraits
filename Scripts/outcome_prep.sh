#!/bin/bash
set -euo pipefail  

# === Prepare GWAS summary statistics file for outcome (with BETA and SE) ===
# Usage: bash outcome_prep.sh outcome_file.txt
# Output: one per-chromosome file with formatted BETA/SE info

input_file=$1
output_prefix="nw_${input_file}"

# Format the input file: add CHR_BP column and select relevant columns
awk -F"\t" '{print $1 "\t" $1 "_" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' "$input_file" > "$output_prefix"

# Remove the first line (header)
sed -i '1d' "$output_prefix"  

# Split by chromosome
for chr in $(seq 1 22); do
  echo "awk -F\t '\$1 == ${chr} {print \$0}' $output_prefix > chr${chr}_$output_prefix" >> Split_outcome.sh
done

bash Split_outcome.sh
rm Split_outcome.sh

# Create header for output files
echo -e "chr\tCHR_BP\tSNP\tA1\tP\tBETA\tSE\tOR" > Header_outcome.txt

# Concatenate header with data for each chromosome
for chr in $(seq 1 22); do
  cat Header_outcome.txt chr${chr}_$output_prefix > BETA_SE_chr${chr}_$output_prefix_Bis
done

# Prepare SNP extraction list for exposure
awk -F"\t" '{print $1 "_" $2}' "$input_file" > "ext_${input_file}_i"
sed -i '1d' "ext_${input_file}_i"  
echo -e "SNP" > extHeader
cat extHeader "ext_${input_file}_i" > extract_outcome.txt

# Clean up intermediate files
for chr in $(seq 1 22); do rm chr${chr}_$output_prefix; done
rm Header_outcome.txt extHeader "ext_${input_file}_i"
