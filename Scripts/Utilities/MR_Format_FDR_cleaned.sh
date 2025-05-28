#!/bin/bash
# === Format MR output files and generate merged FDR table for GIANT exposures ===

input_list=$1
strat=$2

# Clean and concatenate input files
for fi in $(<"$input_list"); do
    sed '1d' "$fi" > "${fi}_i_FDR"
done
cat *_i_FDR > body
echo -e "Exposure\tMethod\tnSNP\tb\tse\tP" > header.txt
awk -F"\t" '{print $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' body > body_i
cat header.txt body_i > "Giant_MR_${strat}.txt"

# Format names depending on stratification
case "$strat" in
    Female) sed -i -e 's/ /_/g' -e 's/.giant-ukbb.meta-analysis.females.23May2018.txt//g' "Giant_MR_${strat}.txt" ;;
    Male) sed -i -e 's/ /_/g' -e 's/.giant-ukbb.meta-analysis.males.23May2018.txt//g' "Giant_MR_${strat}.txt" ;;
    *) sed -i -e 's/ /_/g' -e 's/.giant-ukbb.meta-analysis23May2018.txt//g' "Giant_MR_${strat}.txt" ;;
esac

# Cleanup
rm body* *_i_* header.txt
echo "Done."