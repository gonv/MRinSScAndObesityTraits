#!/bin/bash
set -euo pipefail

# === Multivariable Mendelian Randomization Skeleton Script ===
# Input options:
#   -o  Outcome GWAS file
#   -p  List of PheWAS files (traits to use as exposures)
#   -i  Index SNPs file (optional)
#   -k  Secondary SNPs file (optional)
#   -d  Output directory name

# === Argument parsing ===
usage="$(basename "$0") [-h] -o outcome -p phewas_list -i indexSNP -k secondarySNP -d directory"
while getopts ':h:o:p:i:k:d:' option; do
  case "$option" in
    h) echo "$usage"; exit ;;
    o) outcome=$OPTARG ;;
    p) phewaslist=$OPTARG ;;
    i) indexSNP=$OPTARG ;;
    k) secundarySNP=$OPTARG ;;
    d) directory=$OPTARG ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# === Prepare directory and copy necessary files ===
mkdir "temp_${directory}"
cd "temp_${directory}"
cp ../*.R ../*.sh ../*.py ../"${outcome}" ../"${phewaslist}" .
for fi in $(<"${phewaslist}"); do cp ../"${fi}" .; done

# === Format outcome file and clump all PheWAS ===
bash outcome_prep.sh "${outcome}"
for fi in $(<"${phewaslist}"); do
  bash PheWAS_clump.sh "${fi}" "${outcome}" "${indexSNP}" "${secundarySNP}"
done

# === Merge all SNPs from clumped files ===
cat Clumping_* > UniqueSNPs_i
uniq UniqueSNPs_i > UniqueSNPs

# === Perform pruning with PLINK ===
plink1.9 --bfile /path/to/ref_panel \
  --extract UniqueSNPs \
  --indep-pairwise 50 5 0.05 \
  --out UniqueSNPs_CombineSNPs

# === Annotate PheWAS traits with IDs and standard format ===
echo -e "chr\tposition\teffect_allele\tother_allele\tbeta\tse\tpval\tsamplesize\tid\tSNP" > Header_MR
for fi in $(<"${phewaslist}"); do
  grep -Fwf UniqueSNPs_CombineSNPs.prune.in "${fi}_MRMV_preclump" > "${fi}_pruneSNP_i"
  awk -F',' -v trait="${fi}" '{OFS="\t"; $(NF+1)=trait; print}' "${fi}_pruneSNP_i" > "${fi}_pruneSNP_ii"
  awk -F"_" '{print $1 "\t" $2}' "${fi}_pruneSNP_ii" > "${fi}_pruneSNP_iii"
  awk -F"\t" '{print $0 "\t" $1 "_" $2 }' "${fi}_pruneSNP_iii" > "${fi}_pruneSNP_vi"
  cat Header_MR "${fi}_pruneSNP_vi" > "${fi}_pruneSNP"
  rm "${fi}_pruneSNP_"*
done

rm Header_MR
ls *_pruneSNP > prune.list

# === Run multivariable MR analysis in R ===
Rscript MRMV_fun_v2.R "${outcome}" prune.list "${directory}"

echo "Completed MVMR analysis for directory: ${directory}"