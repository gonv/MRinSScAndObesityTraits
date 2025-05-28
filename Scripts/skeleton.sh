#!/bin/bash
set -euo pipefail  # Exit on error, unset vars, and pipefail for safety

# === Mendelian Randomization Pipeline ===
# This script runs a two-sample MR analysis from one exposure and multiple outcomes.
# It prepares input GWAS data, finds shared SNPs, and executes MR analysis per exposure.

# === Authors ===
# Gonzalo Villanueva-Martin & Lara Bossini-Castillo
# IPBLN-CSIC & University of Granada

# === Example execution ===
# bash skeleton.sh -o META_GWAS_SSc.tsv -p Obesitytraits.txt -s Overall -e OR

# === Required R packages ===
# optparse, dplyr, ggplot2, devtools, TwoSampleMR, na.tools, knitr

# === Parse arguments ===
usage="$(basename "$0") [-h] -o outcome -p phewas_list -s stratification -e effect [-i indexSNP] [-k secundarySNP]"
while getopts ':h:o:p:s:e:i:k:' option; do
  case "$option" in
    h) echo "$usage"; exit ;;
    o) outcome=$OPTARG ;;
    p) phewaslist=$OPTARG ;;
    s) strat=$OPTARG ;;
    e) effect=$OPTARG ;;
    i) indexSNP=$OPTARG ;;
    k) secundarySNP=$OPTARG ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# === Set up working directory ===
mkdir "Temp_${phewaslist}"
cd "Temp_${phewaslist}"

# === Copy scripts and input data ===
cp ../*.R ../*.sh ../${phewaslist} ../${outcome} .
for fi in $(<${phewaslist}); do cp ../${fi} .; done
echo "All required files copied."

# === Format outcome file depending on effect type ===
if [[ "${effect}" == 'OR' ]]; then
  awk -F"	" '{print $1 "_" $2 "	" $4 "	" $7 "	" $8 "	" $9 "	" $3}' ${outcome} > nw_${outcome}
  for chr in $(seq 1 22); do Rscript MR_Formating.R nw_${outcome} ${chr}; done
  mkdir ../${outcome}_BETAS
  cp BETA_SE_chr*_nw_${outcome}_Bis ../${outcome}_BETAS
elif [[ "${effect}" == 'BETA' ]]; then
  bash outcome_prep.sh ${outcome}
  mkdir ../${outcome}_BETAS
  cp BETA_SE_chr*_nw_${outcome}_Bis ../${outcome}_BETAS
else
  echo "Effect must be 'OR' or 'BETA'. Exiting." ; exit 1
fi

# === Loop through each exposure GWAS file ===
for fi in $(<${phewaslist}); do
  echo "Processing ${fi}"
  mkdir temp_${fi}
  [[ ${indexSNP} != 'none' ]] && cp ../${fi}_indexSNP .
  [[ ${secundarySNP} != 'none' ]] && cp ../${fi}_secundarySNP .
  cp ../${fi} .
  cd temp_${fi}
  cp ../${fi} .

  # Run clumping
  bash ../MR_Clumping_Phenos.sh ${fi} ${outcome} ${indexSNP:-none} ${secundarySNP:-none}

  # Prepare headers
  echo -e "chr_pos	alt	beta_EUR	se_EUR	pval_EUR	rsID" > Header.txt

  # Split phenotype file by chromosome
  if [[ "nw_${fi}" == *"giant-ukbb.meta-analysis"* ]]; then
    for chr in $(seq 1 22); do echo "awk -F"_" '\$1 == ${chr} {print \$0}' nw_${fi} > chr${chr}_nw_${fi}" >> Split_Pheno.sh; done
  elif [[ "nw_${fi}" == *".gwas.imputed_v3."* ]]; then
    for chr in $(seq 1 22); do echo "awk -F"	" '\$1 == ${chr} {print \$2 "	" \$3 "	" \$4 "	" \$5 "	" \$6 }' nw_${fi} > chr${chr}_nw_${fi}" >> Split_Pheno.sh; done
  fi
  bash Split_Pheno.sh && rm Split_Pheno.sh

  for chr in $(seq 1 22); do cat Header.txt chr${chr}_nw_${fi} > chr${chr}_nw_${fi}_Bis; done
  rm Header.txt chr*_nw_${fi}

  # Merge exposure and outcome for shared SNPs
  for chr in $(seq 1 22); do
    cp ../BETA_SE_chr${chr}_nw_${outcome}_Bis .
    Rscript ../Merging_GWAS_Phenos.R BETA_SE_chr${chr}_nw_${outcome}_Bis chr${chr}_nw_${fi}_Bis ${chr}
  done

  # Prepare final merged table
  echo -e "CHR_BP	effect_allele.outcome	pval.outcome	beta.outcome	se.outcome	SNP	effect_allele.exposure	beta.exposure	se.exposure	pval.exposure" > Header2.txt
  cat Shared_Beta_SE_BETA_SE_chr* > Body.txt
  sort -n -k1 Body.txt | cut -f2- > Body_sorted.txt
  cat Header2.txt Body_sorted.txt > Shared_Beta_SE_nw_${outcome}_nw_${fi}
  uniq Shared_Beta_SE_nw_${outcome}_nw_${fi} | awk -F"	" '$3 != 0' > Shared_Beta_outcome_${outcome}_exposure_${fi}

  # Run MR
  Rscript ../MR_2SMR_fun.R Shared_Beta_outcome_${outcome}_exposure_${fi} ${fi} > log_MR_fun_${fi}.txt

  cd ..
done

# === Organize final results ===
mkdir ../${phewaslist}_results
cd ..
cp Temp_${phewaslist}/*/Shared_Beta_outcome_* ${phewaslist}_results || true
cp Temp_${phewaslist}/*/mr_* ${phewaslist}_results || true
cp Temp_${phewaslist}/*/presso* ${phewaslist}_results || true
cp Temp_${phewaslist}/*/table_MR_* ${phewaslist}_results || true
cp Temp_${phewaslist}/*/log_MR_fun_* ${phewaslist}_results || true
cp Temp_${phewaslist}/*/*.png ${phewaslist}_results || true
cp Temp_${phewaslist}/*/*.clumped ${phewaslist}_results || true
cp Temp_${phewaslist}/*/*.rds ${phewaslist}_results || true

tar -czf Temp_${phewaslist}.tar.gz Temp_${phewaslist} && cp Temp_${phewaslist}.tar.gz ${phewaslist}_results

# === End of script ===
