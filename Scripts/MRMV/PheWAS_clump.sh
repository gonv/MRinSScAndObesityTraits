#!/bin/bash
set -euo pipefail

# === Prepare and clump PheWAS summary statistics for MVMR ===
# Input: 1 - PheWAS file, 2 - Outcome GWAS, 3 - indexSNP file (optional), 4 - secondarySNP file (optional)

phewas_file="$1"
outcome_file="$2"
indexSNP="${3:-none}"
secondarySNP="${4:-none}"

# === NEALE formatting ===
if [[ "$phewas_file" == *"irnt"* && "$phewas_file" == *"gwas"* ]]; then
    echo "Formatting NEALE PheWAS"
    python PheWAS_prep.py "$phewas_file"

    echo -e "SNP\tA1\tA2\tminor_allele\tn_complete_samples\tbeta\tse\tP" > header_neale_preclump
    awk -F":" '{print $1 " " $2 " " $3 " " $4 " " $5}' "${phewas_file}_MRMV" > "${phewas_file}_MR_i"
    awk -F" " '{print $1 "_" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' "${phewas_file}_MR_i" > "${phewas_file}_MR_ii"
    cat header_neale_preclump "${phewas_file}_MR_ii" > "${phewas_file}_MRMV_preclump"

    awk -F"\t" '{print $1}' "${phewas_file}_MR_ii" > exposure_extract_i
    sed -i '1d' exposure_extract_i
    echo "SNP" > Header_exp.txt
    cat Header_exp.txt exposure_extract_i > exposure_extract

    plink1.9 --bfile /path/to/ref_panel \
        --clump "${phewas_file}_MRMV_preclump" \
        --clump-field P --clump-kb 5000 \
        --clump-p1 5e-9 --clump-p2 5e-8 --clump-r2 0.05 \
        --out "index_temp_${phewas_file}"

    awk '{print $3}' "index_temp_${phewas_file}.clumped" > index_clump_SNPs
    grep -Fwf index_clump_SNPs "${phewas_file}_MR_ii" > "nw_${phewas_file}_index_i"
    sort -V -k1,1 -k2,2 "nw_${phewas_file}_index_i" > "nw_${phewas_file}_index_ii"
    cat header_neale_preclump "nw_${phewas_file}_index_ii" > "nw_${phewas_file}_index"

    plink1.9 --bfile /path/to/ref_panel \
        --clump "nw_${phewas_file}_index,${outcome_file}" \
        --extract exposure_extract \
        --clump-best --clump-replicate --clump-index-first \
        --clump-allow-overlap --clump-p1 1e-4 --clump-p2 1 \
        --clump-kb 250 --clump-r2 0.2 \
        --out "nw_clump_${phewas_file}"

    awk '{print $1}' "nw_clump_${phewas_file}.clumped.best" > index_clumpedSNPs
    sed -i '1d' index_clumpedSNPs
    mv index_clumpedSNPs "Clumping_${phewas_file}_SScextract.clumped.SNP"
    rm -f Header* header* exposure_extract* "${phewas_file}_MR_i"*
fi

# === GIANT formatting ===
if [[ "$phewas_file" == *"giant-ukbb.meta-analysis"* ]]; then
    echo "Formatting GIANT PheWAS"
    awk -F" " '{print $1 "_" $2 "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $9 "\t" $3}' "$phewas_file" > "nw_${phewas_file}_i"
    awk -F":" '{print $1}' "nw_${phewas_file}_i" > "nw_${phewas_file}_ii"
    sed -i '1d' "nw_${phewas_file}_ii"
    echo -e "SNP\tA1\tA2\tBETA\tSE\tP\trsID" > Header_G.txt
    cat Header_G.txt "nw_${phewas_file}_ii" > "${phewas_file}_MRMV_preclump"

    awk -F" " '{print $1 "_" $2}' "$phewas_file" > exposure_extract_i
    sed -i '1d' exposure_extract_i
    echo "SNP" > Header_exp.txt
    cat Header_exp.txt exposure_extract_i > exposure_extract

    if [[ "$indexSNP" != "none" ]]; then
        cp ../"${phewas_file}_indexSNP" .
        awk -F" " '{print $2 "_" $3 "\t" $9}' "${phewas_file}_indexSNP" > "${phewas_file}_index_i"
        sed -i '1d' "${phewas_file}_index_i"
        sort -V -k1,1 -k2,2 "${phewas_file}_index_i" > "${phewas_file}_index_ii"
        echo -e "SNP\tP" > Header_index.txt
        cat Header_index.txt "${phewas_file}_index_ii" > "${phewas_file}_index"

        plink1.9 --bfile /path/to/ref_panel \
            --clump "${phewas_file}_index,${outcome_file}" \
            --clump-best --clump-replicate --clump-index-first \
            --clump-allow-overlap --clump-p1 1e-4 --clump-p2 1 \
            --clump-kb 250 --clump-r2 0.2 \
            --out "nw_clump_index_${phewas_file}"

        awk '{print $1}' "nw_clump_index_${phewas_file}.clumped.best" > index_clumpedSNPs
        sed -i '1d' index_clumpedSNPs

        if [[ "$secondarySNP" != "none" ]]; then
            cp ../"${phewas_file}_secundarySNP" .
            awk -F" " '{print $2 "_" $3 "\t" $9}' "${phewas_file}_secundarySNP" > "${phewas_file}_sec_i"
            sed -i '1d' "${phewas_file}_sec_i"
            sort -V -k1,1 -k2,2 "${phewas_file}_sec_i" > "${phewas_file}_sec_ii"
            echo -e "SNP\tP" > Header_sec.txt
            cat Header_sec.txt "${phewas_file}_sec_ii" > "${phewas_file}_sec"

            plink1.9 --bfile /path/to/ref_panel \
                --clump "${phewas_file}_sec,${outcome_file}" \
                --clump-best --clump-replicate --clump-index-first \
                --clump-allow-overlap --clump-p1 1e-4 --clump-p2 1 \
                --clump-kb 250 --clump-r2 0.2 \
                --out "nw_clump_secundary_${phewas_file}"

            awk '{print $1}' "nw_clump_secundary_${phewas_file}.clumped.best" > secundary_clumpedSNPs
            sed -i '1d' secundary_clumpedSNPs
            cat index_clumpedSNPs secundary_clumpedSNPs > "Clumping_${phewas_file}_SScextract.clumped.SNP"
        else
            cp index_clumpedSNPs "Clumping_${phewas_file}_SScextract.clumped.SNP"
        fi
    fi

    echo "Clumping complete for ${phewas_file}"
    rm -f Header* header* exposure_extract* "${phewas_file}_i"*
fi