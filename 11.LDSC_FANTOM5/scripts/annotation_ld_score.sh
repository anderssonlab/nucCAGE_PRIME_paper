#!/bin/bash

# Load required modules
module load anaconda3/5.3.1
module load ldsc

source activate ldsc

# Directories
plink_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/plink_files"
annot_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/merged_annotations"
output_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/merged_annotations"
snplist="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/w_hm3.snplist"

# Ensure script stops on errors
set -e

# Extract unique score and chromosome pairs from filenames in the annotation directory
score_chr_pairs=$(ls "$annot_dir" | grep -oP 'combined\.\K[0-9]+\.[0-9]+\.[0-9]+' | sort -u)

# Function to calculate LD scores for a given score and chromosome
calculate_ldsc() {
    score_chr=$1
    plink_dir=$2
    annot_dir=$3
    output_dir=$4
    snplist=$5
    score=$(echo "$score_chr" | cut -d. -f1,2)  # Extract score (e.g., 0.75)
    chr=$(echo "$score_chr" | cut -d. -f3)  # Extract chromosome (e.g., 2)

    output_file="${output_dir}/combined.${score}.${chr}.l2.M_5_50"

    # Skip if output file already exists
    if [[ -f "$output_file" ]]; then
        echo "Skipping score ${score}, chromosome ${chr}: Output already exists."
        return
    fi

    echo "Calculating LD scores for score ${score}, chromosome ${chr}..."

    ldsc.py \
        --l2 \
        --bfile "${plink_dir}/1000G.EUR.hg38.${chr}" \
        --annot "${annot_dir}/combined.${score}.${chr}.annot.gz" \
        --out "${output_dir}/combined.${score}.${chr}" \
        --print-snps "$snplist" \
        --ld-wind-cm 1 || {
        echo "Error processing score ${score}, chromosome ${chr}" >&2
        exit 1
    }

    echo "Finished processing score ${score}, chromosome ${chr}"
}

export -f calculate_ldsc

# Run all score-chromosome pairs in parallel, making sure the full path is used for output
parallel --joblog parallel.log -j 22 calculate_ldsc {} "$plink_dir" "$annot_dir" "$output_dir" "$snplist" ::: $score_chr_pairs

echo "LD score calculation completed!"

