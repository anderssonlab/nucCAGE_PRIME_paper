#!/bin/bash

module load anaconda3/5.3.1
module load ldsc

source activate ldsc

# Directories
stat_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/sumstats"
output_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/ldsc_results"
stat_list="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/sumstats.txt"

# Define the score values (update this list if needed)
scores=("0.5" "0.55" "0.6" "0.65" "0.7" "0.75")  

# Function to process each sumstat file
process_sumstat() {
    file_name=$1
    stat_dir=$2
    output_dir=$3
    scores=("${@:4}")  # Capture all scores as an array

    # Paths
    sumstat_file="${stat_dir}/${file_name}"
    base_name=$(basename "$file_name" | cut -d. -f1)  # Extract everything before the first dot

    # Loop over each score
    for score in "${scores[@]}"; do
        ref_ld_chr="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/merged_annotations/combined.${score}."
        output_check="${output_dir}/${base_name}_score${score}.results"  # Check existence
        output_file="${output_dir}/${base_name}_score${score}"  # Output file

        # Skip if output already exists
        if [[ -f "$output_check" ]]; then
            echo "Skipping '$file_name' for score $score - Output already exists: '$output_check'"
            continue
        fi

        echo "Processing '$sumstat_file' -> '$output_file' for score $score"

        # Run LDSC
        ldsc.py \
            --h2 "$sumstat_file" \
            --ref-ld-chr "$ref_ld_chr" \
            --w-ld-chr /projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/weights/weights.hm3_noMHC. \
            --overlap-annot \
            --frqfile-chr /projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/plink_files/1000G.EUR.hg38. \
            --out "$output_file" || {
            echo "Error processing $file_name for score $score" >&2
            exit 1
        }

        echo "Finished processing $file_name for score $score"
    done
}

export -f process_sumstat

# Run in parallel over all sumstats, passing the scores array
parallel --joblog parallel.log -j 6 process_sumstat {} "$stat_dir" "$output_dir" "${scores[@]}" :::: "$stat_list"
