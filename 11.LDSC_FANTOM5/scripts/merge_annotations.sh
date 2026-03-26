#!/bin/bash

# Directories
annot_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/annotation_files"  # Directory containing annotation files
output_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/merged_annotations"    # Output directory for merged files
bim_dir="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/plink_files"                 # Path to .bim files

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Extract unique score and chromosome pairs from filenames
score_chr_pairs=$(ls "$annot_dir" | grep -oP '_\K[0-9]+\.[0-9]+.[0-9]+' | sort -u)

# Function to merge files for a given score and chromosome
merge_annotations() {
    score_chr=$1
    score=$(echo "$score_chr" | cut -d. -f1,2)  # Extract score (e.g., 0.75)
    chr=$(echo "$score_chr" | cut -d. -f3)  # Extract chromosome (e.g., 2)

    echo "Processing score: ${score}, chromosome: ${chr}..."

    # Define input bim file
    bim_file="${bim_dir}/1000G.EUR.hg38.${chr}.bim"
    if [[ ! -f "$bim_file" ]]; then
        echo "BIM file for chromosome ${chr} not found! Skipping..."
        return
    fi

    # Define output merged file (before compression)
    merged_file="${output_dir}/combined.${score}.${chr}.annot"

    # Initialize merged file with metadata columns: CHR, BP, SNP, CM, base
    awk '{print $1, $4, $2, $3, 1}' "$bim_file" > "$merged_file"

    # Initialize header
    header="CHR\tBP\tSNP\tCM\tbase"

    # Find annotation files matching this score and chromosome
    annot_files=$(ls "${annot_dir}"/*_"${score}.${chr}".annot.gz)

    if [[ -z "$annot_files" ]]; then
        echo "No annotation files found for score ${score}, chromosome ${chr}. Skipping..."
        rm -f "$merged_file"
        return
    fi

    # Loop over annotation files
    for annot_file in $annot_files; do
        # Extract feature name
        feature_name=$(basename "$annot_file" | sed -E "s/([a-zA-Z0-9._-]+_[0-9]+\.[0-9]+)\.[0-9]+\.annot.gz/\1/")


        # Create temporary filenames unique for this chromosome
        temp_col="temp_col_"${score}.${chr}".tsv"
        temp_combined="temp_combined_"${score}.${chr}".tsv"

        # Extract annotation values (skip header) into a temporary column file
        zcat "$annot_file" | tail -n +2 > "$temp_col"

        # Paste the current merged file and the new column together
        paste "$merged_file" "$temp_col" > "$temp_combined"
        mv "$temp_combined" "$merged_file"

        # Append the feature name to the header
        header="${header}\t${feature_name}"

        # Remove temporary column file
        rm -f "$temp_col"
    done

    # Prepend the header to the merged file
    sed -i "1i${header}" "$merged_file"

    # Compress the final merged file
    gzip -f "$merged_file"

    echo "Merged file for score ${score}, chromosome ${chr} created at ${merged_file}.gz"
}

export -f merge_annotations
export annot_dir output_dir bim_dir

# Process all score-chromosome pairs in parallel
parallel -j 25 merge_annotations ::: $score_chr_pairs

