#!/bin/bash

module load anaconda3/5.3.1
module load ldsc

source activate ldsc

# Create directory for annotation files if it doesn't exist
mkdir -p /projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/annotation_files

# Define the function that will be run for each BED file
process_bed_file() {
    bed_file=$1
    i=$2
    bim_file="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/plink_files/1000G.EUR.hg38.${i}.bim"
    
    # Extract the basename without the extension
    bed_basename=$(basename "$bed_file" ".bed")
    
    # Define output annotation file path
    annot_file="/projects/ralab/data/projects/nucleiCAGEproject/11.LDSC_FANTOM5/ldsc_files/annotation_files/${bed_basename}.${i}.annot.gz"
    
    # Run the command with the original BED file (no header or prefix to trim anymore)
    make_annot.py \
        --bed-file "$bed_file" \
        --bimfile "$bim_file" \
        --annot-file "$annot_file"
}

export -f process_bed_file  # Export the function to be used by parallel

# Get a list of all BED files
bed_files=(/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/FANTOM5_cellFacet_bedFiles/score_range/*.bed)

# Use parallel to process all combinations of bed files and chromosomes
parallel -j 22 process_bed_file {1} {2} ::: "${bed_files[@]}" ::: {1..22}


