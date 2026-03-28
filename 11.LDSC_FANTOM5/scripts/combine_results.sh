#!/bin/bash

# Directories
results_dir="/projects/ralab/data/projects/nucCAGE_PRIME_paper/11.LDSC_FANTOM5/ldsc_files/ldsc_results"  # Directory containing .results files
output_file="/projects/ralab/data/projects/nucCAGE_PRIME_paper/11.LDSC_FANTOM5/ldsc_files/combined_results.tsv"  # Output file name

# Initialize the output file with the header
echo -e "Category\tProp._SNPs\tProp._h2\tProp._h2_std_error\tEnrichment\tEnrichment_std_error\tEnrichment_p\tTrait" > "$output_file"

# Process each .results file
for file in "$results_dir"/*.results; do
    # Get the base name (without path and extension)
    trait=$(basename "$file" .results)
    
    # Append the file content to the output file with the new 'trait' column
    # Skip the first line (header) from each file
    tail -n +2 "$file" | awk -v trait="$trait" '{print $0 "\t" trait}' >> "$output_file"
    
    echo "Processed $file"
done

echo "All files have been concatenated into $output_file"