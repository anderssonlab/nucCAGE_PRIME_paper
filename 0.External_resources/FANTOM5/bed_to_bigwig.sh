#!/bin/bash

# Define the necessary paths
OUT="/projects/ralab/data/projects/nucCAGE_PRIME_paper/0.External_resources/FANTOM5/primary_bw"
WD="/projects/ralab/data/projects/nucCAGE_PRIME_paper/0.External_resources/FANTOM5/primary_ctss"
TEMPDIR="/projects/ralab/data/projects/nucCAGE_PRIME_paper/0.External_resources/FANTOM5/primary_bw/tmp"
SCRIPTDIR="/projects/ralab/apps/"
GENOME_PATH="/projects/ralab/data/genome/hg38"
GENOME="hg38"

# Input file (gzipped)
INPUT_FILE="$1"  # Pass the gzipped file as an argument

# Extract the PREFIX (everything before ".hg38.nobarcode.ctss.bed.gz")
PREFIX=$(basename "$INPUT_FILE" .hg38.nobarcode.ctss.bed.gz)

# Unzip the input file to a temporary location
zcat ${WD}/${PREFIX}.hg38.nobarcode.ctss.bed.gz > ${TEMPDIR}/${PREFIX}.ctss.bed


# Create the minus and plus bedgraph files
grep -F ",-" ${TEMPDIR}/${PREFIX}.ctss.bed | cut -f 1,2,3,5 | sort -k1,1 -k2,2n > ${TEMPDIR}/${PREFIX}.minus.bedgraph
grep -F ",+" ${TEMPDIR}/${PREFIX}.ctss.bed | cut -f 1,2,3,5 | sort -k1,1 -k2,2n > ${TEMPDIR}/${PREFIX}.plus.bedgraph

# Convert bedgraph files to BigWig format
${SCRIPTDIR}/bin/bedGraphToBigWig ${TEMPDIR}/${PREFIX}.minus.bedgraph ${GENOME_PATH}/chrom_size/${GENOME}.chrom.sizes ${OUT}/${PREFIX}.minus.bw
${SCRIPTDIR}/bin/bedGraphToBigWig ${TEMPDIR}/${PREFIX}.plus.bedgraph ${GENOME_PATH}/chrom_size/${GENOME}.chrom.sizes ${OUT}/${PREFIX}.plus.bw

# Clean up the temporary bed file if necessary
rm ${TEMPDIR}/${PREFIX}.ctss.bed
