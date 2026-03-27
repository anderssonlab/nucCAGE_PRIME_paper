ls /projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/NET_CAGE/fastq/*.fastq | parallel -P 6 '/projects/ralab/apps/scripts/CAGE_STAR_mapping.sh -g hg38 -b 0 -i true -f {}'
