ls ./primary_ctss/*.hg38.nobarcode.ctss.bed.gz | parallel -P 20 "./bed_to_bigwig.sh {}" 
