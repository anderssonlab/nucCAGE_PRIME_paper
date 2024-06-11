#!/bin/bash

module load STAR
module load bedops
module load samtools
module load bedtools

## Generate genome reads
#./kmer_fasta.pl fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta 73 100000000 genome_reads/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/

#STAR_CMD="STAR --readFilesType Fastx --alignEndsType Extend5pOfRead1 --readQualityScoreBase 33 --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --genomeDir /projects/ralab/data/genome/hg38/STAR --runThreadN 30 --outSAMtype SAM --outSAMattributes NH HI AS MD nM"

#for FILE in $( ls genome_reads/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/bin.*.fa ); do
#    OUT=$( basename $FILE | sed -e "s/.fa//" )
    
#    CMD="${STAR_CMD} --readFilesIn ${FILE} --outFileNamePrefix mapped/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/${OUT}"	

#    eval " $CMD" >> GRCh38_no_alt_analysis_set_GCA_000001405.15.73.log 2>&1
#done

## Convert to bed
#for FILE in $( ls mapped/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/bin.*.sam ); do
#	OUT=$( echo $FILE | sed -e "s/.sam/.bed/" )
#	convert2bed --input=sam --do-not-sort < $FILE > $OUT
#done

TEMPDIR=$( mktemp -d -p /projects/ralab/scratch/ )

## Extract 5' ends
cat mapped/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/bin.*.bed | awk 'BEGIN{OFS="\t"}{print $1, $2, $2+1, $4, 1, "+"}' > $TEMPDIR/plus &
cat mapped/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/bin.*.bed | awk 'BEGIN{OFS="\t"}{print $1, $3-1, $3, $4, 1, "-"}' > $TEMPDIR/minus &

wait

sort -k 1,1 -k 2,2n -o $TEMPDIR/plus $TEMPDIR/plus --parallel=20 --temporary-directory=$TEMPDIR &
sort -k 1,1 -k 2,2n -o $TEMPDIR/minus $TEMPDIR/minus --parallel=20 --temporary-directory=$TEMPDIR &

wait

bedtools merge -i $TEMPDIR/plus | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, NR, 1, "+"}' > tracks/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/GRCh38_no_alt_analysis_set_GCA_000001405.15.73.uniquely.mappable.plus.bed &
bedtools merge -i $TEMPDIR/minus | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, NR, 1, "-"}' > tracks/GRCh38_no_alt_analysis_set_GCA_000001405.15.73/GRCh38_no_alt_analysis_set_GCA_000001405.15.73.uniquely.mappable.minus.bed &

wait

rm -rf $TEMPDIR
