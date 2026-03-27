

# consider only pip>=0.1 variants
gzip -dc UKBB_94traits_release1.bed.gz | awk 'BEGIN{FS="\t"} $18>=0.1' > UKBB_94traits_release1_pip01.bed

# clean up formatting errors (notation and missing fields)
sed -e 's/8.1e+07/81000000/g' < UKBB_94traits_release1_pip01.bed | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, 1, "*", $5, $12, $18}' | sort -k 1,1 -k 2,2n > UKBB_94traits_release1_pip01_clean.bed

# make unique
module load bedtools
bedtools groupby -i UKBB_94traits_release1_pip01_clean.bed -g 1,2,3,4,5,6,7,8 -c 9 -o max > UKBB_94traits_release1_pip01_clean_uniq.bed

# lift over to hg38 coordinates
./liftOver -bedPlus=3 UKBB_94traits_release1_pip01_clean_uniq.bed hg19ToHg38.over.chain UKBB_94traits_release1_pip01_clean_uniq_hg38.bed UKBB_94traits_release1_pip01_clean_uniq_unmapped.bed
