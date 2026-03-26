
#!/usr/bin/env bash

set -euo pipefail

INPUT="20260302_variant_summary.txt.gz"

OUT_VUS="ClinVar_GRCh38_VUS_SNV.txt"
OUT_LP="ClinVar_GRCh38_LikelyPathogenic_SNV.txt"
OUT_P="ClinVar_GRCh38_Pathogenic_SNV.txt"

echo "Filtering ClinVar variants..."

zcat "$INPUT" | awk '
BEGIN{
    FS=OFS="\t"
}

NR==1 {
    print > "'$OUT_VUS'"
    print > "'$OUT_LP'"
    print > "'$OUT_P'"
    next
}

$2=="single nucleotide variant" &&
$17=="GRCh38" &&
$7 !~ /Conflicting/ {

    if ($7 ~ /Uncertain significance/)
        print >> "'$OUT_VUS'"

    else if ($7 ~ /Likely pathogenic/)
        print >> "'$OUT_LP'"

    else if ($7 ~ /Pathogenic/)
        print >> "'$OUT_P'"
}
'

echo "Processing ClinVar variants in R..."

Rscript filter_clinvar.R
