
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

############################################
# Directories
############################################

SRC_DIR="../0.External_resources/FANTOM5_cellFacet_bedFiles"
FACET_DIR="./PRIME_FANTOM5_facets"
AGNOSTIC_DIR="./PRIME_FANTOM5_agnostic"

mkdir -p "$FACET_DIR" "$AGNOSTIC_DIR"

############################################
# 1. Produce FANTOM facet BED6 files
############################################

echo "Reformating bed files..."

NPROC=$(nproc)

for FILE_PATH in "$SRC_DIR"/*_qn.PL.score0.5.bed; do
(    
    base=$(basename "$FILE_PATH")
    prefix=${base%%_qn.PL.score0.5.bed}

    OUT_FILE="$FACET_DIR/PRIME_FANTOM5_${prefix}_0.5.bed"

    awk 'BEGIN{OFS="\t"}
    {
        score=$4
        name=$1":"$2"-"$3
        print $1,$2,$3,name,score,"*"
    }' "$FILE_PATH" | sort -k 1,1 -k 2,2n > "$OUT_FILE"
) &

while (( $(jobs -r | wc -l) >= NPROC )); do
    sleep 0.2
done
done

wait

# ############################################
# # 2. Create pooled ≥0.75 dataset
# ############################################

INPUT_FILE="../8.Genomewide_prediction/FANTOM5_rmSingletons/PRIMEloci_pred_0_75_FANTOM5_rmSingletons_combined_coreovlwith-d.bed"
OUTPUT_FILE="${AGNOSTIC_DIR}/PRIME_FANTOM5_pooled_0.75.bed"

sed 1d "$INPUT_FILE" | awk 'BEGIN{FS="\t"; OFS="\t"} {
    name = $1":"$2"-"$3;    
    print $1, $2, $3, name, $5, $6, $7, $8
}' | sort --parallel=$(nproc) -k1,1 -k2,2n > "$OUTPUT_FILE"

############################################
# 3. Merge facet GREs
############################################

cd "$FACET_DIR" || exit

merge_facets() {
    local pattern=$1
    local out_name=$2
    local min_score=$3

    echo "Generating $out_name (min Score: $min_score)..."

    awk -v min="$min_score" '
    BEGIN{FS=OFS="\t"}
    {
        if ($5 < min) next;

        split(FILENAME, a, "_");
        facet = a[3];

        name = $4;

        if (!(name in seen)) {
            chr[name]=$1; start[name]=$2; end[name]=$3;
            score[name]=$5; strand[name]=$6;
            facets[name]=facet;
            seen[name]=1;
        } else {
            if ($5 > score[name]) score[name]=$5;

            if (!match(";"facets[name]";", ";"facet";")) {
                facets[name] = facets[name]";"facet;
            }
        }
    }
    END {
        for (n in chr) {
            print chr[n], start[n], end[n], n, score[n], strand[n], facets[n]
        }
    }' $pattern | sort -k1,1 -k2,2n > "${out_name}.bed"
}

merge_facets "*_proximal_0.5.bed" "../${AGNOSTIC_DIR}/PRIME_FANTOM5_agnostic_proximal_0.5" 0.0 &
merge_facets "*_proximal_0.5.bed" "../${AGNOSTIC_DIR}/PRIME_FANTOM5_agnostic_proximal_0.75" 0.75 &

merge_facets "*_distal_0.5.bed" "../${AGNOSTIC_DIR}/PRIME_FANTOM5_agnostic_distal_0.5" 0.0 &
merge_facets "*_distal_0.5.bed" "../${AGNOSTIC_DIR}/PRIME_FANTOM5_agnostic_distal_0.75" 0.75 &

merge_facets "*_0.5.bed" "../${AGNOSTIC_DIR}/PRIME_FANTOM5_agnostic_all_0.5" 0.0 &
merge_facets "*_0.5.bed" "../${AGNOSTIC_DIR}/PRIME_FANTOM5_agnostic_all_0.75" 0.75 &

wait

############################################
# 4. Generate facet statistics
############################################

gawk '
BEGIN {
    FS=OFS="\t"
    print "facet", "proximal_0.5", "proximal_0.75", "distal_0.5", "distal_0.75", "all_0.5", "all_0.75"
}
{
    split(FILENAME, a, "_")
    facet = a[3]; type = a[4]; score = $5 + 0
    facets[facet] = 1
    
    if(score >= 0.5)  all50[facet]++
    if(score >= 0.75) all75[facet]++

    if(type == "proximal") {
        if(score >= 0.5)  prox50[facet]++
        if(score >= 0.75) prox75[facet]++
    } else if(type == "distal") {
        if(score >= 0.5)  dist50[facet]++
        if(score >= 0.75) dist75[facet]++
    }
}
END {
    for(f in all50) {
        p50 = prox50[f] + 0; p75 = prox75[f] + 0
        d50 = dist50[f] + 0; d75 = dist75[f] + 0
        a50 = all50[f] + 0; a75 = all75[f] + 0

        print f, p50, p75, d50, d75, a50, a75

        count++
        list[2][count] = p50; sum[2] += p50
        list[3][count] = p75; sum[3] += p75
        list[4][count] = d50; sum[4] += d50
        list[5][count] = d75; sum[5] += d75
        list[6][count] = a50; sum[6] += a50
        list[7][count] = a75; sum[7] += a75
    }

    if (count > 0) {
        printf "mean"
        for(i=2; i<=7; i++) printf "\t%.2f", sum[i]/count
        print ""

        printf "median"
        for(i=2; i<=7; i++) {
            n = asort(list[i])
            if (n % 2 == 1) {
                med = list[i][(n+1)/2]
            } else {
                med = (list[i][n/2] + list[i][n/2+1]) / 2
            }
            printf "\t%g", med
        }
        print ""
    }
}
' PRIME_FANTOM5_*_0.5.bed > "../PRIME_FANTOM5_agnostic/facet_statistics.tsv"

