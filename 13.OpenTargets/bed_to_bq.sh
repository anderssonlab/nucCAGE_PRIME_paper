
#!/bin/bash

# Usage: ./bed_to_bq.sh input.bed output.csv
INPUT_BED=$1
OUTPUT_CSV=$2

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bed output.csv"
    exit 1
fi

echo "chromosome,start,end,name,score,strand,facets" > "$OUTPUT_CSV"

awk -F'\t' -v OFS=',' '
{
    gsub(/^chr/, "", $1); 
    $2 = $2 + 1; 
    print $1, $2, $3, $4, $5, $6, "\""$7"\""
}' "$INPUT_BED" >> "$OUTPUT_CSV"
