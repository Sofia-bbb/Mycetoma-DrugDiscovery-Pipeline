#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"
RMSD="$ROOT/06_md_simulation/analysis/results/RMSD"
OUT="$ROOT/06_md_simulation/analysis/results/Rg"

mkdir -p "$OUT"

SYSTEMS=(
NQLMM002110
NQLMM007485
NQLMM002852
NQLMM000788
NQLMM003435
NQLMM002327
NQLMM003691
NQLMM005370
NQLMM000191
Carpropamid
Fenoxanil
)

SUMMARY="$OUT/Corrected_MD11_Rg_Summary.tsv"

printf "System\tRg_Mean_nm\tRg_SD_nm\tRg_Min_nm\tRg_Max_nm\n" > "$SUMMARY"

for x in "${SYSTEMS[@]}"
do
    echo "========== $x =========="

    TPR="$BASE/$x/Production/md_200ns.tpr"
    XTC="$RMSD/${x}_nojump.xtc"
    XVG="$OUT/${x}_Protein_Rg.xvg"

    echo "1" | gmx gyrate \
        -s "$TPR" \
        -f "$XTC" \
        -o "$XVG" \
        >/dev/null 2>&1

    vals=$(awk '
        !/^[@#]/{
            s+=$2
            ss+=$2*$2
            n++

            if(n==1 || $2<min) min=$2
            if(n==1 || $2>max) max=$2
        }
        END{
            m=s/n
            sd=sqrt(ss/n-m*m)

            printf "%.4f %.4f %.4f %.4f",
                   m,sd,min,max
        }
    ' "$XVG")

    read mean sd min max <<< "$vals"

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$x" "$mean" "$sd" "$min" "$max" \
        >> "$SUMMARY"

    echo "Rg = $mean ± $sd nm"
done

echo
echo "FINAL Rg SUMMARY:"
column -t "$SUMMARY"


