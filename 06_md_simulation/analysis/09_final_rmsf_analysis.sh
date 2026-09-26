#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"
RMSD="$ROOT/06_md_simulation/analysis/results/RMSD"
OUT="$ROOT/06_md_simulation/analysis/results/RMSF"

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

SUMMARY="$OUT/Corrected_MD11_Catalytic_RMSF_Summary.tsv"

printf "System\tTYR28_nm\tTYR48_nm\tHIS83_nm\tHIS109_nm\tSER128_nm\tASN130_nm\tCatalytic_Mean_nm\n" > "$SUMMARY"

for x in "${SYSTEMS[@]}"
do
    echo "========== $x =========="

    TPR="$BASE/$x/Production/md_200ns.tpr"
    XTC="$RMSD/${x}_nojump.xtc"
    XVG="$OUT/${x}_CA_RMSF.xvg"

    echo "3" | gmx rmsf \
        -s "$TPR" \
        -f "$XTC" \
        -o "$XVG" \
        -res \
        >/dev/null 2>&1

    vals=$(awk '
        !/^[@#]/ {
            if ($1==28)  a=$2
            if ($1==48)  b=$2
            if ($1==83)  c=$2
            if ($1==109) d=$2
            if ($1==128) e=$2
            if ($1==130) f=$2
        }
        END {
            mean=(a+b+c+d+e+f)/6
            printf "%.4f %.4f %.4f %.4f %.4f %.4f %.4f",
                   a,b,c,d,e,f,mean
        }
    ' "$XVG")

    read a b c d e f mean <<< "$vals"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$x" "$a" "$b" "$c" "$d" "$e" "$f" "$mean" \
        >> "$SUMMARY"

    echo "Catalytic mean RMSF = $mean nm"
done

echo
echo "FINAL CATALYTIC RMSF SUMMARY:"
column -t "$SUMMARY"


