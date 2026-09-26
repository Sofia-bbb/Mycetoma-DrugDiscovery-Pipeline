#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"
OUT="$ROOT/06_md_simulation/analysis/results/RMSD"

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

SUMMARY="$OUT/Corrected_MD11_RMSD_Summary.tsv"

printf "System\tProtein_Mean_nm\tProtein_SD_nm\tLigand_Mean_nm\tLigand_SD_nm\n" > "$SUMMARY"

for x in "${SYSTEMS[@]}"
do
    echo "=================================================="
    echo "$x"
    echo "=================================================="

    PROD="$BASE/$x/Production"
    TPR="$PROD/md_200ns.tpr"
    XTC="$PROD/md_200ns.xtc"

    NDX="$OUT/${x}.ndx"
    NOJUMP="$OUT/${x}_nojump.xtc"

    PROT="$OUT/${x}_Protein_Backbone_RMSD.xvg"
    LIG="$OUT/${x}_Ligand_RMSD_nojump.xvg"

    # Create Protein_MOL = group 18
    printf "1 | 13\nname 18 Protein_MOL\nq\n" | \
    gmx make_ndx \
        -f "$TPR" \
        -o "$NDX" \
        >/dev/null 2>&1

    # PBC no-jump trajectory containing protein + ligand
    echo "18" | \
    gmx trjconv \
        -s "$TPR" \
        -f "$XTC" \
        -n "$NDX" \
        -o "$NOJUMP" \
        -pbc nojump \
        >/dev/null 2>&1

    # Protein backbone RMSD
    printf "4\n4\n" | \
    gmx rms \
        -s "$TPR" \
        -f "$NOJUMP" \
        -n "$NDX" \
        -o "$PROT" \
        -tu ns \
        >/dev/null 2>&1

    # Ligand RMSD after backbone fitting
    printf "4\n13\n" | \
    gmx rms \
        -s "$TPR" \
        -f "$NOJUMP" \
        -n "$NDX" \
        -o "$LIG" \
        -tu ns \
        >/dev/null 2>&1

    read PM PSD <<< $(awk '
        !/^[@#]/{
            s+=$2;
            ss+=$2*$2;
            n++
        }
        END{
            m=s/n;
            sd=sqrt(ss/n-m*m);
            printf "%.4f %.4f",m,sd
        }
    ' "$PROT")

    read LM LSD <<< $(awk '
        !/^[@#]/{
            s+=$2;
            ss+=$2*$2;
            n++
        }
        END{
            m=s/n;
            sd=sqrt(ss/n-m*m);
            printf "%.4f %.4f",m,sd
        }
    ' "$LIG")

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$x" "$PM" "$PSD" "$LM" "$LSD" \
        >> "$SUMMARY"

    echo "Protein RMSD: $PM ± $PSD nm"
    echo "Ligand RMSD : $LM ± $LSD nm"
done

echo
echo "FINAL SUMMARY:"
column -t "$SUMMARY"



