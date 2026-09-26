#!/bin/bash
#SBATCH --job-name=HBOND_MD11
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=06_md_simulation/analysis/results/HBonds/HBOND_%j.out
#SBATCH --error=06_md_simulation/analysis/results/HBonds/HBOND_%j.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gromacs/2024.6

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"
RMSD="$ROOT/06_md_simulation/analysis/results/RMSD"
OUT="$ROOT/06_md_simulation/analysis/results/HBonds"

cd "$ROOT"
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

SUMMARY="$OUT/Corrected_MD11_HBond_Summary.tsv"

printf "System\tMean_HBonds\tSD_HBonds\tMax_HBonds\tFrames_With_HBond\tOccupancy_Percent\n" > "$SUMMARY"

for x in "${SYSTEMS[@]}"
do
    echo "========== $x =========="

    TPR="$BASE/$x/Production/md_200ns.tpr"
    XTC="$RMSD/${x}_nojump.xtc"
    XVG="$OUT/${x}_Protein_Ligand_HBonds.xvg"

    # Protein <-> ligand hydrogen bonds
    gmx hbond \
        -s "$TPR" \
        -f "$XTC" \
        -r "group Protein" \
        -t "group MOL" \
        -num "$XVG" \
        >/dev/null 2>&1

    vals=$(awk '
        !/^[@#]/{
            n++
            s+=$2
            ss+=$2*$2

            if(n==1 || $2>max) max=$2
            if($2>0) occupied++
        }
        END{
            mean=s/n
            sd=sqrt(ss/n-mean*mean)
            occ=(occupied/n)*100

            printf "%.4f %.4f %.0f %d %.2f",
                   mean,sd,max,occupied,occ
        }
    ' "$XVG")

    read mean sd max occupied occupancy <<< "$vals"

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$x" "$mean" "$sd" "$max" "$occupied" "$occupancy" \
        >> "$SUMMARY"

    echo "Mean H-bonds = $mean"
    echo "Frames with >=1 H-bond = $occupancy %"
done

echo
echo "FINAL H-BOND SUMMARY:"
column -t "$SUMMARY"


