#!/bin/bash
#SBATCH --job-name=RMSD_NP_MD11
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-11
#SBATCH --output=06_md_simulation/analysis/results/RMSD_NoPBC/RMSD_%A_%a.out
#SBATCH --error=06_md_simulation/analysis/results/RMSD_NoPBC/RMSD_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load gromacs/2024.6

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"
RMSD="$ROOT/06_md_simulation/analysis/results/RMSD"
OUT="$ROOT/06_md_simulation/analysis/results/RMSD_NoPBC"

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

x="${SYSTEMS[$((SLURM_ARRAY_TASK_ID-1))]}"

TPR="$BASE/$x/Production/md_200ns.tpr"
XTC="$RMSD/${x}_nojump.xtc"
NDX="$RMSD/${x}.ndx"

RMSDOUT="$OUT/${x}_Ligand_RMSD_NoPBC.xvg"
SUMMARY="$OUT/${x}_Ligand_RMSD_NoPBC_Summary.tsv"

echo "SYSTEM: $x"

printf "4\n13\n" | \
gmx rms \
    -s "$TPR" \
    -f "$XTC" \
    -n "$NDX" \
    -o "$RMSDOUT" \
    -nopbc \
    -fit rot+trans \
    -tu ns

awk -v sys="$x" '
BEGIN{
    OFS="\t"
}
!/^[@#]/{
    n++
    s+=$2
    ss+=$2*$2

    if(n==1 || $2<min) min=$2
    if(n==1 || $2>max) max=$2
}
END{
    mean=s/n
    sd=sqrt((ss-s*s/n)/(n-1))

    print "System","N","Mean_RMSD_nm","SD_nm","Min_nm","Max_nm"
    printf "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
           sys,n,mean,sd,min,max
}' "$RMSDOUT" > "$SUMMARY"

echo "Finished: $x"
cat "$SUMMARY"


