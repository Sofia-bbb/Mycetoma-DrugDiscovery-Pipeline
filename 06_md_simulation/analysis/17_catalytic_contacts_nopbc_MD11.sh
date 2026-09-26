#!/bin/bash
#SBATCH --job-name=CATNOPBC_MD11
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-11
#SBATCH --output=06_md_simulation/analysis/results/ActiveSite_NoPBC/CATNOPBC_%A_%a.out
#SBATCH --error=06_md_simulation/analysis/results/ActiveSite_NoPBC/CATNOPBC_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load gromacs/2024.6

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"
RMSD="$ROOT/06_md_simulation/analysis/results/RMSD"
OUT="$ROOT/06_md_simulation/analysis/results/ActiveSite_NoPBC"

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

echo "=================================================="
echo "SYSTEM: $x"
echo "NO-PBC ACTIVE-SITE CONTACT ANALYSIS"
echo "=================================================="

TPR="$BASE/$x/Production/md_200ns.tpr"
XTC="$RMSD/${x}_nojump.xtc"

SYSOUT="$OUT/$x"
mkdir -p "$SYSOUT"

NDX="$SYSOUT/${x}_active_site.ndx"

printf \
"r 28\nname 18 TYR28\nr 48\nname 19 TYR48\nr 83\nname 20 HIS83\nr 109\nname 21 HIS109\nr 128\nname 22 SER128\nr 130\nname 23 ASN130\nq\n" \
| gmx make_ndx \
    -f "$TPR" \
    -o "$NDX" \
    >/dev/null 2>&1

declare -A RESGROUP
RESGROUP[TYR28]=18
RESGROUP[TYR48]=19
RESGROUP[HIS83]=20
RESGROUP[HIS109]=21
RESGROUP[SER128]=22
RESGROUP[ASN130]=23

CUTOFF=0.45

SUMMARY="$SYSOUT/${x}_ActiveSite_Contact_Occupancy_NoPBC.tsv"

printf \
"System\tResidue\tCutoff_nm\tFrames_Total\tFrames_In_Contact\tOccupancy_Percent\tMean_Min_Distance_nm\tMin_Distance_nm\tMax_Distance_nm\n" \
> "$SUMMARY"

for residue in TYR28 TYR48 HIS83 HIS109 SER128 ASN130
do

    grp="${RESGROUP[$residue]}"
    DIST="$SYSOUT/${x}_${residue}_mindist_NoPBC.xvg"

    echo "Analyzing $x -> $residue"

    printf "13\n%s\n" "$grp" | \
    gmx mindist \
        -s "$TPR" \
        -f "$XTC" \
        -n "$NDX" \
        -od "$DIST" \
        -d "$CUTOFF" \
        -nopbc \
        -tu ns \
        >/dev/null 2>&1

    awk -v sys="$x" \
        -v residue="$residue" \
        -v cutoff="$CUTOFF" '
        !/^[@#]/ {
            n++
            d=$2
            sum+=d

            if(n==1 || d<min) min=d
            if(n==1 || d>max) max=d

            if(d<=cutoff) contact++
        }

        END {
            if(n>0) {
                printf "%s\t%s\t%.2f\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\n",
                       sys,
                       residue,
                       cutoff,
                       n,
                       contact,
                       100*contact/n,
                       sum/n,
                       min,
                       max
            }
        }
    ' "$DIST" >> "$SUMMARY"

done

echo "Finished: $x"
echo "Output: $SUMMARY"


