#!/bin/bash
#SBATCH --job-name=MMGBSA10
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-11
#SBATCH --output=07_mmgbsa/results/MMGBSA_%A_%a.out
#SBATCH --error=07_mmgbsa/results/MMGBSA_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load ambertools/25.0

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation"
ANALYSIS="$BASE/analysis/results"
MMPBSA="$ROOT/07_mmgbsa/results"

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

WORK="$MMPBSA/$x"
mkdir -p "$WORK"
cd "$WORK"

echo "=========================================="
echo "SYSTEM: $x"
echo "=========================================="

TOP="$BASE/systems/$x/Topology/${x}_corrected_unsolvated.prmtop"
XTC="$ANALYSIS/RMSD/${x}_nojump.xtc"

# ---------------------------------------------------------
# 1. Prepare 1001 snapshots from last 100 ns
# ---------------------------------------------------------

cat > prepare_snapshots.in <<EOF
parm $TOP
trajin $XTC 10001 20001 10
trajout ${x}_last100ns.nc netcdf
run
quit
EOF

cpptraj -i prepare_snapshots.in

# ---------------------------------------------------------
# 2. Prepare complex / receptor / ligand topologies
# Protein = residues 1-172
# Ligand  = residue 173
# ---------------------------------------------------------

cp "$TOP" complex.prmtop

ante-MMPBSA.py \
-p complex.prmtop \
-r receptor.prmtop \
-l ligand.prmtop \
-m ':1-172'

# ---------------------------------------------------------
# 3. MM/GBSA input
# ---------------------------------------------------------

cat > mmgbsa.in <<EOF
&general
  startframe=1,
  endframe=1001,
  interval=1,
  verbose=1,
/

&gb
  igb=5,
  saltcon=0.150,
/
EOF

# ---------------------------------------------------------
# 4. Run MM/GBSA
# ---------------------------------------------------------

MMPBSA.py -O \
-i mmgbsa.in \
-o FINAL_RESULTS_MMGBSA.dat \
-eo FINAL_RESULTS_MMGBSA.csv \
-cp complex.prmtop \
-rp receptor.prmtop \
-lp ligand.prmtop \
-y ${x}_last100ns.nc

echo
echo "MM/GBSA COMPLETE: $x"

grep -A15 "Differences (Complex - Receptor - Ligand)" \
FINAL_RESULTS_MMGBSA.dat || true


