#!/bin/bash
#SBATCH --job-name=NVT_CORR_MD11
#SBATCH --array=1-11
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=06_md_simulation/logs/NVT_%A_%a.out
#SBATCH --error=06_md_simulation/logs/NVT_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"

SYSTEM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" \
"$ROOT/06_md_simulation/MD11_corrected.txt")

SYS="$BASE/$SYSTEM"
NVT="$SYS/NVT"

mkdir -p "$NVT" "$ROOT/06_md_simulation/logs"

module purge
module load StdEnv/2023
module load gromacs/2024.6

echo "============================================"
echo "SYSTEM: $SYSTEM"
echo "CORRECTED NVT EQUILIBRATION"
echo "START: $(date)"
echo "============================================"

# ------------------------------------------------------------
# NVT parameters: same protocol as original MD
# ------------------------------------------------------------

cat > "$NVT/nvt.mdp" <<EOF
title                   = ${SYSTEM}_NVT

integrator              = md
dt                      = 0.002
nsteps                  = 50000

nstxout-compressed      = 500
nstenergy               = 500
nstlog                  = 500

continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
nstlist                 = 20

coulombtype             = PME
rcoulomb                = 1.0

vdwtype                 = Cut-off
rvdw                    = 1.0

tcoupl                  = V-rescale
tc-grps                 = Protein_MOL Water_and_ions
tau_t                   = 0.1 0.1
ref_t                   = 300 300

pcoupl                  = no

pbc                     = xyz

DispCorr                 = EnerPres

gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
EOF

cd "$NVT"

# ------------------------------------------------------------
# Build index groups
# ------------------------------------------------------------

printf \
"1 | 13\nname 18 Protein_MOL\n15 | 14\nname 19 Water_and_ions\nq\n" | \
gmx make_ndx \
-f ../Minimization/em.gro \
-o index.ndx \
>/dev/null 2>&1

grep -q '\[ Protein_MOL \]' index.ndx || {
    echo "ERROR: Protein_MOL group missing"
    exit 1
}

grep -q '\[ Water_and_ions \]' index.ndx || {
    echo "ERROR: Water_and_ions group missing"
    exit 1
}

echo "===== GROMPP ====="

gmx grompp \
-f nvt.mdp \
-c ../Minimization/em.gro \
-r ../Minimization/em.gro \
-p ../Solvation/topol.top \
-n index.ndx \
-o nvt.tpr \
-maxwarn 0

echo
echo "===== MDRUN ====="

gmx mdrun \
-deffnm nvt \
-ntomp ${SLURM_CPUS_PER_TASK} \
-v

test -s nvt.gro || {
    echo "ERROR: nvt.gro missing"
    exit 1
}

# ------------------------------------------------------------
# Post-NVT active-site QC
# ------------------------------------------------------------

echo
echo "===== POST-NVT ACTIVE-SITE QC ====="

printf \
"r MOL\nname 18 Ligand\nr 28 | r 48 | r 83 | r 109 | r 128 | r 130\nname 19 ActiveSite\nq\n" | \
gmx make_ndx \
-f nvt.tpr \
-o qc.ndx \
>/dev/null 2>&1

printf "18\n19\n" | \
gmx mindist \
-s nvt.tpr \
-f nvt.gro \
-n qc.ndx \
-od postNVT_active_distance.xvg \
-pbc yes \
>/dev/null 2>&1

DIST=$(awk '!/^[@#]/{print $2; exit}' postNVT_active_distance.xvg)

python3 - "$DIST" <<'PY'
import sys

d_nm=float(sys.argv[1])
d_A=d_nm*10

print(f"Post-NVT minimum distance = {d_A:.2f} Å")

if d_A <= 5.0:
    print("POST-NVT QC = PASS")
else:
    print("POST-NVT QC = FAIL")
    raise SystemExit(1)
PY

echo
echo "===== $SYSTEM NVT COMPLETE ====="
echo "End: $(date)"


