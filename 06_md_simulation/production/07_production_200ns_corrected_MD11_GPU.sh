#!/bin/bash
#SBATCH --job-name=MD200_CORR11
#SBATCH --array=1-11
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --gpus=a100:1
#SBATCH --output=06_md_simulation/logs/MD200_%A_%a.out
#SBATCH --error=06_md_simulation/logs/MD200_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"

SYSTEM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" \
"$ROOT/06_md_simulation/MD11_corrected.txt")

SYS="$BASE/$SYSTEM"
PROD="$SYS/Production"

mkdir -p "$PROD" "$ROOT/06_md_simulation/logs"

module --force purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load cuda/12.6
module load gromacs/2024.6

echo "================================================"
echo "SYSTEM: $SYSTEM"
echo "CORRECTED 200-ns PRODUCTION MD"
echo "START: $(date)"
echo "NODE: $(hostname)"
echo "GPU:"
nvidia-smi -L
echo "================================================"

cat > "$PROD/md_200ns.mdp" <<EOF
integrator              = md
dt                      = 0.002
nsteps                  = 100000000

nstlog                  = 5000
nstenergy               = 5000
nstxout-compressed      = 5000
compressed-x-precision  = 1000

continuation            = yes
gen_vel                 = no

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

pcoupl                  = C-rescale
pcoupltype               = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

pbc                     = xyz

DispCorr                = EnerPres
EOF

cd "$PROD"

# ------------------------------------------------------------
# Coupling groups
# ------------------------------------------------------------

printf \
"1 | 13\nname 18 Protein_MOL\n15 | 14\nname 19 Water_and_ions\nq\n" | \
gmx make_ndx \
-f ../NPT/npt2.gro \
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

# ------------------------------------------------------------
# Final pre-production active-site QC
# ------------------------------------------------------------

printf \
"r MOL\nname 18 Ligand\nr 28 | r 48 | r 83 | r 109 | r 128 | r 130\nname 19 ActiveSite\nq\n" | \
gmx make_ndx \
-f ../NPT/npt2.tpr \
-o qc_preproduction.ndx \
>/dev/null 2>&1

printf "18\n19\n" | \
gmx mindist \
-s ../NPT/npt2.tpr \
-f ../NPT/npt2.gro \
-n qc_preproduction.ndx \
-od preproduction_active_distance.xvg \
-pbc yes \
>/dev/null 2>&1

DIST=$(awk '!/^[@#]/{print $2; exit}' preproduction_active_distance.xvg)

python3 - "$DIST" <<'PY'
import sys

d_nm=float(sys.argv[1])
d_A=d_nm*10

print(f"Pre-production minimum distance = {d_A:.2f} Å")

if d_A <= 5.0:
    print("PRE-PRODUCTION QC = PASS")
else:
    print("PRE-PRODUCTION QC = FAIL")
    raise SystemExit(1)
PY

# ------------------------------------------------------------
# Production TPR
# ------------------------------------------------------------

echo
echo "===== GROMPP PRODUCTION ====="

gmx grompp \
-f md_200ns.mdp \
-c ../NPT/npt2.gro \
-t ../NPT/npt2.cpt \
-p ../Solvation/topol.top \
-n index.ndx \
-o md_200ns.tpr \
-maxwarn 0

# ------------------------------------------------------------
# 200 ns production
# ------------------------------------------------------------

echo
echo "===== STARTING 200 NS MD ====="

gmx mdrun \
-deffnm md_200ns \
-ntmpi 1 \
-ntomp ${SLURM_CPUS_PER_TASK} \
-nb gpu \
-pme gpu \
-cpt 15 \
-maxh 167.5 \
-v

STATUS=$?

echo
echo "================================================"
echo "SYSTEM: $SYSTEM"
echo "END: $(date)"
echo "EXIT STATUS: $STATUS"
echo "================================================"

exit $STATUS


