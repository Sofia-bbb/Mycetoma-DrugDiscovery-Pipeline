#!/bin/bash
#SBATCH --job-name=EM_CORR_MD11
#SBATCH --array=1-11
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=06_md_simulation/logs/EM_%A_%a.out
#SBATCH --error=06_md_simulation/logs/EM_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$ROOT/06_md_simulation/systems"

SYSTEM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" \
"$ROOT/06_md_simulation/MD11_corrected.txt")

SYS="$BASE/$SYSTEM"
EM="$SYS/Minimization"

mkdir -p "$EM" "$ROOT/06_md_simulation/logs"

module purge
module load StdEnv/2023
module load gromacs/2024.6

echo "============================================"
echo "SYSTEM: $SYSTEM"
echo "CORRECTED ENERGY MINIMIZATION"
echo "START: $(date)"
echo "============================================"

cat > "$EM/em.mdp" <<'EOF'
integrator      = steep

emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000

cutoff-scheme   = Verlet
nstlist         = 20

coulombtype     = PME
rcoulomb        = 1.0

vdwtype         = Cut-off
rvdw            = 1.0

pbc             = xyz

nstenergy       = 100
EOF

cd "$EM"

echo "===== GROMPP ====="

gmx grompp \
-f em.mdp \
-c ../Solvation/${SYSTEM}_solvated.gro \
-p ../Solvation/topol.top \
-o em.tpr \
-maxwarn 0

echo
echo "===== MDRUN ====="

gmx mdrun \
-deffnm em \
-ntomp ${SLURM_CPUS_PER_TASK} \
-v

test -s em.gro || {
    echo "ERROR: em.gro missing"
    exit 1
}

echo
echo "===== POST-EM ACTIVE-SITE QC ====="

# Create explicit ligand and active-site groups.
printf \
"r MOL\nname 18 Ligand\nr 28 | r 48 | r 83 | r 109 | r 128 | r 130\nname 19 ActiveSite\nq\n" \
| gmx make_ndx \
-f em.tpr \
-o qc.ndx \
>/dev/null 2>&1

# Measure minimum ligand-active-site distance.
printf "18\n19\n" | \
gmx mindist \
-s em.tpr \
-f em.gro \
-n qc.ndx \
-od postEM_active_distance.xvg \
-pbc yes \
>/dev/null 2>&1

DIST=$(awk '!/^[@#]/{print $2; exit}' postEM_active_distance.xvg)

echo "Post-EM minimum distance = $DIST nm"

python3 - "$DIST" <<'PY'
import sys

d_nm=float(sys.argv[1])
d_A=d_nm*10

print(f"Post-EM minimum distance = {d_A:.2f} Å")

if d_A <= 5.0:
    print("POST-EM QC = PASS")
else:
    print("POST-EM QC = FAIL")
    raise SystemExit(1)
PY

echo
echo "===== $SYSTEM COMPLETE ====="
echo "End: $(date)"


