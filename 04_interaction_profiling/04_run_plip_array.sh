#!/bin/bash
#SBATCH --job-name=PLIP1188
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=04_interaction_profiling/Logs/plip_%A_%a.out
#SBATCH --error=04_interaction_profiling/Logs/plip_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load python/3.11

cd "$SLURM_SUBMIT_DIR"
source .venv/bin/activate

ROOT="04_interaction_profiling"
COMPLEXDIR="$ROOT/Complexes"
OUTROOT="$ROOT/PLIP_Results"
LIST="$ROOT/complex_list.txt"

mkdir -p "$OUTROOT" "$ROOT/Logs"

complex=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${complex:-}" ]]; then
    echo "ERROR: no complex for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

name=$(basename "$complex" _complex.pdb)
out="$OUTROOT/$name"

mkdir -p "$out"

# Skip safely if XML already exists
if ls "$out"/*report.xml >/dev/null 2>&1; then
    echo "Already completed: $name"
    exit 0
fi

echo "Running PLIP: $name"

plip \
  -f "$complex" \
  -o "$out" \
  -x \
  -t \
  --silent

if ls "$out"/*report.xml >/dev/null 2>&1; then
    echo "SUCCESS: $name"
else
    echo "FAILED: $name"
    exit 1
fi

