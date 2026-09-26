#!/bin/bash
#SBATCH --job-name=MmSCD_CorrectedVS
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=03_docking/logs/slurm_%A_%a.out
#SBATCH --error=03_docking/logs/slurm_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load autodock-vina/1.2.6

ROOT="$SLURM_SUBMIT_DIR"

RECEPTOR="$ROOT/02_target_preparation/receptor/MmSCD_receptor.pdbqt"
LIST="$ROOT/03_docking/ligand_list.txt"

OUTDIR="$ROOT/03_docking/poses"
LOGDIR="$ROOT/03_docking/logs"

mkdir -p "$OUTDIR" "$LOGDIR"

# One ligand per SLURM array task
LIGAND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${LIGAND:-}" ]]; then
    echo "ERROR: No ligand found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

NAME=$(basename "$LIGAND" .pdbqt)

OUT="$OUTDIR/${NAME}_docked.pdbqt"
LOG="$LOGDIR/${NAME}_vina.log"

echo "============================================"
echo "Ligand:     $NAME"
echo "Task ID:    $SLURM_ARRAY_TASK_ID"
echo "Host:       $(hostname)"
echo "Start:      $(date)"
echo "============================================"

# Resume safely
if [[ -s "$OUT" ]] && grep -q '^MODEL' "$OUT"; then
    echo "Already completed: $NAME"
    exit 0
fi

vina \
    --receptor "$RECEPTOR" \
    --ligand "$LIGAND" \
    --center_x 0.0 \
    --center_y -1.5 \
    --center_z -3.5 \
    --size_x 22 \
    --size_y 22 \
    --size_z 22 \
    --exhaustiveness 32 \
    --num_modes 20 \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --out "$OUT" \
    > "$LOG" 2>&1

if [[ -s "$OUT" ]] && grep -q '^MODEL' "$OUT"; then
    echo "SUCCESS: $NAME"
else
    echo "FAILED: $NAME"
    exit 1
fi

echo "End: $(date)"


