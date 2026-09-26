#!/bin/bash
#SBATCH --job-name=NQLMM
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=01_library_preparation/logs/NQLMM_%j.out

cd $SLURM_SUBMIT_DIR

module load StdEnv/2023
module load python/3.11
module load rdkit/2025.09.4

source .venv/bin/activate

python 01_library_preparation/005_retrieve_natural_quinones.py


