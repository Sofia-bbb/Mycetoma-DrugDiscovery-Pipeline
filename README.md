# Reproducible Structure-Based Drug-Discovery Pipeline

This repository provides the computational workflow developed for systematic screening, prioritization, and molecular-dynamics-based evaluation of compounds from large chemical libraries.

The workflow integrates library preparation, target and binding-site assessment, molecular docking, protein–ligand interaction profiling, medicinal-chemistry and ADMET assessment, integrated candidate prioritization, molecular dynamics simulation, and MM/GBSA binding-energy estimation.

The pipeline is demonstrated using natural quinone compounds targeting scytalone dehydratase from *Madurella mycetomatis* (MmSCD) as a case study.

The outputs of this workflow represent **computational prioritization predictions** and should not be interpreted as evidence of biochemical inhibition, antifungal activity, or therapeutic efficacy.

## Workflow

The complete computational workflow comprised the following major stages:

```text
1. Natural-product library retrieval and inventory
   ↓
2. Quinone identification and library curation
   ↓
3. 3D structure generation and ligand preparation
   ↓
4. Target structure preparation and validation
   ↓
5. Binding-site mapping and receptor preparation
   ↓
6. Molecular docking and spatial quality control
   ↓
7. Protein–ligand interaction profiling
   ↓
8. Drug-likeness, medicinal-chemistry, and ADMET assessment
   ↓
9. Integrated candidate prioritization
   ↓
10. Molecular dynamics simulation
    ├── System preparation
    ├── Energy minimization
    ├── NVT/NPT equilibration
    ├── 200-ns production MD
    └── Trajectory analysis
   ↓
11. MM/GBSA binding-energy estimation
   ↓
12. Final computational candidate prioritization
## Repository Structure

```text
01_library_preparation/     Compound retrieval, curation, and 3D generation
02_target_preparation/      Target validation, binding-site mapping, and receptor preparation
03_docking/                 Ligand preparation and molecular docking
04_interaction_profiling/   PLIP interaction analysis
05_prioritization/          Drug-likeness, medicinal chemistry, and ADMET prioritization
06_md_simulation/           Molecular dynamics preparation, production, and analysis
07_mmgbsa/                  MM/GBSA binding-energy estimation
config/                     Project configuration
data/                       Input and processed data structure
examples/                   Small example inputs and outputs
```

## Case Study

For the MmSCD case study, the workflow progressed from 738,823 COCONUT records to 8,160 quinone hits and 8,093 unique compounds. A total of 8,044 compounds successfully underwent 3D structure generation, 8,023 were prepared for docking, and 8,021 produced valid docking poses.

Docking and spatial filtering retained 1,188 candidates. PLIP interaction filtering reduced this set to 119 candidates, followed by medicinal-chemistry and ADMET assessment and molecular-dynamics-based prioritization.

## Software Requirements

The workflow uses the following principal software:

- Python 3
- RDKit
- Biopython
- AutoDock Vina 1.2.5
- PLIP 3.0.1
- GROMACS 2024.6
- AmberTools 25.0
- Meeko

Python dependencies are listed in `requirements.txt`. Cluster-specific execution details and SLURM examples are documented separately.
