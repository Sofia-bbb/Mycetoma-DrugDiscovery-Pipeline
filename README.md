
# Reproducible Structure-Based Drug-Discovery Pipeline

This repository provides the computational workflow developed for systematic screening, prioritization, and molecular-dynamics-based evaluation of compounds from large chemical libraries.

The workflow integrates library preparation, target and binding-site assessment, molecular docking, protein–ligand interaction profiling, medicinal-chemistry and ADMET assessment, integrated candidate prioritization, molecular dynamics simulation, and MM/GBSA binding-energy estimation.

The pipeline is demonstrated using natural quinone compounds targeting scytalone dehydratase from *Madurella mycetomatis* (MmSCD) as a case study.

The outputs of this workflow represent **computational prioritization predictions** and should not be interpreted as evidence of biochemical inhibition, antifungal activity, or therapeutic efficacy.

## Workflow

The computational workflow includes the following sequential operations, grouped into seven documented computational stages:

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
12. Final qualitative computational candidate prioritization
```

## Repository Structure

```text
01_library_preparation/     Compound retrieval, curation, and 3D generation
02_target_preparation/      Target assessment, validation, and binding-site preparation
03_docking/                 Ligand preparation and molecular docking
04_interaction_profiling/    PLIP interaction analysis
05_prioritization/          Medicinal-chemistry and ADMET prioritization
06_md_simulation/           MD preparation, production, and analysis
07_mmgbsa/                 MM/GBSA binding-energy estimation
config/                    Project configuration
data/                      Data directories and any files provided therein
docs/                      Workflow documentation
examples/                  Representative ligand input and docking outputs
README.md
requirements.txt
```

The presence of a directory does not imply that all inputs, intermediate files, or generated outputs from the original study are archived in this repository.

## Case Study

For the MmSCD case study, the workflow progressed from 738,823 COCONUT records to 8,160 quinone hits and 8,093 unique compounds. A total of 8,044 compounds successfully underwent 3D structure generation, 8,023 were prepared for docking, and 8,021 produced valid docking poses.

Docking and spatial filtering retained 1,188 candidates, and PLIP interaction filtering reduced this set to 119 compounds. Of these, 88 passed the medicinal-chemistry progression criteria and underwent ADMET assessment. Nine compounds classified as `ADVANCE`, together with two reference inhibitors, underwent 200-ns molecular dynamics simulations. An integrated qualitative assessment of the resulting computational profiles prioritized five compounds for subsequent experimental evaluation.

Final prioritization was based on the overall balance of ligand RMSD, hydrogen-bond occupancy, binding-site contact persistence, MM/GBSA estimates, QED, and predicted ADMET properties rather than a single metric. The remaining four compounds were assigned lower relative priority; they were not classified as inactive.

## Software Requirements

The following software and versions are documented in the study methods. Versions not documented for individual Python packages are intentionally left unspecified rather than estimated.

```text
Python 3.11.4
├── RDKit (version not documented)
├── Biopython (version not documented)
├── Requests (version not documented)
└── PyYAML (version not documented)

Additional software
├── R 4.6.1
└── GNU Bash 5.1.16

Molecular docking
├── AutoDock Vina 1.2.5
└── Meeko 0.7.1

Protein–ligand interaction analysis
└── PLIP 3.0.1

ADMET prediction
└── ADMET-AI (version not documented)

Molecular dynamics
└── GROMACS 2024.6

MM/GBSA analysis and parameterization
├── AmberTools 25.0
├── MMPBSA.py
└── cpptraj

Computational environment
└── SLURM-based high-performance computing (HPC)
```

**Data resources:** COCONUT natural-product database (case-study input: `coconut_sdf_2d_lite-08-2026.sdf`), AlphaFold Protein Structure Database, UniProt, Protein Data Bank (PDB; reference structure 2STD), and PubChem (reference-inhibitor structures). These are external resources; access and versions of the underlying databases should be checked when reproducing the analyses.

Python package names are listed in [`requirements.txt`](requirements.txt); this file is not a fully version-locked environment specification. Cluster-specific paths, software modules, and SLURM resources must be adapted to the target computing environment.

## Installation and Execution

Clone the repository:

```bash
git clone https://github.com/Sofia-bbb/Mycetoma-DrugDiscovery-Pipeline.git
cd Mycetoma-DrugDiscovery-Pipeline
```

Prepare a compatible Python environment and install the listed Python dependencies as appropriate for that environment:

```bash
python -m pip install -r requirements.txt
```

The installation command above installs Python packages only. AutoDock Vina, Meeko, PLIP, GROMACS, AmberTools, and the required HPC tools must also be available in compatible environments for the corresponding workflow stages. The scripts are **not** presented as a single-command, end-to-end installation or execution system.

The seven documented computational stages, in execution order, are:

1. `01_library_preparation/`
2. `02_target_preparation/`
3. `03_docking/`
4. `04_interaction_profiling/`
5. `05_prioritization/`
6. `06_md_simulation/`
7. `07_mmgbsa/`

See [`docs/workflow.md`](docs/workflow.md) for reported parameters, filtering criteria, and stage descriptions. Inspect the scripts and their expected file paths and inputs before execution. Several stages use SLURM scripts and require environment-specific adjustments; computationally intensive work should be submitted to compute nodes rather than run on HPC login nodes.

The receptor used for docking is provided at [`02_target_preparation/receptor/MmSCD_receptor.pdbqt`](02_target_preparation/receptor/MmSCD_receptor.pdbqt), and the associated docking configuration is at [`02_target_preparation/receptor/vina_config.txt`](02_target_preparation/receptor/vina_config.txt). The [`examples/`](examples/) directory provides representative ligand inputs and docking outputs.

ADMET predictions were obtained using **ADMET-AI** and then integrated into the prioritization workflow. The ADMET prediction step and import of prediction outputs may require external tool execution and matching the expected input/output schema; this repository should not be assumed to provide a fully automated ADMET-AI service invocation. Likewise, final selection of five of the nine simulated candidates was an **integrated qualitative assessment**, not an automated hard-threshold or weighted-ranking algorithm.

## Reproducibility and Data Availability

The repository contains curated scripts and configuration files corresponding to the reported computational workflow. Historical, exploratory, and superseded scripts are not included. Large docking outputs, molecular-dynamics trajectories, temporary files, and cluster logs are excluded; therefore, this repository is **not a complete archive of all generated results**.

The repository provides selected receptor structures, docking configurations, molecular-dynamics ligand inputs, and representative example files. Some original study inputs and intermediate or output datasets may need to be obtained separately or regenerated. The documentation records the workflow and principal study parameters but does not establish that every stage can be run without adaptation on a fresh system. See [`docs/workflow.md`](docs/workflow.md) for stage-level documentation. Reproducibility of the complete numerical results depends on the original inputs, software and parameter versions, environment configuration, and any documented manual processing steps.

## Authors and Affiliations

**Sofia Ali**  
Department of Chemistry and Biochemistry  
University of Northern British Columbia (UNBC)  
Prince George, British Columbia, Canada

**Moaaz M. Saadaldin**  
Department of Bioinformatics and Biostatistics  
National University Biomedical Research Institute  
National University-Sudan  
Khartoum, Sudan

**Jacob Walsh**
Department of Chemistry and Biochemistry  
University of Northern British Columbia (UNBC)
Prince George, British Columbia, Canada
