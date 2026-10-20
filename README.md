
# Reproducible Structure-Based Drug-Discovery Pipeline

This repository provides a computational workflow for systematic screening, prioritization, and molecular-dynamics-based evaluation of compounds from large chemical libraries.

The framework integrates chemical-library preparation, target-structure assessment, molecular docking, protein–ligand interaction profiling, medicinal-chemistry and ADMET evaluation, candidate prioritization, molecular dynamics (MD) simulations, and MM/GBSA binding-energy estimation.

The workflow is demonstrated using natural quinone compounds targeting scytalone dehydratase from *Madurella mycetomatis* (MmSCD) as a case study.

The results represent computational predictions and candidate prioritization, not experimental evidence of biochemical inhibition, antifungal activity, or therapeutic efficacy.

## Workflow

The workflow comprises the following sequential operations. Implementation scripts are organized into seven numbered computational-stage directories, while final candidate prioritization is performed through an integrated qualitative assessment.

1. Natural-product library retrieval and inventory
2. Quinone identification and library curation
3. Three-dimensional structure generation and ligand preparation
4. Target-structure preparation and validation
5. Binding-site mapping and receptor preparation
6. Molecular docking and spatial quality control
7. Protein–ligand interaction profiling
8. Drug-likeness, medicinal-chemistry, and ADMET assessment
9. Integrated candidate prioritization
10. Molecular dynamics simulations and trajectory analysis
11. MM/GBSA binding-energy estimation
12. Final qualitative computational candidate prioritization

Detailed descriptions of the seven computational stages, associated parameters, and filtering criteria are provided in [`docs/workflow.md`](docs/workflow.md).

## Repository Structure

The repository is organized into seven computational stages, with supporting directories for configuration, documentation, and example data.

| Directory | Description |
|---|---|
| `01_library_preparation/` | Compound retrieval, curation, and 3D structure generation |
| `02_target_preparation/` | Target assessment, validation, and binding-site preparation |
| `03_docking/` | Ligand preparation and molecular docking |
| `04_interaction_profiling/` | Protein–ligand interaction analysis using PLIP |
| `05_prioritization/` | Medicinal-chemistry and ADMET-based prioritization |
| `06_md_simulation/` | MD system preparation, production simulations, and trajectory analysis |
| `07_mmgbsa/` | MM/GBSA binding-energy estimation |
| `config/` | Project configuration |
| `data/` | Available data directories and files |
| `docs/` | Workflow documentation |
| `examples/` | Representative ligand inputs and docking outputs |

Additional repository files include `README.md` and `requirements.txt`.

## Case Study

For the MmSCD case study, the workflow progressed from 738,823 COCONUT records to 8,160 quinone hits and 8,093 unique compounds. A total of 8,044 compounds successfully underwent 3D structure generation, 8,023 were prepared for docking, and 8,021 produced valid docking poses.

Docking and spatial quality control retained 1,188 candidates, and PLIP interaction filtering reduced this set to 119 compounds. Of these, 88 passed the medicinal-chemistry progression criteria and underwent ADMET assessment.

Nine compounds classified as `ADVANCE`, together with two reference inhibitors, carpropamid and fenoxanil, underwent 200-ns molecular dynamics simulations. An integrated qualitative assessment of their computational profiles prioritized five compounds for subsequent experimental evaluation:

- NQLMM000191
- NQLMM003691
- NQLMM002852
- NQLMM002327
- NQLMM002110

Final prioritization considered the overall balance of ligand RMSD, hydrogen-bond occupancy, binding-site contact persistence, MM/GBSA binding-energy estimates, quantitative estimate of drug-likeness (QED), and predicted ADMET properties rather than any single parameter.

The remaining four compounds were assigned lower relative priority and were not classified as inactive or unsuitable for future investigation.

## Software Requirements

The following software and versions are documented in the study methods. Versions that were not recorded are left unspecified rather than estimated.

| Software or package | Version | Application |
|---|---|---|
| Python | 3.11.4 | Workflow scripting and data processing |
| RDKit | Not documented | Chemical-library curation, molecular descriptors, and 3D preparation |
| Biopython | Not documented | Sequence alignment and structural analysis |
| Requests | Not documented | Data retrieval |
| PyYAML | Not documented | Configuration handling |
| R | 4.6.1 | Computational data processing |
| GNU Bash | 5.1.16 | Workflow execution and shell scripting |
| AutoDock Vina | 1.2.5 | Molecular docking |
| Meeko | 0.7.1 | Ligand and receptor preparation |
| PLIP | 3.0.1 | Protein–ligand interaction profiling |
| ADMET-AI | Not documented | ADMET prediction |
| GROMACS | 2024.6 | Molecular dynamics simulations |
| AmberTools | 25.0 | Molecular parameterization and MM/GBSA analysis |
| MMPBSA.py | AmberTools component | MM/GBSA calculations |
| cpptraj | AmberTools component | Molecular trajectory processing |

### Databases and Structural Resources

| Database or resource | Purpose |
|---|---|
| COCONUT | Natural-product chemical-library retrieval |
| AlphaFold Protein Structure Database | Predicted MmSCD structure |
| UniProt | Protein sequence retrieval |
| Protein Data Bank (PDB) | Experimental reference SCD structure |
| PubChem | Reference-inhibitor chemical structures |

### COCONUT Input Dataset

The original COCONUT SDF file used in this study was
`coconut_sdf_2d_lite-08-2026.sdf`, containing 738,823 molecular records.

To reproduce the library-preparation stage, obtain the corresponding
COCONUT dataset and place the SDF file at:

`data/raw/coconut_sdf_2d_lite-08-2026.sdf`

This large source dataset is not included in the repository.
Different COCONUT releases may produce different screening results.

The COCONUT source file used in the case study was `coconut_sdf_2d_lite-08-2026.sdf`. The experimental structural reference was PDB `2STD`.

The computational work was conducted using a SLURM-based high-performance computing environment, including the Narval system of the Digital Research Alliance of Canada.

Python package dependencies are listed in [`requirements.txt`](requirements.txt). The file is not a fully version-locked environment specification.

## Installation and Execution

### 1. Clone the repository

```bash
git clone https://github.com/Sofia-bbb/Mycetoma-DrugDiscovery-Pipeline.git
cd Mycetoma-DrugDiscovery-Pipeline
```

### 2. Prepare the Python environment

A compatible Python environment is required. The listed Python dependencies may be installed using:

```bash
python -m pip install -r requirements.txt
```

Additional scientific software, including AutoDock Vina, Meeko, PLIP, GROMACS, and AmberTools, must be available for their respective workflow stages.

Installing the Python requirements alone does not establish a complete computational environment.

### 3. Execute the computational stages

The implementation is organized in the following order:

| Stage | Directory |
|---|---|
| 1. Library preparation | `01_library_preparation/` |
| 2. Target preparation | `02_target_preparation/` |
| 3. Molecular docking | `03_docking/` |
| 4. Interaction profiling | `04_interaction_profiling/` |
| 5. Candidate prioritization | `05_prioritization/` |
| 6. Molecular dynamics simulations | `06_md_simulation/` |
| 7. MM/GBSA analysis | `07_mmgbsa/` |

Workflow parameters, filtering criteria, and stage descriptions are available in [`docs/workflow.md`](docs/workflow.md).

Scripts should be reviewed for their required input files, output locations, external dependencies, and environment-specific settings before execution.

Computationally intensive stages use SLURM scripts. Resource requests, software modules, and filesystem paths must be configured for the target HPC environment. Intensive computations should be submitted to compute nodes rather than executed on login nodes.

### 4. Reference Inputs and Examples

The repository provides selected prepared structures, docking configurations, and example files.

The prepared docking receptor is available at:

`02_target_preparation/receptor/MmSCD_receptor.pdbqt`

The docking configuration is available at:

`02_target_preparation/receptor/vina_config.txt`

Representative input and output files are provided in:

`examples/`

### 5. ADMET Assessment and Final Prioritization

ADMET predictions were generated using ADMET-AI and incorporated into the medicinal-chemistry and candidate-prioritization workflow.

The repository includes scripts for preparing ADMET inputs and integrating prediction outputs. It does not provide a documented automated invocation of ADMET-AI; prediction generation must therefore be handled separately when reproducing this stage.

Following molecular dynamics and MM/GBSA analysis, final selection of five compounds was based on an integrated qualitative assessment of complementary computational metrics.

This final decision was not implemented as an automated weighted-score or hard-threshold selection algorithm.

## Reproducibility and Data Availability

The repository contains curated scripts, configuration files, selected receptor structures, molecular-dynamics ligand inputs, and representative examples corresponding to the reported computational workflow.

Historical, exploratory, and superseded scripts, large docking outputs, molecular-dynamics trajectories, temporary files, and cluster logs are excluded.

Accordingly, the repository documents the computational implementation but does not constitute a complete archive of all original input datasets, intermediate files, or generated results.

Reproduction of the numerical results requires the relevant input datasets, compatible software versions, reported computational parameters, environment-specific configuration, and applicable manual processing steps.

The workflow and its main computational settings are described in [`docs/workflow.md`](docs/workflow.md).

## Authors and Affiliations

**Sofia B Mohamed**
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
