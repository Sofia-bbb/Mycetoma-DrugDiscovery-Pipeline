# Workflow Documentation

This document describes the reproducible computational workflow implemented in this repository.

The pipeline is organized into seven stages:

1. Library preparation
2. Target preparation
3. Molecular docking
4. Protein–ligand interaction profiling
5. Candidate prioritization
6. Molecular dynamics simulation
7. MM/GBSA analysis

The repository contains the cleaned final workflow corresponding to the reported analysis. Historical, exploratory, and superseded scripts are not included.

## 1. Library Preparation

The starting chemical library is the COCONUT natural-products database. The workflow performs:

- inventory of the source SDF library;
- retrieval of quinone-containing compounds;
- removal of duplicate structures and assignment of internal NQLMM identifiers;
- generation of 3D molecular structures using RDKit;
- MMFF94 geometry optimization, with UFF used as fallback when required;
- explicit recording of unsuccessful structure-generation attempts.

In the MmSCD case study, 738,823 COCONUT records yielded 8,160 quinone hits and 8,093 unique compounds. Successful 3D structures were generated for 8,044 compounds.

Scripts are located in `01_library_preparation/`.

## 2. Target Preparation

The target-preparation workflow includes:

- retrieval of the MmSCD sequence and AlphaFold structural model;
- assessment of AlphaFold residue-level confidence (pLDDT);
- retrieval of the experimental reference SCD structure (PDB 2STD);
- pairwise sequence alignment and mapping of reference catalytic and binding-site residues to MmSCD;
- structural superposition of MmSCD and the experimental reference structure;
- identification and mapping of the reference ligand-binding pocket;
- preparation of the final receptor structure for molecular docking.

The MmSCD target is encoded by `MMYC01_208144` and corresponds to UniProt accession `A0A175VTY3`. The experimental reference structure is PDB `2STD`.

Scripts and reference structures are located in `02_target_preparation/`.

## 3. Molecular Docking

Ligands successfully converted to docking-ready PDBQT format were screened against the prepared MmSCD receptor using AutoDock Vina 1.2.5.

Final docking parameters were:

- grid center: `(0.0, -1.5, -3.5)`;
- grid size: `22 × 22 × 22 Å`;
- exhaustiveness: `32`;
- number of poses: `20`;
- CPU allocation: `4` per ligand.

A total of 8,023 ligands were prepared for docking, of which 8,021 produced valid docking poses.

Candidates were retained after docking and spatial quality control when they had a Vina score ≤ −10.0 kcal/mol and at least 3 of the 6 predefined MmSCD binding-site residues within 5 Å of the docked ligand. This retained 1,188 candidates.

Scripts and docking configuration are located in `03_docking/`.

## 4. Protein–Ligand Interaction Profiling

Protein–ligand interactions were characterized using PLIP 3.0.1.

The top-ranked docking pose was extracted for each retained ligand and combined with the receptor to generate protein–ligand complexes for interaction analysis.

Candidates were advanced when they retained the docking criterion (Vina score ≤ −10.0 kcal/mol) and interacted with at least 4 of the 6 predefined MmSCD binding-site residues.

This interaction-based filtering reduced the 1,188 docking/spatial candidates to 119 candidates for subsequent medicinal-chemistry and ADMET assessment.

Scripts are located in `04_interaction_profiling/`.

## 5. Candidate Prioritization

The 119 interaction-filtered candidates underwent sequential medicinal-chemistry and ADMET assessment.

Drug-likeness assessment included Lipinski, Veber, and additional physicochemical criteria. PAINS, Brenk, quinone-reactivity, and developability alerts were retained as descriptive information rather than treated as direct evidence of biological activity or inactivity.

After medicinal-chemistry filtering, 88 candidates proceeded to ADMET assessment. ADMET properties were integrated using a directional scoring framework, including evaluation of AMES, clinical toxicity, drug-induced liver injury, and hERG-related predictions.

Nine candidates classified as `ADVANCE` proceeded to molecular dynamics simulation together with the two reference compounds.

Scripts are located in `05_prioritization/`.

## 6. Molecular Dynamics Simulation

Nine prioritized candidates and two reference compounds were subjected to molecular dynamics simulation using GROMACS 2024.6.

Protein parameters used the ff14SB force field, ligand parameters used GAFF2 with AM1-BCC charges, and systems were solvated with TIP3P water using a 12 Å buffer and neutralizing ions.

Following energy minimization, systems underwent 100-ps NVT equilibration, 100-ps NPT equilibration, and an additional 500-ps NPT equilibration. Production simulations were then performed for 200 ns using a 2-fs timestep at 300 K and 1 bar.

Trajectory analyses included protein RMSD, RMSF, radius of gyration, SASA, hydrogen bonding, ligand RMSD, and binding-site contact persistence. A protein–ligand contact cutoff of ≤0.45 nm was used for contact-persistence analysis.

Scripts are located in `06_md_simulation/`.

## 7. MM/GBSA Analysis

MM/GBSA binding-energy estimates were calculated using MMPBSA.py from AmberTools 25.0.

For each complex, 1,001 trajectory snapshots were sampled from the final 100 ns of the 200-ns production trajectory. Calculations used the GB model `igb=5` with a salt concentration of 0.150 M. Entropic contributions were not included.

MM/GBSA values are interpreted as comparative computational binding-energy estimates rather than experimental binding affinities.

Scripts are located in `07_mmgbsa/`.

## HPC and Portability Notes

Core Python processing scripts use repository-relative paths and are intended to be portable across compatible environments.

SLURM scripts are provided for computationally intensive docking, PLIP, molecular dynamics, and MM/GBSA stages. Resource requests and software-module commands may require adjustment for the user's HPC environment.

Large docking outputs, molecular dynamics trajectories, temporary files, and cluster logs are intentionally excluded from the GitHub repository.
