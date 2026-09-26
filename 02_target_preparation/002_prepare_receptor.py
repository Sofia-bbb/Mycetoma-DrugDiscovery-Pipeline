#!/usr/bin/env python3

"""
003_prepare_receptor.py

Mycetoma Drug Discovery Platform

Purpose
-------
Prepare the validated MmSCD receptor for AutoDock Vina.

Input
-----
02_target_preparation/receptor/MmSCD_CE_superposed_on_2STD.pdb

Outputs
-------
06_Benchmark_Docking/Protein/MmSCD_receptor.pdb
06_Benchmark_Docking/Protein/MmSCD_receptor.pdbqt
06_Benchmark_Docking/Validation/Receptor_Preparation_Report.tsv
"""
from pathlib import Path
import subprocess
import csv
import shutil

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "02_target_preparation"

SOURCE = ROOT / "02_target_preparation" / "receptor" / "MmSCD_CE_superposed_on_2STD.pdb"

PROTEIN_DIR = PROJECT / "receptor"
VALIDATION_DIR = PROJECT / "validation"
PDBQT_FILE = PROTEIN_DIR / "MmSCD_receptor.pdbqt"

PROTEIN_DIR.mkdir(parents=True, exist_ok=True)
VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

PDB_FILE = PROTEIN_DIR / "MmSCD_receptor.pdb"
PDBQT_FILE = PROTEIN_DIR / "MmSCD_receptor.pdbqt"

REPORT = VALIDATION_DIR / "Receptor_Preparation_Report.tsv"

print("Preparing receptor...")

shutil.copy2(SOURCE, PDB_FILE)

output_prefix = PROTEIN_DIR / "MmSCD_receptor"

subprocess.run(
    [
        "mk_prepare_receptor.py",
        "--read_pdb",
        str(PDB_FILE),
        "-o",
        str(output_prefix),
        "-p",
    ],
    check=True,
)

with REPORT.open("w", newline="", encoding="utf-8") as handle:

    writer = csv.writer(handle, delimiter="\t")

    writer.writerow(
        [
            "Input",
            "Output",
            "Status",
        ]
    )

    writer.writerow(
        [
            PDB_FILE.name,
            PDBQT_FILE.name,
            "Prepared",
        ]
    )

print("\n========================================")
print("Receptor preparation completed.")
print(f"PDBQT : {PDBQT_FILE}")
print("========================================")


