#!/usr/bin/env python3

"""
NQL-MM 3D Conformer Generation

Project:
Mycetoma Drug Discovery Platform

Purpose:
Generate one optimized 3D conformer for every molecule in the
curated NQL-MM master library.

Method:
1. Read NQL-MM_Master.sdf
2. Add explicit hydrogens
3. Generate a conformer with ETKDGv3
4. Optimize with MMFF94
5. Use UFF when MMFF parameters are unavailable
6. Write optimized structures and QC reports

Notes:
- Progress is printed every 100 molecules.
- Output is flushed immediately for SLURM monitoring.
- Geometry optimization is limited to 200 iterations for
  high-throughput library preparation.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

INPUT_SDF = DATA / "processed" / "NQL-MM_Master.sdf"

OUTPUT_DIR = DATA / "processed"
REPORT_DIR = ROOT / "01_library_preparation" / "reports"
VALIDATION_DIR = ROOT / "01_library_preparation" / "validation"

OUTPUT_SDF = OUTPUT_DIR / "NQL-MM_3D.sdf"
REPORT_FILE = REPORT_DIR / "Optimization_Report.tsv"
QC_FILE = VALIDATION_DIR / "Geometry_QC.tsv"

for directory in [
    OUTPUT_DIR,
    REPORT_DIR,
    VALIDATION_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)


def generate_conformer(
    mol: Chem.Mol,
    random_seed: int,
) -> tuple[Chem.Mol | None, str, int, float]:
    """
    Generate and optimize one 3D conformer.

    Returns:
        prepared molecule or None
        optimization method
        optimization status code
        elapsed time in seconds
    """

    molecule_start = time.time()

    working_mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.useRandomCoords = True
    params.enforceChirality = True

    embed_status = AllChem.EmbedMolecule(
        working_mol,
        params,
    )

    if embed_status != 0:
        elapsed = time.time() - molecule_start
        return None, "Embedding_Failed", embed_status, elapsed

    if AllChem.MMFFHasAllMoleculeParams(working_mol):
        status = AllChem.MMFFOptimizeMolecule(
            working_mol,
            mmffVariant="MMFF94",
            maxIters=200,
        )

        elapsed = time.time() - molecule_start
        return working_mol, "MMFF94", status, elapsed

    try:
        status = AllChem.UFFOptimizeMolecule(
            working_mol,
            maxIters=200,
        )

        elapsed = time.time() - molecule_start
        return working_mol, "UFF", status, elapsed

    except Exception:
        elapsed = time.time() - molecule_start
        return None, "Optimization_Failed", -1, elapsed


def main() -> int:
    """Run 3D conformer generation for the full NQL-MM library."""

    if not INPUT_SDF.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_SDF}"
        )

    start_time = time.time()

    supplier = Chem.SDMolSupplier(
        str(INPUT_SDF),
        removeHs=False,
        sanitize=True,
    )

    writer = Chem.SDWriter(str(OUTPUT_SDF))

    total_records = 0
    valid_records = 0
    successful_records = 0
    mmff_records = 0
    uff_records = 0
    embedding_failures = 0
    optimization_failures = 0
    nonconverged_records = 0

    with REPORT_FILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as report_handle:

        report_writer = csv.writer(
            report_handle,
            delimiter="\t",
        )

        report_writer.writerow(
            [
                "NQL_ID",
                "Optimization_Method",
                "Optimization_Status",
                "Final_Status",
                "Elapsed_Seconds",
            ]
        )

        for mol in supplier:
            total_records += 1

            if mol is None:
                report_writer.writerow(
                    [
                        "",
                        "Invalid_Input",
                        "",
                        "Failed",
                        "",
                    ]
                )
                report_handle.flush()
                continue

            valid_records += 1

            if mol.HasProp("NQL_ID"):
                nql_id = mol.GetProp("NQL_ID")
            else:
                nql_id = f"NQLMM{valid_records:06d}"

            prepared_mol, method, status, elapsed = generate_conformer(
                mol,
                random_seed=20260801 + valid_records,
            )

            if prepared_mol is None:
                if method == "Embedding_Failed":
                    embedding_failures += 1
                else:
                    optimization_failures += 1

                report_writer.writerow(
                    [
                        nql_id,
                        method,
                        status,
                        "Failed",
                        f"{elapsed:.4f}",
                    ]
                )
                report_handle.flush()
                continue

            successful_records += 1

            if method == "MMFF94":
                mmff_records += 1
            elif method == "UFF":
                uff_records += 1

            if status != 0:
                nonconverged_records += 1

            prepared_mol.SetProp("NQL_ID", nql_id)
            prepared_mol.SetProp("3D_Method", method)
            prepared_mol.SetProp(
                "Optimization_Status",
                str(status),
            )
            prepared_mol.SetProp(
                "Optimization_Elapsed_Seconds",
                f"{elapsed:.4f}",
            )
            prepared_mol.SetProp("_Name", nql_id)

            writer.write(prepared_mol)
            writer.flush()

            report_writer.writerow(
                [
                    nql_id,
                    method,
                    status,
                    "Passed",
                    f"{elapsed:.4f}",
                ]
            )
            report_handle.flush()

            if total_records % 100 == 0:
                print(
                    f"Processed {total_records} molecules | "
                    f"Success={successful_records} | "
                    f"Failed="
                    f"{embedding_failures + optimization_failures}",
                    flush=True,
                )

    writer.close()

    elapsed_seconds = time.time() - start_time

    with QC_FILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        qc_writer = csv.writer(
            handle,
            delimiter="\t",
        )

        qc_writer.writerow(["Metric", "Value"])
        qc_writer.writerow(
            ["Input_records", total_records]
        )
        qc_writer.writerow(
            ["Valid_records", valid_records]
        )
        qc_writer.writerow(
            ["Successful_3D_records", successful_records]
        )
        qc_writer.writerow(
            ["MMFF94_records", mmff_records]
        )
        qc_writer.writerow(
            ["UFF_records", uff_records]
        )
        qc_writer.writerow(
            ["Embedding_failures", embedding_failures]
        )
        qc_writer.writerow(
            ["Optimization_failures", optimization_failures]
        )
        qc_writer.writerow(
            ["Nonconverged_records", nonconverged_records]
        )
        qc_writer.writerow(
            ["Elapsed_seconds", round(elapsed_seconds, 2)]
        )
        qc_writer.writerow(
            ["Output_SDF", str(OUTPUT_SDF)]
        )

    print("3D conformer generation completed.", flush=True)
    print(f"Input records: {total_records}", flush=True)
    print(
        f"Successful 3D structures: {successful_records}",
        flush=True,
    )
    print(f"MMFF94 optimized: {mmff_records}", flush=True)
    print(f"UFF optimized: {uff_records}", flush=True)
    print(
        f"Embedding failures: {embedding_failures}",
        flush=True,
    )
    print(
        f"Optimization failures: {optimization_failures}",
        flush=True,
    )
    print(
        f"Nonconverged optimizations: {nonconverged_records}",
        flush=True,
    )
    print(
        f"Elapsed time: {elapsed_seconds:.2f} seconds",
        flush=True,
    )
    print(f"Output: {OUTPUT_SDF}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


