#!/usr/bin/env python3

"""
002_prepare_screening_ligands.py

Mycetoma Drug Discovery Platform

Purpose
-------
Prepare the complete NQL-MM 3D library for AutoDock Vina using
Meeko. Each valid molecule is written to an individual SDF file
and converted into a docking-ready PDBQT file.

Input
-----
05_NQL-MM/3D/NQL-MM_3D.sdf

Outputs
-------
08_Virtual_Screening/Ligands/SDF/*.sdf
08_Virtual_Screening/Ligands/PDBQT/*.pdbqt
08_Virtual_Screening/Validation/
    Screening_Ligand_Preparation.tsv
    Screening_Ligand_Preparation_Summary.txt
    Failed_Ligands.tsv

Features
--------
- Preserves permanent NQL-MM identifiers.
- Supports restart after interruption.
- Skips valid PDBQT files already generated.
- Validates every output file.
- Records failed ligands without stopping the full workflow.
- Prints progress every 100 molecules.
"""

from __future__ import annotations

import csv
import subprocess
import time
from pathlib import Path

from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_SDF = PROJECT_ROOT / "data" / "processed" / "NQL-MM_3D.sdf"

MODULE_DIR = PROJECT_ROOT / "03_docking"

SDF_DIR = MODULE_DIR / "Ligands" / "SDF"
PDBQT_DIR = MODULE_DIR / "Ligands" / "PDBQT"
VALIDATION_DIR = MODULE_DIR / "Validation"

REPORT_FILE = (
    VALIDATION_DIR
    / "Screening_Ligand_Preparation.tsv"
)

FAILED_FILE = VALIDATION_DIR / "Failed_Ligands.tsv"

SUMMARY_FILE = (
    VALIDATION_DIR
    / "Screening_Ligand_Preparation_Summary.txt"
)

SDF_DIR.mkdir(parents=True, exist_ok=True)
PDBQT_DIR.mkdir(parents=True, exist_ok=True)
VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

PROGRESS_INTERVAL = 100
MINIMUM_PDBQT_SIZE = 100


def safe_identifier(text: str) -> str:
    """Convert a molecule identifier into a safe filename."""

    cleaned = "".join(
        character
        if character.isalnum() or character in {"-", "_"}
        else "_"
        for character in text.strip()
    )

    return cleaned or "UNNAMED"


def get_molecule_id(
    molecule: Chem.Mol,
    record_number: int,
) -> str:
    """
    Retrieve the permanent NQL-MM identifier when available.

    Property priority:
    1. NQL_ID
    2. _Name
    3. Generated fallback identifier
    """

    if molecule.HasProp("NQL_ID"):
        value = molecule.GetProp("NQL_ID").strip()

        if value:
            return safe_identifier(value)

    if molecule.HasProp("_Name"):
        value = molecule.GetProp("_Name").strip()

        if value:
            return safe_identifier(value)

    return f"NQLMM{record_number:07d}"


def write_single_sdf(
    molecule: Chem.Mol,
    output_file: Path,
    molecule_id: str,
) -> None:
    """Write one molecule to an individual SDF file."""

    molecule.SetProp("_Name", molecule_id)
    molecule.SetProp("NQL_ID", molecule_id)

    writer = Chem.SDWriter(str(output_file))

    if writer is None:
        raise RuntimeError(
            f"Could not create SDF writer for {output_file}"
        )

    writer.write(molecule)
    writer.close()

    if (
        not output_file.exists()
        or output_file.stat().st_size < 100
    ):
        raise RuntimeError(
            f"Individual SDF output is invalid: {output_file}"
        )


def prepare_with_meeko(
    input_sdf: Path,
    output_pdbqt: Path,
) -> tuple[str, str]:
    """
    Convert one ligand from SDF to PDBQT using Meeko.

    Returns
    -------
    status, message
    """

    command = [
        "mk_prepare_ligand.py",
        "-i",
        str(input_sdf),
        "-o",
        str(output_pdbqt),
    ]

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )

    if result.returncode != 0:
        message = (
            result.stderr.strip()
            or result.stdout.strip()
            or f"Meeko exit code {result.returncode}"
        )

        return "Failed", message

    if (
        not output_pdbqt.exists()
        or output_pdbqt.stat().st_size < MINIMUM_PDBQT_SIZE
    ):
        return (
            "Failed",
            "Meeko completed but no valid PDBQT was created.",
        )

    return "Prepared", ""


def main() -> None:
    if not INPUT_SDF.exists():
        raise FileNotFoundError(
            f"Final NQL-MM library not found: {INPUT_SDF}"
        )

    supplier = Chem.SDMolSupplier(
        str(INPUT_SDF),
        removeHs=False,
        sanitize=False,
    )

    start_time = time.perf_counter()

    total_records = 0
    valid_inputs = 0
    prepared = 0
    skipped = 0
    failed = 0
    invalid_inputs = 0

    report_rows: list[list[str]] = []
    failed_rows: list[list[str]] = []

    print("=" * 72)
    print("Preparing NQL-MM ligands for virtual screening")
    print(f"Input library : {INPUT_SDF}")
    print(f"Meeko output  : {PDBQT_DIR}")
    print("=" * 72)

    for record_number, molecule in enumerate(
        supplier,
        start=1,
    ):
        total_records += 1

        if molecule is None:
            invalid_inputs += 1
            failed += 1

            fallback_id = f"NQLMM{record_number:07d}"

            report_rows.append(
                [
                    str(record_number),
                    fallback_id,
                    "",
                    "",
                    "Invalid_Input",
                    "RDKit could not read the molecule.",
                    "",
                ]
            )

            failed_rows.append(
                [
                    fallback_id,
                    "Invalid_Input",
                    "RDKit could not read the molecule.",
                ]
            )

            continue

        valid_inputs += 1

        molecule_id = get_molecule_id(
            molecule,
            record_number,
        )

        sdf_file = SDF_DIR / f"{molecule_id}.sdf"
        pdbqt_file = PDBQT_DIR / f"{molecule_id}.pdbqt"

        molecule_start = time.perf_counter()

        # Restart support: keep a previously generated valid PDBQT.
        if (
            pdbqt_file.exists()
            and pdbqt_file.stat().st_size
            >= MINIMUM_PDBQT_SIZE
        ):
            skipped += 1

            elapsed = time.perf_counter() - molecule_start

            report_rows.append(
                [
                    str(record_number),
                    molecule_id,
                    str(sdf_file),
                    str(pdbqt_file),
                    "Skipped_Existing",
                    "",
                    f"{elapsed:.3f}",
                ]
            )

        else:
            try:
                if (
                    not sdf_file.exists()
                    or sdf_file.stat().st_size < 100
                ):
                    write_single_sdf(
                        molecule,
                        sdf_file,
                        molecule_id,
                    )

                status, message = prepare_with_meeko(
                    sdf_file,
                    pdbqt_file,
                )

                elapsed = (
                    time.perf_counter()
                    - molecule_start
                )

                if status == "Prepared":
                    prepared += 1
                else:
                    failed += 1

                    failed_rows.append(
                        [
                            molecule_id,
                            status,
                            message,
                        ]
                    )

                report_rows.append(
                    [
                        str(record_number),
                        molecule_id,
                        str(sdf_file),
                        str(pdbqt_file),
                        status,
                        message,
                        f"{elapsed:.3f}",
                    ]
                )

            except Exception as error:
                failed += 1

                elapsed = (
                    time.perf_counter()
                    - molecule_start
                )

                error_message = str(error)

                report_rows.append(
                    [
                        str(record_number),
                        molecule_id,
                        str(sdf_file),
                        str(pdbqt_file),
                        "Failed",
                        error_message,
                        f"{elapsed:.3f}",
                    ]
                )

                failed_rows.append(
                    [
                        molecule_id,
                        "Failed",
                        error_message,
                    ]
                )

        if total_records % PROGRESS_INTERVAL == 0:
            print(
                f"Processed {total_records:,} molecules | "
                f"Prepared={prepared:,} | "
                f"Skipped={skipped:,} | "
                f"Failed={failed:,}",
                flush=True,
            )

    total_elapsed = time.perf_counter() - start_time

    with REPORT_FILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "Record_Number",
                "NQL_ID",
                "Individual_SDF",
                "PDBQT_File",
                "Status",
                "Message",
                "Elapsed_Seconds",
            ]
        )

        writer.writerows(report_rows)

    with FAILED_FILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "NQL_ID",
                "Failure_Type",
                "Message",
            ]
        )

        writer.writerows(failed_rows)

    final_valid_pdbqt = sum(
        1
        for path in PDBQT_DIR.glob("*.pdbqt")
        if path.stat().st_size >= MINIMUM_PDBQT_SIZE
    )

    # Create the ordered ligand list used by the docking array.
    ligand_list_file = MODULE_DIR / "ligand_list.txt"
    valid_pdbqt_files = sorted(
        path for path in PDBQT_DIR.glob("*.pdbqt")
        if path.stat().st_size >= MINIMUM_PDBQT_SIZE
    )

    with ligand_list_file.open("w", encoding="utf-8") as handle:
        for path in valid_pdbqt_files:
            handle.write(str(path.resolve()) + "\n")

    with SUMMARY_FILE.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "NQL-MM Screening Ligand Preparation Summary\n"
        )
        handle.write(
            "===========================================\n\n"
        )
        handle.write(f"Input library: {INPUT_SDF}\n")
        handle.write(
            f"Total records: {total_records}\n"
        )
        handle.write(
            f"Valid input molecules: {valid_inputs}\n"
        )
        handle.write(
            f"Invalid input molecules: {invalid_inputs}\n"
        )
        handle.write(
            f"Newly prepared ligands: {prepared}\n"
        )
        handle.write(
            f"Existing ligands skipped: {skipped}\n"
        )
        handle.write(
            f"Failed preparations: {failed}\n"
        )
        handle.write(
            f"Final valid PDBQT files: {final_valid_pdbqt}\n"
        )
        handle.write(
            f"Elapsed seconds: {total_elapsed:.2f}\n"
        )
        handle.write(f"Detailed report: {REPORT_FILE}\n")
        handle.write(f"Failure report: {FAILED_FILE}\n")

    print("\n" + "=" * 72)
    print("Screening-ligand preparation finished.")
    print(f"Total records        : {total_records:,}")
    print(f"Valid inputs         : {valid_inputs:,}")
    print(f"Newly prepared       : {prepared:,}")
    print(f"Skipped existing     : {skipped:,}")
    print(f"Failed               : {failed:,}")
    print(f"Final valid PDBQT    : {final_valid_pdbqt:,}")
    print(f"Elapsed time         : {total_elapsed:.2f} seconds")
    print(f"Detailed report      : {REPORT_FILE}")
    print(f"Failure report       : {FAILED_FILE}")
    print("=" * 72)


if __name__ == "__main__":
    main()



