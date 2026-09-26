#!/usr/bin/env python3

"""
NQL-MM Library Curation Engine

Project:
Mycetoma Drug Discovery Platform

Library:
Natural Quinone Library for Madurella mycetomatis

Version:
1.0

Purpose:
Create the curated master quinone library for downstream
3D preparation, PDBQT conversion, and virtual screening.
"""

from pathlib import Path
import csv

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

INPUT_SDF = DATA / "processed" / "Natural_Quinones.sdf"

MASTER_DIR = DATA / "processed"
CURATED_DIR = DATA / "processed"
REPORTS_DIR = ROOT / "01_library_preparation" / "reports"
VALIDATION_DIR = ROOT / "01_library_preparation" / "validation"
STATISTICS_DIR = ROOT / "01_library_preparation" / "statistics"

for directory in [
    MASTER_DIR,
    CURATED_DIR,
    REPORTS_DIR,
    VALIDATION_DIR,
    STATISTICS_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)


MASTER_SDF = MASTER_DIR / "NQL-MM_Master.sdf"
MASTER_TSV = CURATED_DIR / "NQL-MM_Master.tsv"
MASTER_CSV = CURATED_DIR / "NQL-MM_Master.csv"
MASTER_SMI = CURATED_DIR / "NQL-MM_Master.smi"

QC_REPORT = VALIDATION_DIR / "Library_QC.tsv"
DUPLICATE_REPORT = VALIDATION_DIR / "Duplicate_Report.tsv"
STATISTICS_FILE = STATISTICS_DIR / "Library_Statistics.tsv"


OUTPUT_COLUMNS = [
    "NQL_ID",
    "Source_ID",
    "Source_Database",
    "Library_Version",
    "Quinone_Class",
    "Canonical_SMILES",
    "Isomeric_SMILES",
    "InChI",
    "InChIKey",
    "Molecular_Formula",
    "Molecular_Weight",
    "Exact_Mass",
    "LogP",
    "TPSA",
    "HBD",
    "HBA",
    "Rotatable_Bonds",
    "Heavy_Atoms",
    "Aromatic_Rings",
    "Total_Rings",
]


def get_property(mol: Chem.Mol, name: str) -> str:
    """Return an SDF property if present."""
    if mol.HasProp(name):
        return mol.GetProp(name).strip()
    return ""


def calculate_record(
    mol: Chem.Mol,
    nql_id: str,
) -> dict[str, str | int | float]:
    """Calculate identifiers and physicochemical descriptors."""

    canonical_smiles = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=False,
    )

    isomeric_smiles = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True,
    )

    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.MolToInchiKey(mol)

    return {
        "NQL_ID": nql_id,
        "Source_ID": get_property(mol, "identifier"),
        "Source_Database": "COCONUT",
        "Library_Version": "1.0",
        "Quinone_Class": get_property(mol, "Quinone_Class"),
        "Canonical_SMILES": canonical_smiles,
        "Isomeric_SMILES": isomeric_smiles,
        "InChI": inchi,
        "InChIKey": inchikey,
        "Molecular_Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular_Weight": round(Descriptors.MolWt(mol), 4),
        "Exact_Mass": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 4),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 4),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Rotatable_Bonds": Lipinski.NumRotatableBonds(mol),
        "Heavy_Atoms": mol.GetNumHeavyAtoms(),
        "Aromatic_Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Total_Rings": rdMolDescriptors.CalcNumRings(mol),
    }


def main() -> int:
    """Curate, deduplicate, annotate, and export the master library."""

    if not INPUT_SDF.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_SDF}")

    supplier = Chem.SDMolSupplier(
        str(INPUT_SDF),
        removeHs=False,
        sanitize=True,
    )

    sdf_writer = Chem.SDWriter(str(MASTER_SDF))

    unique_records: list[dict[str, str | int | float]] = []
    duplicate_rows: list[list[str]] = []

    seen_inchikeys: dict[str, str] = {}

    total_records = 0
    valid_records = 0
    invalid_records = 0
    duplicate_records = 0

    for mol in supplier:
        total_records += 1

        if mol is None:
            invalid_records += 1
            continue

        valid_records += 1

        temporary_id = f"NQLMM{valid_records:06d}"
        record = calculate_record(mol, temporary_id)

        inchikey = str(record["InChIKey"])

        if not inchikey:
            invalid_records += 1
            continue

        if inchikey in seen_inchikeys:
            duplicate_records += 1
            duplicate_rows.append(
                [
                    get_property(mol, "identifier"),
                    inchikey,
                    seen_inchikeys[inchikey],
                ]
            )
            continue

        final_id = f"NQLMM{len(unique_records) + 1:06d}"
        record["NQL_ID"] = final_id
        seen_inchikeys[inchikey] = final_id

        for key, value in record.items():
            mol.SetProp(key, str(value))

        mol.SetProp("_Name", final_id)

        sdf_writer.write(mol)
        unique_records.append(record)

    sdf_writer.close()

    with MASTER_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(unique_records)

    with MASTER_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
        )
        writer.writeheader()
        writer.writerows(unique_records)

    with MASTER_SMI.open("w", encoding="utf-8") as handle:
        handle.write("Canonical_SMILES\tNQL_ID\n")
        for record in unique_records:
            handle.write(
                f"{record['Canonical_SMILES']}\t{record['NQL_ID']}\n"
            )

    with DUPLICATE_REPORT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Duplicate_Source_ID",
                "InChIKey",
                "Retained_NQL_ID",
            ]
        )
        writer.writerows(duplicate_rows)

    with QC_REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Input_records", total_records])
        writer.writerow(["Valid_records", valid_records])
        writer.writerow(["Invalid_records", invalid_records])
        writer.writerow(["Duplicate_records", duplicate_records])
        writer.writerow(["Final_unique_records", len(unique_records)])
        writer.writerow(["Master_SDF", MASTER_SDF])
        writer.writerow(["Master_TSV", MASTER_TSV])
        writer.writerow(["Master_CSV", MASTER_CSV])
        writer.writerow(["Master_SMILES", MASTER_SMI])

    with STATISTICS_FILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Retrieved_quinones", total_records])
        writer.writerow(["Unique_quinones", len(unique_records)])
        writer.writerow(["Duplicates_removed", duplicate_records])

    print("NQL-MM master-library curation completed.")
    print(f"Input records: {total_records}")
    print(f"Valid records: {valid_records}")
    print(f"Invalid records: {invalid_records}")
    print(f"Duplicates removed: {duplicate_records}")
    print(f"Final unique records: {len(unique_records)}")
    print(f"Master SDF: {MASTER_SDF}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


