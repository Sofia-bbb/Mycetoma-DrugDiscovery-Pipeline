#!/usr/bin/env python3

from pathlib import Path
import csv

ROOT = Path("05_prioritization")

DEVELOPABILITY = ROOT / "Results" / "PLIP121_Developability_Profile.csv"
OUTDIR = ROOT / "ADMET" / "Input"

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT = (
    OUTDIR
    / "PLIP121_ADMET_Input.csv"
)


if not DEVELOPABILITY.exists():
    raise FileNotFoundError(
        f"Missing input: {DEVELOPABILITY}"
    )


rows = []

with DEVELOPABILITY.open(
    newline="",
    encoding="utf-8-sig",
) as handle:

    reader = csv.DictReader(handle)

    required = {
        "Ligand",
        "Compound_Type",
        "Canonical_SMILES",
    }

    if not required.issubset(
        reader.fieldnames or []
    ):
        raise ValueError(
            "Missing required columns. "
            f"Found: {reader.fieldnames}"
        )

    for row in reader:

        compound = row[
            "Ligand"
        ].strip()

        smiles = row[
            "Canonical_SMILES"
        ].strip()

        source_type = row[
            "Compound_Type"
        ].strip()

        if not compound or not smiles:
            continue

        if source_type == "Candidate":
            compound_type = "Candidate"

        elif source_type == "Reference":
            compound_type = "Benchmark"

        else:
            compound_type = source_type

        rows.append(
            {
                "Compound": compound,
                "Type": compound_type,
                "smiles": smiles,
            }
        )


with OUTPUT.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "Compound",
            "Type",
            "smiles",
        ],
    )

    writer.writeheader()
    writer.writerows(rows)


candidate_count = sum(
    row["Type"] == "Candidate"
    for row in rows
)

benchmark_count = sum(
    row["Type"] == "Benchmark"
    for row in rows
)


print("=" * 60)
print("PLIP121 ADMET INPUT PREPARATION COMPLETED")
print("=" * 60)

print(
    f"Total      : {len(rows)}"
)

print(
    f"Candidates : {candidate_count}"
)

print(
    f"Benchmarks : {benchmark_count}"
)

print(
    f"Output     : {OUTPUT}"
)

print("=" * 60)


if len(rows) != 90:
    print(
        "WARNING: expected 90 compounds "
        f"but found {len(rows)}"
    )


