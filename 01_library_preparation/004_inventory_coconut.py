#!/usr/bin/env python3

"""
Project: Mycetoma Drug Discovery
Module: NQL-MM
Purpose:
Inventory the COCONUT database before retrieval.

Outputs:
1. Database_Inventory.tsv
2. Property_List.tsv
"""

from rdkit import Chem
from collections import Counter
import csv
from pathlib import Path

INPUT = Path(
    str(Path(__file__).resolve().parents[1] / "data" / "raw" / "coconut_sdf_2d_lite-08-2026.sdf")
)

OUTDIR = Path(__file__).resolve().parents[1] / "01_library_preparation" / "validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

inventory = OUTDIR / "Database_Inventory.tsv"
properties = OUTDIR / "Property_List.tsv"

supplier = Chem.SDMolSupplier(str(INPUT), removeHs=False)

total = 0
valid = 0
invalid = 0

prop_counter = Counter()

for mol in supplier:

    total += 1

    if mol is None:
        invalid += 1
        continue

    valid += 1

    for p in mol.GetPropNames():
        prop_counter[p] += 1

with open(inventory, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["Metric", "Value"])

    writer.writerow(["Total_records", total])
    writer.writerow(["Valid_records", valid])
    writer.writerow(["Invalid_records", invalid])
    writer.writerow(["Unique_properties", len(prop_counter)])

with open(properties, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["Property", "Occurrences"])

    for p, c in sorted(prop_counter.items()):
        writer.writerow([p, c])

print("Inventory complete")
print("Total:", total)
print("Valid:", valid)
print("Invalid:", invalid)
print("Unique properties:", len(prop_counter))


