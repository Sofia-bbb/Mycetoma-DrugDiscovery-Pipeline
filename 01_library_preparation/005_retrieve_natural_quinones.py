#!/usr/bin/env python3

"""
Project: Mycetoma Drug Discovery
Module: NQL-MM

Purpose
-------
Retrieve natural quinones from the COCONUT database using SMARTS patterns.

Outputs
-------
Natural_Quinones.sdf
Natural_Quinones.tsv
Retrieval_Summary.tsv
"""

from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger
import csv

RDLogger.DisableLog("rdApp.*")

INPUT = Path(
    str(Path(__file__).resolve().parents[1] / "data" / "raw" / "coconut_sdf_2d_lite-08-2026.sdf")
)

OUTDIR = Path(__file__).resolve().parents[1] / "data" / "processed"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTSDF = OUTDIR / "Natural_Quinones.sdf"
OUTTSV = OUTDIR / "Natural_Quinones.tsv"
OUTSUM = OUTDIR / "Retrieval_Summary.tsv"

patterns = {

    "1,4-Naphthoquinone":
    Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12"),

    "1,2-Naphthoquinone":
    Chem.MolFromSmarts("O=C1C(=O)C=Cc2ccccc12"),

    "Anthraquinone":
    Chem.MolFromSmarts("O=C1c2ccccc2C(=O)c2ccccc12")

}

supplier = Chem.SDMolSupplier(str(INPUT), removeHs=False)

writer = Chem.SDWriter(str(OUTSDF))

table = open(OUTTSV, "w", newline="")
csvwriter = csv.writer(table, delimiter="\t")

csvwriter.writerow([
    "Identifier",
    "Subclass"
])

counts = {}

total = 0
hits = 0

for mol in supplier:

    total += 1

    if mol is None:
        continue

    ident = ""

    if mol.HasProp("identifier"):
        ident = mol.GetProp("identifier")

    subclasses = []

    for name, patt in patterns.items():

        if mol.HasSubstructMatch(patt):

            subclasses.append(name)

            counts[name] = counts.get(name, 0) + 1

    if len(subclasses) == 0:
        continue

    hits += 1

    mol.SetProp("Quinone_Class",
                ";".join(subclasses))

    writer.write(mol)

    csvwriter.writerow([
        ident,
        ";".join(subclasses)
    ])

writer.close()
table.close()

with open(OUTSUM, "w", newline="") as f:

    w = csv.writer(f, delimiter="\t")

    w.writerow(["Metric", "Value"])

    w.writerow(["Total molecules", total])

    w.writerow(["Natural quinones", hits])

    for k in counts:

        w.writerow([k, counts[k]])

print("Finished")
print("Total:", total)
print("Natural quinones:", hits)


