#!/usr/bin/env python3

from pathlib import Path
import csv

from rdkit import Chem
from rdkit.Chem import (
    Crippen,
    Descriptors,
    Lipinski,
    QED,
    rdMolDescriptors,
)

ROOT = Path("05_prioritization")

INPUT = ROOT / "Results" / "PLIP121_Quinone_Chemistry.csv"
OUTDIR = ROOT / "Results"
LOGDIR = ROOT / "Logs"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "PLIP121_Developability_Profile.csv"
SUMMARY = OUTDIR / "PLIP121_Developability_Summary.txt"
FAILURE = LOGDIR / "PLIP121_developability_failures.tsv"


def esol_logS(mol):

    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rotb = Lipinski.NumRotatableBonds(mol)

    aromatic_atoms = sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetIsAromatic()
    )

    heavy_atoms = mol.GetNumHeavyAtoms()

    aromatic_proportion = (
        aromatic_atoms / heavy_atoms
        if heavy_atoms > 0
        else 0
    )

    return (
        0.16
        - 1.5 * logp
        - 0.01 * (mw - 40)
        + 0.066 * rotb
        + 0.066 * aromatic_proportion
    )


def solubility_class(logS):

    if logS >= -1:
        return "High"

    if logS >= -3:
        return "Moderate"

    if logS >= -5:
        return "Low"

    return "Very_Low"


with INPUT.open(
    newline="",
    encoding="utf-8-sig"
) as f:

    rows = list(csv.DictReader(f))


print("Input compounds:", len(rows))

results = []
failures = []


for row in rows:

    ligand = row["Ligand"].strip()
    smiles = row["Canonical_SMILES"].strip()
    ctype = row["Compound_Type"].strip()

    try:

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            raise ValueError("Invalid SMILES")

        mol = Chem.RemoveHs(mol)

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)

        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rotb = Lipinski.NumRotatableBonds(mol)

        fraction_csp3 = (
            rdMolDescriptors.CalcFractionCSP3(mol)
        )

        qed = QED.qed(mol)

        logs = esol_logS(mol)

        favorable_physchem = int(
            mw <= 500
            and logp <= 5
            and tpsa <= 140
            and hbd <= 5
            and hba <= 10
            and rotb <= 10
        )

        result = dict(row)

        result.update({
            "Developability_Molecular_Weight":
                round(mw, 3),

            "Developability_LogP":
                round(logp, 3),

            "Developability_TPSA":
                round(tpsa, 3),

            "Developability_HBD":
                hbd,

            "Developability_HBA":
                hba,

            "Developability_Rotatable_Bonds":
                rotb,

            "Fraction_CSP3":
                round(fraction_csp3, 3),

            "QED":
                round(qed, 3),

            "ESOL_LogS":
                round(logs, 3),

            "ESOL_Solubility_Class":
                solubility_class(logs),

            "Favorable_Physicochemical_Profile":
                favorable_physchem,

            # Keep all references regardless of result.
            # Candidates are still not eliminated here.
            "Retain_For_ADMET":
                1,
        })

        results.append(result)

    except Exception as exc:

        failures.append(
            (ligand, str(exc), smiles)
        )


if not results:
    raise RuntimeError(
        "No compounds successfully processed."
    )


fields = list(results[0].keys())

with OUTPUT.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(results)


with FAILURE.open(
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "Compound\tReason\tStructure\n"
    )

    for item in failures:
        f.write(
            "\t".join(item) + "\n"
        )


candidates = [
    r for r in results
    if r["Compound_Type"] == "Candidate"
]

references = [
    r for r in results
    if r["Compound_Type"] == "Reference"
]


def count_class(group, cls):

    return sum(
        r["ESOL_Solubility_Class"] == cls
        for r in group
    )


summary = [
    "===== PLIP121 DEVELOPABILITY =====",
    f"Total input: {len(rows)}",
    f"Candidates processed: {len(candidates)}",
    f"References processed: {len(references)}",
    f"Failures: {len(failures)}",
    "",
    "CANDIDATE ESOL SOLUBILITY",
    f"High: {count_class(candidates, 'High')}",
    f"Moderate: {count_class(candidates, 'Moderate')}",
    f"Low: {count_class(candidates, 'Low')}",
    f"Very low: {count_class(candidates, 'Very_Low')}",
    "",
    (
        "Candidates with favorable physicochemical profile: "
        f"{sum(int(r['Favorable_Physicochemical_Profile']) for r in candidates)}"
    ),
    "",
    "REFERENCES",
]


for r in references:

    summary.append(
        f'{r["Ligand"]}: '
        f'QED={r["QED"]}; '
        f'ESOL={r["ESOL_LogS"]}; '
        f'Solubility={r["ESOL_Solubility_Class"]}; '
        f'FavorablePhyschem={r["Favorable_Physicochemical_Profile"]}'
    )


SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("\n".join(summary))
print()
print("Saved:", OUTPUT)


