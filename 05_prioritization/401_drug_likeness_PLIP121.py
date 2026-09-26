#!/usr/bin/env python3

from pathlib import Path
import csv
import math

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


ROOT = Path("05_prioritization")

INPUT = ROOT / "Input" / "PLIP121_Candidates_Plus_References.csv"
OUTDIR = ROOT / "Results"
LOGDIR = ROOT / "Logs"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "PLIP121_Drug_Likeness.csv"
PASSED = OUTDIR / "PLIP121_Drug_Likeness_Progression.csv"
SUMMARY = OUTDIR / "PLIP121_Drug_Likeness_Summary.txt"
FAILURES = LOGDIR / "PLIP121_drug_likeness_failures.tsv"


def safe_round(value, digits=3):
    try:
        value = float(value)
        if math.isfinite(value):
            return round(value, digits)
    except:
        pass
    return ""


def lipinski(mw, logp, hbd, hba):
    violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10,
    ])
    return int(violations == 0), violations


def veber(rot, tpsa):
    violations = sum([
        rot > 10,
        tpsa > 140,
    ])
    return int(violations == 0), violations


def ghose(mw, logp, mr, atoms):
    violations = sum([
        not (160 <= mw <= 480),
        not (-0.4 <= logp <= 5.6),
        not (40 <= mr <= 130),
        not (20 <= atoms <= 70),
    ])
    return int(violations == 0), violations


def egan(logp, tpsa):
    violations = sum([
        logp > 5.88,
        tpsa > 131.6,
    ])
    return int(violations == 0), violations


def muegge(mw, logp, tpsa, rings, carbon, hetero, rot, hbd, hba):
    violations = sum([
        not (200 <= mw <= 600),
        not (-2 <= logp <= 5),
        tpsa > 150,
        rings > 7,
        carbon <= 4,
        hetero <= 1,
        rot > 15,
        hbd > 5,
        hba > 10,
    ])
    return int(violations == 0), violations


with open(INPUT, encoding="utf-8-sig") as f:
    source_rows = list(csv.DictReader(f))

print("Input compounds:", len(source_rows))

results = []
failures = []

for row in source_rows:

    ligand = row["Ligand"].strip()
    smiles = row["Canonical_SMILES"].strip()
    ctype = row["Compound_Type"].strip()

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        failures.append((ligand, "Invalid SMILES", smiles))
        continue

    mol = Chem.RemoveHs(mol)

    mw = Descriptors.MolWt(mol)
    exact_mw = Descriptors.ExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    mr = Crippen.MolMR(mol)

    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = Lipinski.NumRotatableBonds(mol)

    tpsa = rdMolDescriptors.CalcTPSA(mol)
    csp3 = rdMolDescriptors.CalcFractionCSP3(mol)

    heavy = mol.GetNumHeavyAtoms()
    atoms = mol.GetNumAtoms()

    hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
    aliphatic = rdMolDescriptors.CalcNumAliphaticRings(mol)

    carbon = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6
    )

    charge = Chem.GetFormalCharge(mol)

    lp, lpv = lipinski(mw, logp, hbd, hba)
    vp, vpv = veber(rot, tpsa)
    gp, gpv = ghose(mw, logp, mr, atoms)
    ep, epv = egan(logp, tpsa)
    mp, mpv = muegge(
        mw, logp, tpsa,
        rings, carbon, hetero,
        rot, hbd, hba
    )

    passed_rules = lp + vp + gp + ep + mp

    candidate_pass = (
        lp == 1
        and vp == 1
        and passed_rules >= 3
    )

    retain = (
        candidate_pass
        if ctype == "Candidate"
        else True
    )

    result = dict(row)

    result.update({
        "Molecular_Weight": safe_round(mw),
        "Exact_Molecular_Weight": safe_round(exact_mw),
        "LogP": safe_round(logp),
        "Molar_Refractivity": safe_round(mr),
        "TPSA": safe_round(tpsa),
        "H_Bond_Donors": hbd,
        "H_Bond_Acceptors": hba,
        "Rotatable_Bonds": rot,
        "Heavy_Atom_Count": heavy,
        "Total_Atom_Count": atoms,
        "Carbon_Atom_Count": carbon,
        "Heteroatom_Count": hetero,
        "Ring_Count": rings,
        "Aromatic_Ring_Count": aromatic,
        "Aliphatic_Ring_Count": aliphatic,
        "Fraction_CSP3": safe_round(csp3),
        "Formal_Charge": charge,

        "Lipinski_Pass": lp,
        "Lipinski_Violations": lpv,

        "Veber_Pass": vp,
        "Veber_Violations": vpv,

        "Ghose_Pass": gp,
        "Ghose_Violations": gpv,

        "Egan_Pass": ep,
        "Egan_Violations": epv,

        "Muegge_Pass": mp,
        "Muegge_Violations": mpv,

        "DrugLikeness_Rules_Passed": passed_rules,
        "All_Five_Rules_Pass": int(passed_rules == 5),

        "Candidate_DrugLikeness_Pass": int(candidate_pass),
        "Retain_For_Progression": int(retain),
    })

    results.append(result)


progression = [
    r for r in results
    if int(r["Retain_For_Progression"]) == 1
]


fields = list(results[0].keys())

with open(OUTPUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(results)

with open(PASSED, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(progression)

with open(FAILURES, "w") as f:
    f.write("Ligand\tReason\tValue\n")
    for x in failures:
        f.write("\t".join(x) + "\n")


def count_pass(name, ctype=None):
    return sum(
        int(r[name]) == 1
        and (
            ctype is None
            or r["Compound_Type"] == ctype
        )
        for r in results
    )


candidates = [
    r for r in results
    if r["Compound_Type"] == "Candidate"
]

references = [
    r for r in results
    if r["Compound_Type"] == "Reference"
]

candidate_progression = [
    r for r in candidates
    if int(r["Candidate_DrugLikeness_Pass"]) == 1
]


summary = [
    "===== PLIP121 DRUG-LIKENESS =====",
    f"Total input: {len(source_rows)}",
    f"Candidates: {len(candidates)}",
    f"References: {len(references)}",
    f"Successfully processed: {len(results)}",
    f"Failures: {len(failures)}",
    "",
    "CANDIDATES",
    f"Lipinski pass: {count_pass('Lipinski_Pass', 'Candidate')}",
    f"Veber pass: {count_pass('Veber_Pass', 'Candidate')}",
    f"Ghose pass: {count_pass('Ghose_Pass', 'Candidate')}",
    f"Egan pass: {count_pass('Egan_Pass', 'Candidate')}",
    f"Muegge pass: {count_pass('Muegge_Pass', 'Candidate')}",
    f"All five pass: {count_pass('All_Five_Rules_Pass', 'Candidate')}",
    f"Candidate progression: {len(candidate_progression)}",
    "",
    "REFERENCES",
]

for r in references:
    summary.append(
        f'{r["Ligand"]}: '
        f'{r["DrugLikeness_Rules_Passed"]}/5 rules; '
        f'Lipinski={r["Lipinski_Pass"]}; '
        f'Veber={r["Veber_Pass"]}; '
        f'Retained=YES'
    )

summary.append("")
summary.append(f"Total retained including references: {len(progression)}")

SUMMARY.write_text(
    "\n".join(summary) + "\n"
)

print("\n".join(summary))
print()
print("Saved:", OUTPUT)
print("Progression:", PASSED)


