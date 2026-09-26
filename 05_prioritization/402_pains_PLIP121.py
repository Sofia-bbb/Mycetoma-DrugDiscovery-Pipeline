#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import Counter

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


ROOT = Path("05_prioritization")

INPUT = ROOT / "Results" / "PLIP121_Drug_Likeness_Progression.csv"
OUTDIR = ROOT / "Results"
LOGDIR = ROOT / "Logs"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "PLIP121_PAINS.csv"
SUMMARY = OUTDIR / "PLIP121_PAINS_Summary.txt"
FAILURES = LOGDIR / "PLIP121_pains_failures.tsv"


params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)

catalog = FilterCatalog(params)


def get_alerts(mol):
    alerts = []

    for match in catalog.GetMatches(mol):
        desc = match.GetDescription().strip()

        if desc and desc not in alerts:
            alerts.append(desc)

    return sorted(alerts)


with open(INPUT, encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

print("Input compounds:", len(rows))

results = []
failures = []
alert_counter = Counter()

for row in rows:

    ligand = row["Ligand"].strip()
    smiles = row["Canonical_SMILES"].strip()

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        failures.append((ligand, "Invalid SMILES", smiles))
        continue

    alerts = get_alerts(mol)

    for alert in alerts:
        alert_counter[alert] += 1

    result = dict(row)

    result.update({
        "PAINS_Count": len(alerts),
        "PAINS_Alerts": ";".join(alerts),
        "PAINS_No_Alert": int(len(alerts) == 0),
        "Retain_For_Further_Analysis": 1,
    })

    results.append(result)


fields = list(results[0].keys())

with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)


with open(FAILURES, "w") as f:
    f.write("Ligand\tReason\tValue\n")

    for x in failures:
        f.write("\t".join(x) + "\n")


candidates = [
    r for r in results
    if r["Compound_Type"] == "Candidate"
]

references = [
    r for r in results
    if r["Compound_Type"] == "Reference"
]


def no_alert_count(group):
    return sum(
        int(r["PAINS_No_Alert"]) == 1
        for r in group
    )


summary = [
    "===== PLIP121 PAINS ANALYSIS =====",
    f"Total input: {len(rows)}",
    f"Successfully processed: {len(results)}",
    f"Failures: {len(failures)}",
    "",
    "CANDIDATES",
    f"No PAINS alert: {no_alert_count(candidates)}",
    f"One or more PAINS alerts: {len(candidates) - no_alert_count(candidates)}",
    "",
    "REFERENCES",
]

for r in references:
    summary.append(
        f'{r["Ligand"]}: '
        f'PAINS_Count={r["PAINS_Count"]}; '
        f'Alerts={r["PAINS_Alerts"] or "NONE"}'
    )

summary.append("")
summary.append("Most frequent PAINS alerts:")

for name, count in alert_counter.most_common(20):
    summary.append(f"{name}: {count}")


SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("\n".join(summary))
print()
print("Saved:", OUTPUT)


