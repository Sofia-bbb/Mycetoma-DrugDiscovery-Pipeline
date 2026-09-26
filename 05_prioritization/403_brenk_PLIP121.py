#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import Counter

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

ROOT = Path("05_prioritization")

INPUT = ROOT / "Results" / "PLIP121_PAINS.csv"

OUTDIR = ROOT / "Results"
LOGDIR = ROOT / "Logs"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "PLIP121_Brenk.csv"
SUMMARY = OUTDIR / "PLIP121_Brenk_Summary.txt"
FAILURES = LOGDIR / "PLIP121_brenk_failures.tsv"

# ---------------------------------------------------------
# Brenk structural-alert catalog
# ---------------------------------------------------------

params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog(params)


def get_alerts(mol):
    alerts = []

    for match in catalog.GetMatches(mol):
        desc = match.GetDescription().strip()

        if desc and desc not in alerts:
            alerts.append(desc)

    return sorted(alerts)


# ---------------------------------------------------------
# Load compounds
# ---------------------------------------------------------

with INPUT.open(encoding="utf-8-sig") as f:
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
        "Brenk_Count": len(alerts),
        "Brenk_Alerts": ";".join(alerts),
        "Brenk_No_Alert": int(len(alerts) == 0),
        # Annotation only at this stage — no automatic exclusion.
        "Brenk_Retain_For_Analysis": 1,
    })

    results.append(result)


# ---------------------------------------------------------
# Save full results
# ---------------------------------------------------------

if not results:
    raise RuntimeError("No compounds were successfully processed.")

fields = list(results[0].keys())

with OUTPUT.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)


with FAILURES.open("w", encoding="utf-8") as f:
    f.write("Ligand\tReason\tValue\n")

    for item in failures:
        f.write("\t".join(item) + "\n")


# ---------------------------------------------------------
# Separate candidates and references
# ---------------------------------------------------------

candidates = [
    r for r in results
    if r.get("Compound_Type", "") == "Candidate"
]

references = [
    r for r in results
    if r.get("Compound_Type", "") == "Reference"
]


def no_alert_count(group):
    return sum(
        int(r["Brenk_No_Alert"]) == 1
        for r in group
    )


# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

summary = [
    "===== PLIP121 BRENK ANALYSIS =====",
    f"Total input: {len(rows)}",
    f"Successfully processed: {len(results)}",
    f"Failures: {len(failures)}",
    "",
    "CANDIDATES",
    f"Total candidates: {len(candidates)}",
    f"No Brenk alert: {no_alert_count(candidates)}",
    (
        "One or more Brenk alerts: "
        f"{len(candidates) - no_alert_count(candidates)}"
    ),
    "",
    "REFERENCES",
]

for r in references:
    summary.append(
        f'{r["Ligand"]}: '
        f'Brenk_Count={r["Brenk_Count"]}; '
        f'Alerts={r["Brenk_Alerts"] or "NONE"}'
    )

summary.append("")
summary.append("Most frequent Brenk alerts:")

for name, count in alert_counter.most_common(25):
    summary.append(f"{name}: {count}")

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("\n".join(summary))
print()
print("Saved:", OUTPUT)


