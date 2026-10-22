import csv

base = "08_sensitivity_analysis/"
metrics = {
    "RMSD": "min",
    "Hbond_occupancy": "max",
    "Contact_occupancy": "max",
    "MMGBSA": "min",
    "QED": "max",
    "ADMET": "max",
}

original_top5 = {
    "NQLMM000191", "NQLMM003691", "NQLMM002852",
    "NQLMM002110", "NQLMM002327"
}

with open(base + "ligand_metrics.csv", newline="") as f:
    rows = list(csv.DictReader(f))

normalized = {}
for metric, direction in metrics.items():
    values = [float(r[metric]) for r in rows]
    low, high = min(values), max(values)
    normalized[metric] = {
        r["Ligand"]: (
            1.0 if high == low else
            (high - float(r[metric])) / (high - low)
            if direction == "min" else
            (float(r[metric]) - low) / (high - low)
        )
        for r in rows
    }

results = []
for excluded in metrics:
    included = [m for m in metrics if m != excluded]
    ranked = sorted(
        rows,
        key=lambda r: (
            -sum(normalized[m][r["Ligand"]] for m in included),
            r["Ligand"]
        )
    )
    top5 = [r["Ligand"] for r in ranked[:5]]
    overlap = len(set(top5) & original_top5)
    results.append([excluded, overlap, ";".join(top5)])
    print(f"Excluded: {excluded:20s} | Overlap: {overlap}/5")
    print("Top 5:", ", ".join(top5))

with open(base + "02_leave_one_out_results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Excluded_metric", "Original_top5_overlap", "Top5"])
    writer.writerows(results)

print("\nSaved:", base + "02_leave_one_out_results.csv")
