import csv

base = "08_sensitivity_analysis/"

metrics = {
    "RMSD": "min",
    "Hbond_occupancy": "max",
    "Contact_occupancy": "max",
    "MMGBSA": "min",
    "QED": "max",
    "ADMET": "max"
}

with open(base + "ligand_metrics.csv", newline="") as f:
    rows = list(csv.DictReader(f))

# Rank each metric; tied values receive average ranks
for metric, direction in metrics.items():
    values = sorted(
        [float(r[metric]) for r in rows],
        reverse=(direction == "max")
    )

    for row in rows:
        value = float(row[metric])
        positions = [
            i + 1 for i, v in enumerate(values) if v == value
        ]
        row[metric + "_rank"] = sum(positions) / len(positions)

# Calculate mean rank across six metrics
for row in rows:
    row["Mean_rank"] = sum(
        row[m + "_rank"] for m in metrics
    ) / len(metrics)

rows.sort(key=lambda r: (r["Mean_rank"], r["Ligand"]))

output = base + "03_rank_aggregation_results.csv"

with open(output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Final_rank", "Ligand", "Mean_rank"])
    for i, row in enumerate(rows, 1):
        writer.writerow([
            i, row["Ligand"], f'{row["Mean_rank"]:.6f}'
        ])

print("Rank  Ligand        Mean rank")
for i, row in enumerate(rows, 1):
    print(f'{i:>2}    {row["Ligand"]}   {row["Mean_rank"]:.4f}')

print("\nSaved:", output)
