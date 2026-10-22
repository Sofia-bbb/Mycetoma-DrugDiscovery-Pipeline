import csv

base = "08_sensitivity_analysis/"

with open(base + "ligand_metrics.csv", newline="") as f:
    rows = list(csv.DictReader(f))

metrics = {
    "RMSD": "min",
    "Hbond_occupancy": "max",
    "Contact_occupancy": "max",
    "MMGBSA": "min",
    "QED": "max",
    "ADMET": "max",
}

for metric, direction in metrics.items():
    values = [float(r[metric]) for r in rows]
    low, high = min(values), max(values)

    for row in rows:
        value = float(row[metric])
        if high == low:
            normalized = 1.0
        elif direction == "min":
            normalized = (high - value) / (high - low)
        else:
            normalized = (value - low) / (high - low)

        row[metric + "_normalized"] = normalized

for row in rows:
    row["Equal_weight_score"] = sum(
        row[m + "_normalized"] for m in metrics
    ) / len(metrics)

rows.sort(key=lambda r: r["Equal_weight_score"], reverse=True)

output = base + "01_equal_weight_results.csv"
with open(output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Rank", "Ligand", "Equal_weight_score"])
    for rank, row in enumerate(rows, 1):
        writer.writerow([
            rank,
            row["Ligand"],
            f'{row["Equal_weight_score"]:.6f}'
        ])

print("Rank  Ligand        Score")
for rank, row in enumerate(rows, 1):
    print(f'{rank:>2}    {row["Ligand"]}   {row["Equal_weight_score"]:.6f}')

print("\nSaved:", output)
