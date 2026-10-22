import csv
import random
from collections import Counter

base = "08_sensitivity_analysis/"
random.seed(42)

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

# Normalize all six metrics to the 0–1 range
normalized = {}

for metric, direction in metrics.items():
    values = [float(r[metric]) for r in rows]
    low, high = min(values), max(values)

    normalized[metric] = {}
    for row in rows:
        value = float(row[metric])
        if high == low:
            score = 1.0
        elif direction == "min":
            score = (high - value) / (high - low)
        else:
            score = (value - low) / (high - low)

        normalized[metric][row["Ligand"]] = score

# Monte Carlo sensitivity analysis
n_iterations = 100000
selection_counts = Counter()
exact_matches = 0

for _ in range(n_iterations):
    # Dirichlet(1,...,1): uniform over the weight simplex
    raw_weights = [random.expovariate(1.0) for _ in metrics]
    total = sum(raw_weights)
    weights = [w / total for w in raw_weights]

    scores = {}
    for row in rows:
        ligand = row["Ligand"]
        scores[ligand] = sum(
            weight * normalized[metric][ligand]
            for weight, metric in zip(weights, metrics)
        )

    top5 = set(sorted(
        scores,
        key=lambda ligand: (-scores[ligand], ligand)
    )[:5])

    selection_counts.update(top5)

    if top5 == original_top5:
        exact_matches += 1

output = base + "04_random_weights_results.csv"

with open(output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Ligand", "Selection_frequency_percent"])
    for ligand, count in selection_counts.most_common():
        writer.writerow([ligand, f"{100 * count / n_iterations:.4f}"])

print("Iterations:", n_iterations)
print("Random seed: 42")
print("Exact original top-five matches:", exact_matches)
print(f"Exact match percentage: {100 * exact_matches / n_iterations:.4f}%")

print("\nLigand selection frequencies:")
for ligand, count in selection_counts.most_common():
    print(f"{ligand}: {100 * count / n_iterations:.4f}%")

print("\nSaved:", output)
