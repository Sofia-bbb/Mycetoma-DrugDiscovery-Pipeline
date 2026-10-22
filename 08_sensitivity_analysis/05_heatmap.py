import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = "08_sensitivity_analysis/"

df = pd.read_csv(base + "ligand_metrics.csv").set_index("Ligand")

metrics = [
    "RMSD",
    "Hbond_occupancy",
    "Contact_occupancy",
    "MMGBSA",
    "QED",
    "ADMET"
]

selected = [
    "NQLMM000191",
    "NQLMM003691",
    "NQLMM002852",
    "NQLMM002110",
    "NQLMM002327"
]

remaining = [x for x in df.index if x not in selected]
df = df.loc[selected + remaining]

# Min-Max normalization
normalized = pd.DataFrame(index=df.index)

for metric in metrics:
    values = df[metric].astype(float)
    low, high = values.min(), values.max()

    if high == low:
        normalized[metric] = 1.0
    elif metric in ["RMSD", "MMGBSA"]:
        normalized[metric] = (high - values) / (high - low)
    else:
        normalized[metric] = (values - low) / (high - low)

fig, ax = plt.subplots(figsize=(10.5, 6.2))

im = ax.imshow(
    normalized.values,
    cmap="YlGnBu",
    vmin=0,
    vmax=1,
    aspect="auto"
)

labels = [
    "Ligand RMSD",
    "H-bond occupancy",
    "Contact occupancy",
    "MM/GBSA",
    "QED",
    "ADMET score"
]

ax.set_xticks(range(6))
ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=10)

ax.set_yticks(range(9))
ax.set_yticklabels(df.index, fontsize=10)

# Separate final five from remaining four
ax.axhline(4.5, color="#333333", linewidth=1.8)

# Annotate normalized values
for i in range(9):
    for j in range(6):
        value = normalized.iloc[i, j]
        text_color = "white" if value > 0.60 else "#222222"
        ax.text(
            j, i, f"{value:.2f}",
            ha="center", va="center",
            fontsize=9, color=text_color
        )

# Color scale
cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.035)
cbar.set_label("Normalized desirability (0–1)", fontsize=10)

# Clean scientific layout
ax.tick_params(axis="both", length=0)
plt.tight_layout()

fig.savefig(
    base + "05_normalized_heatmap.png",
    dpi=400, bbox_inches="tight"
)

fig.savefig(
    base + "05_normalized_heatmap.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("Improved heatmap PNG and PDF created successfully.")
