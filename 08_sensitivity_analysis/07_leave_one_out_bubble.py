import pandas as pd
import matplotlib.pyplot as plt

# input
csv_file = "08_sensitivity_analysis/02_leave_one_out_results.csv"
df = pd.read_csv(csv_file)

# desired ligand order
ligands = [
    "NQLMM002852",
    "NQLMM002110",
    "NQLMM000191",
    "NQLMM003691",
    "NQLMM002327",
    "NQLMM000788",
    "NQLMM007485",
    "NQLMM005370",
    "NQLMM003435",
]

# row labels
metric_labels = {
    "RMSD": "Exclude RMSD",
    "Hbond_occupancy": "Exclude H-bond occupancy",
    "Contact_occupancy": "Exclude contact occupancy",
    "MMGBSA": "Exclude MM/GBSA",
    "QED": "Exclude QED",
    "ADMET": "Exclude ADMET",
}

# build plotting table
rows = []
for _, row in df.iterrows():
    excluded = row["Excluded_metric"]
    overlap = int(row["Original_top5_overlap"])
    top5 = str(row["Top5"]).split(";")
    top5_set = set(top5)

    for lig in ligands:
        rows.append({
            "Excluded_metric": excluded,
            "Label": metric_labels.get(excluded, excluded),
            "Ligand": lig,
            "Selected": 1 if lig in top5_set else 0,
            "Overlap": overlap,
        })

plot_df = pd.DataFrame(rows)

# order rows top to bottom
row_order = [
    "RMSD",
    "Hbond_occupancy",
    "Contact_occupancy",
    "MMGBSA",
    "QED",
    "ADMET",
]
row_positions = {m: i for i, m in enumerate(row_order)}
col_positions = {lig: i for i, lig in enumerate(ligands)}

plot_df["x"] = plot_df["Ligand"].map(col_positions)
plot_df["y"] = plot_df["Excluded_metric"].map(row_positions)

fig, ax = plt.subplots(figsize=(12, 6.5))

# colors and sizes
selected_color = "#2c7fb8"
not_selected_color = "#d9d9d9"

for _, r in plot_df.iterrows():
    ax.scatter(
        r["x"],
        r["y"],
        s=420 if r["Selected"] == 1 else 90,
        c=selected_color if r["Selected"] == 1 else not_selected_color,
        edgecolors="black" if r["Selected"] == 1 else "#888888",
        linewidths=0.8,
        zorder=3
    )

# axes formatting
ax.set_xticks(range(len(ligands)))
ax.set_xticklabels(ligands, rotation=45, ha="right", fontsize=10)

ax.set_yticks(range(len(row_order)))
ax.set_yticklabels([metric_labels[m] for m in row_order], fontsize=10)
ax.invert_yaxis()

# grid
ax.set_axisbelow(True)
ax.grid(True, linestyle="--", alpha=0.3)

# title
ax.set_title("Leave-One-Metric-Out Sensitivity Analysis", fontsize=15, fontweight="bold", pad=14)

# extend x-limit for overlap annotation
ax.set_xlim(-0.5, len(ligands) + 1.2)

# add overlap text on right
overlap_map = dict(zip(df["Excluded_metric"], df["Original_top5_overlap"]))
for metric in row_order:
    y = row_positions[metric]
    overlap = overlap_map[metric]
    ax.text(
        len(ligands) + 0.15,
        y,
        f"{overlap}/5",
        va="center",
        ha="left",
        fontsize=10,
        fontweight="bold"
    )

ax.text(
    len(ligands) + 0.15,
    -0.75,
    "Top-5 overlap",
    va="center",
    ha="left",
    fontsize=10,
    fontweight="bold"
)

# legend-like guide
ax.scatter([], [], s=420, c=selected_color, edgecolors="black", linewidths=0.8, label="Selected in top 5")
ax.scatter([], [], s=90, c=not_selected_color, edgecolors="#888888", linewidths=0.8, label="Not selected")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)

plt.tight_layout()

out_png = "08_sensitivity_analysis/07_leave_one_out_bubble.png"
out_pdf = "08_sensitivity_analysis/07_leave_one_out_bubble.pdf"

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
print("Panel C created successfully.")
