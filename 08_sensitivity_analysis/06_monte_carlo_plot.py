import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np

base = "08_sensitivity_analysis/"

df = pd.read_csv(base + "04_random_weights_results.csv")

selected = {
    "NQLMM002852", "NQLMM002110", "NQLMM000191",
    "NQLMM003691", "NQLMM002327"
}

df = df.sort_values(
    "Selection_frequency_percent", ascending=False
).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(10, 6.5))
y = np.arange(len(df))
values = df["Selection_frequency_percent"].to_numpy()

# Gentle background shading for the two groups
ax.axhspan(-0.5, 4.5, color="#EFF7FD", zorder=0)
ax.axhspan(4.5, 8.5, color="#F7F7F7", zorder=0)

# Gradient-filled horizontal bars
for i, (ligand, value) in enumerate(
    zip(df["Ligand"], values)
):
    height = 0.58

    if ligand in selected:
        start, end = "#29A8F2", "#073C94"
    else:
        start, end = "#C7CDD3", "#8A939C"

    gradient = np.linspace(
        mcolors.to_rgb(start),
        mcolors.to_rgb(end),
        256
    ).reshape(1, 256, 3)

    if value >= 0.2:
        ax.imshow(
            gradient,
            extent=[0, value, i-height/2, i+height/2],
            aspect="auto",
            interpolation="bicubic",
            zorder=3
        )
    else:
        ax.barh(i, max(value, 0.08), height=height,
                color=start, zorder=3)

    ax.text(
        min(value + 1.3, 103), i,
        f"{value:.2f}%",
        va="center", fontsize=10,
        fontweight="medium"
    )

ax.set_yticks(y)
ax.set_yticklabels(df["Ligand"], fontsize=10)
ax.invert_yaxis()

# Separate prioritized and lower-priority candidates
ax.axhline(4.5, color="#8A939C", linestyle="--", linewidth=1.2)

ax.set_xlim(0, 110)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("Selection frequency (%)", fontsize=11)

ax.grid(axis="x", linestyle="--", alpha=0.25)
ax.set_axisbelow(True)

for side in ["top", "right", "left"]:
    ax.spines[side].set_visible(False)

ax.tick_params(axis="y", length=0)

# Group names placed above their respective sections
ax.text(
    0.01, 1.025, "Top five prioritized candidates",
    transform=ax.transAxes, fontsize=10,
    color="#073C94", fontweight="bold"
)

plt.tight_layout(rect=[0, 0, 1, 0.97])

fig.savefig(
    base + "06_monte_carlo_plot.png",
    dpi=400, bbox_inches="tight"
)
fig.savefig(
    base + "06_monte_carlo_plot.pdf",
    bbox_inches="tight"
)

plt.close(fig)
print("Improved Panel B created successfully.")
