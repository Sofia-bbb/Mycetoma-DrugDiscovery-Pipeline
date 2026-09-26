from pathlib import Path
import csv

ROOT = Path(".")

DOCK = ROOT / "03_docking/Advanced_Candidates.tsv"
PLIP = ROOT / "04_interaction_profiling/Analysis/Supplementary_Table_S9_PLIP_Master_Interactions.tsv"
OUT = ROOT / "04_interaction_profiling/Analysis"

OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Load docking
# ------------------------------------------------------------

docking = {}

with open(DOCK) as f:
    reader = csv.DictReader(f, delimiter="\t")
    docking_fields = reader.fieldnames

    for r in reader:
        # identify ligand ID column robustly
        lid = (
            r.get("Ligand_ID")
            or r.get("NQL_ID")
            or r.get("ID")
            or r.get("Ligand")
        )

        if lid:
            docking[lid] = r

print("Docking candidates loaded:", len(docking))
print("Docking columns:", docking_fields)

# ------------------------------------------------------------
# Load PLIP
# ------------------------------------------------------------

plip = {}

with open(PLIP) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        plip[r["Ligand_ID"]] = r

print("PLIP rows loaded:", len(plip))

# ------------------------------------------------------------
# Detect important docking columns
# ------------------------------------------------------------

def find_col(possible):
    for c in docking_fields:
        lc = c.lower()
        for p in possible:
            if p in lc:
                return c
    return None

score_col = find_col(["vina", "score", "affinity"])
protein_col = find_col(["protein_distance", "protein distance", "protein_dist"])
active_col = find_col(["active", "catalytic-site", "catalytic_site"])
coverage_col = find_col(["residues_within_5", "catalytic_residues", "within_5a"])

print()
print("Detected docking columns:")
print("Score:", score_col)
print("Protein distance:", protein_col)
print("Active-site distance:", active_col)
print("Spatial coverage:", coverage_col)

# ------------------------------------------------------------
# Integrate candidates
# ------------------------------------------------------------

rows = []

for lid, d in docking.items():

    p = plip.get(lid)

    if p is None:
        continue

    def num(x, default=0.0):
        try:
            return float(x)
        except:
            return default

    score = num(d.get(score_col, "")) if score_col else ""
    pdist = num(d.get(protein_col, "")) if protein_col else ""
    adist = num(d.get(active_col, "")) if active_col else ""
    spatial = d.get(coverage_col, "") if coverage_col else ""

    keycount = int(p["Key_Residue_Count"])

    # Strong interaction gate
    pass_gate = (
        score != ""
        and score <= -10.0
        and keycount >= 4
    )

    rows.append({
        "Ligand_ID": lid,
        "Vina_Score": score,
        "Protein_Distance_A": pdist,
        "ActiveSite_Distance_A": adist,
        "Spatial_Catalytic_Coverage": spatial,

        "PLIP_Key_Residue_Count": keycount,
        "PLIP_Key_Residues": p["Key_MmSCD_Residues"],

        "Hydrogen_Bonds": int(p["Hydrogen_Bonds"]),
        "Hydrophobic": int(p["Hydrophobic"]),
        "Water_Bridges": int(p["Water_Bridges"]),
        "Salt_Bridges": int(p["Salt_Bridges"]),
        "Pi_Stacks": int(p["Pi_Stacks"]),
        "Pi_Cation": int(p["Pi_Cation"]),
        "Halogen_Bonds": int(p["Halogen_Bonds"]),
        "Metal_Complexes": int(p["Metal_Complexes"]),
        "Total_PLIP_Interactions": int(p["Total_PLIP_Interactions"]),

        "Advance_PLIP_Gate": "YES" if pass_gate else "NO"
    })

# ------------------------------------------------------------
# Add the two references separately
# ------------------------------------------------------------

reference_scores = {
    "Carpropamid": -8.805,
    "Fenoxanil": -8.592,
}

for ref, score in reference_scores.items():

    p = plip.get(ref)

    if p is None:
        continue

    rows.append({
        "Ligand_ID": ref,
        "Vina_Score": score,
        "Protein_Distance_A": "",
        "ActiveSite_Distance_A": "",
        "Spatial_Catalytic_Coverage": "6/6",

        "PLIP_Key_Residue_Count": int(p["Key_Residue_Count"]),
        "PLIP_Key_Residues": p["Key_MmSCD_Residues"],

        "Hydrogen_Bonds": int(p["Hydrogen_Bonds"]),
        "Hydrophobic": int(p["Hydrophobic"]),
        "Water_Bridges": int(p["Water_Bridges"]),
        "Salt_Bridges": int(p["Salt_Bridges"]),
        "Pi_Stacks": int(p["Pi_Stacks"]),
        "Pi_Cation": int(p["Pi_Cation"]),
        "Halogen_Bonds": int(p["Halogen_Bonds"]),
        "Metal_Complexes": int(p["Metal_Complexes"]),
        "Total_PLIP_Interactions": int(p["Total_PLIP_Interactions"]),

        "Advance_PLIP_Gate": "REFERENCE"
    })

# ------------------------------------------------------------
# Ranking
# ------------------------------------------------------------

candidates = [r for r in rows if r["Advance_PLIP_Gate"] != "REFERENCE"]

advanced = [
    r for r in candidates
    if r["Advance_PLIP_Gate"] == "YES"
]

advanced = sorted(
    advanced,
    key=lambda r: (
        -r["PLIP_Key_Residue_Count"],
        r["Vina_Score"],
        -r["Hydrogen_Bonds"],
        -r["Total_PLIP_Interactions"],
    )
)

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

fields = list(rows[0].keys())

master = OUT / "Integrated_Docking_PLIP_Master.tsv"
with open(master, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

advout = OUT / "Integrated_Docking_PLIP_Advanced.tsv"
with open(advout, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(advanced)

print()
print("===== INTEGRATED DOCKING + PLIP =====")
print("Candidates integrated:", len(candidates))
print("Passed PLIP gate:", len(advanced))
print("Gate: Vina <= -10.0 and PLIP key residues >=4/6")

print()
print("TOP 25 INTEGRATED")

for r in advanced[:25]:
    print(
        f'{r["Ligand_ID"]:14s} '
        f'Score={r["Vina_Score"]:7.3f} '
        f'Key={r["PLIP_Key_Residue_Count"]}/6 '
        f'HB={r["Hydrogen_Bonds"]:2d} '
        f'Hyd={r["Hydrophobic"]:2d} '
        f'Total={r["Total_PLIP_Interactions"]:2d} '
        f'{r["PLIP_Key_Residues"]}'
    )

print()
print("REFERENCES")

for r in rows:
    if r["Advance_PLIP_Gate"] == "REFERENCE":
        print(
            f'{r["Ligand_ID"]:12s} '
            f'Score={r["Vina_Score"]:7.3f} '
            f'Key={r["PLIP_Key_Residue_Count"]}/6 '
            f'HB={r["Hydrogen_Bonds"]:2d} '
            f'Hyd={r["Hydrophobic"]:2d} '
            f'Total={r["Total_PLIP_Interactions"]:2d}'
        )

print()
print("Saved:")
print(master)
print(advout)


