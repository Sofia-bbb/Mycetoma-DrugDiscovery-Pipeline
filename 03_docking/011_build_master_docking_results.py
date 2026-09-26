from pathlib import Path
import math
import csv

ROOT = Path("03_docking")
POSES = ROOT / "poses"
LIST = ROOT / "ligand_list.txt"
REC = Path("02_target_preparation/receptor/MmSCD_receptor.pdbqt")

OUT = ROOT / "Master_Corrected_Docking_Results.tsv"

targets = {
    28: "TYR28",
    48: "TYR48",
    83: "HIS83",
    109: "HIS109",
    128: "SER128",
    130: "ASN130"
}

def dist(a, b):
    return math.sqrt(
        (a[0]-b[0])**2 +
        (a[1]-b[1])**2 +
        (a[2]-b[2])**2
    )

# -----------------------------
# Read receptor coordinates
# -----------------------------

protein = []
catalytic = {}

for line in open(REC):
    if not line.startswith("ATOM"):
        continue

    try:
        rn = int(line[22:26])
        xyz = tuple(
            float(line[i:j])
            for i,j in [(30,38),(38,46),(46,54)]
        )
    except:
        continue

    protein.append(xyz)

    if rn in targets:
        catalytic.setdefault(rn, []).append(xyz)

print("Protein atoms:", len(protein))
print("Catalytic residues found:", sorted(catalytic))

# -----------------------------
# Analyze completed poses
# -----------------------------

records = []

ligands = LIST.read_text().splitlines()

for task, lig in enumerate(ligands, start=1):

    name = Path(lig).stem
    pose = POSES / f"{name}_docked.pdbqt"

    if not pose.exists() or pose.stat().st_size == 0:
        continue

    score = None
    coords = []
    inside = False

    for line in open(pose):

        if line.startswith("MODEL"):
            try:
                inside = (line.split()[1] == "1")
            except:
                inside = False
            continue

        if line.startswith("ENDMDL") and inside:
            break

        if inside and line.startswith("REMARK VINA RESULT:"):
            try:
                score = float(line.split()[3])
            except:
                pass

        if inside and line.startswith("ATOM"):
            try:
                xyz = tuple(
                    float(line[i:j])
                    for i,j in [(30,38),(38,46),(46,54)]
                )
                coords.append(xyz)
            except:
                pass

    if score is None or not coords:
        continue

    protein_min = min(
        dist(l,p)
        for l in coords
        for p in protein
    )

    catdist = {}

    for rn, atoms in catalytic.items():

        catdist[rn] = min(
            dist(l,r)
            for l in coords
            for r in atoms
        )

    active_min = min(catdist.values())

    cat5 = sum(
        d <= 5.0
        for d in catdist.values()
    )

    records.append({
        "Task": task,
        "Ligand": name,
        "Vina_Score": score,
        "Protein_Min_A": protein_min,
        "ActiveSite_Min_A": active_min,
        "Catalytic_Residues_Within_5A": cat5,
        "TYR28_A": catdist.get(28, ""),
        "TYR48_A": catdist.get(48, ""),
        "HIS83_A": catdist.get(83, ""),
        "HIS109_A": catdist.get(109, ""),
        "SER128_A": catdist.get(128, ""),
        "ASN130_A": catdist.get(130, "")
    })

# -----------------------------
# Rank by corrected docking
# -----------------------------

records.sort(key=lambda x: x["Vina_Score"])

fields = [
    "Task",
    "Ligand",
    "Vina_Score",
    "Protein_Min_A",
    "ActiveSite_Min_A",
    "Catalytic_Residues_Within_5A",
    "TYR28_A",
    "TYR48_A",
    "HIS83_A",
    "HIS109_A",
    "SER128_A",
    "ASN130_A"
]

with open(OUT, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(records)

print()
print("===== MASTER CORRECTED DOCKING =====")
print("Completed poses analyzed:", len(records))
print("Output:", OUT)

print()
print("TOP 20 CORRECTED ACTIVE-SITE DOCKING")

for r in records[:20]:

    print(
        f"{r['Ligand']:14s} "
        f"Score={r['Vina_Score']:7.3f} "
        f"ActiveSite={r['ActiveSite_Min_A']:5.2f}A "
        f"Cat<=5A={r['Catalytic_Residues_Within_5A']}"
    )



# ---------------------------------------------------------
# Select candidates for downstream interaction profiling
# Criteria:
#   Vina score <= -10.0 kcal/mol
#   >= 3 of 6 predefined binding-site residues within 5 Å
# ---------------------------------------------------------

advanced = [
    r for r in records
    if r["Vina_Score"] <= -10.0
    and r["Catalytic_Residues_Within_5A"] >= 3
]

advanced_out = ROOT / "Advanced_Candidates.tsv"

with open(advanced_out, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(advanced)

print("Advanced candidates:", len(advanced))
print("Output:", advanced_out)
