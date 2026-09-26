#!/usr/bin/env python3

from pathlib import Path
import subprocess
import shutil
import math
import sys

ROOT = Path(__file__).resolve().parents[2]

OLD = ROOT / "06_md_simulation/input/ligands"
NEW = ROOT / "06_md_simulation/systems"

RECEPTOR = (
    ROOT / "02_target_preparation/receptor/MmSCD_dockingframe_noH.pdb"
)

SYSTEMS = [
    "NQLMM002110",
    "NQLMM007485",
    "NQLMM002852",
    "NQLMM000788",
    "NQLMM003435",
    "NQLMM002327",
    "NQLMM003691",
    "NQLMM005370",
    "NQLMM000191",
    "Carpropamid",
    "Fenoxanil",
]

ACTIVE = {28,48,83,109,128,130}

if not RECEPTOR.exists():
    sys.exit(f"ERROR: receptor missing: {RECEPTOR}")

# ------------------------------------------------------------
# Distance QC
# ------------------------------------------------------------

def xyz(line):
    return (
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
    )


def active_site_distance(pdb):

    active = []
    ligand = []

    with pdb.open() as fh:

        for line in fh:

            if not line.startswith(("ATOM","HETATM")):
                continue

            try:
                resid = int(line[22:26])
            except:
                continue

            resname = line[17:20].strip()

            if resid in ACTIVE:
                active.append(xyz(line))

            if resname == "MOL":
                ligand.append(xyz(line))

    if not active:
        raise RuntimeError("No active-site atoms identified")

    if not ligand:
        raise RuntimeError("No MOL atoms identified")

    d = min(
        math.dist(a,b)
        for a in active
        for b in ligand
    )

    return d, len(active), len(ligand)


# ------------------------------------------------------------
# Build systems
# ------------------------------------------------------------

summary = []

for system in SYSTEMS:

    print()
    print("=" * 70)
    print("SYSTEM:", system)
    print("=" * 70)

    old = OLD / system
    new = NEW / system

    ligand_dir = new / "Ligand"
    protein_dir = new / "Protein"
    topology_dir = new / "Topology"

    ligand_dir.mkdir(parents=True, exist_ok=True)
    protein_dir.mkdir(parents=True, exist_ok=True)
    topology_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Copy receptor in correct docking coordinate frame
    # --------------------------------------------------------

    receptor_out = protein_dir / "MmSCD_dockingframe_noH.pdb"

    shutil.copy2(RECEPTOR, receptor_out)

    # --------------------------------------------------------
    # Existing ligand parameters
    # --------------------------------------------------------

    mol2_src = old / f"{system}_gaff2_neutral.mol2"
    frc_src  = old / f"{system}_correct.frcmod"

    if not mol2_src.exists():
        print("FAIL: missing MOL2")
        summary.append((system,"FAIL_MOL2","",""))
        continue

    if not frc_src.exists():
        print("FAIL: missing FRCMOD")
        summary.append((system,"FAIL_FRCMOD","",""))
        continue

    mol2_dst = ligand_dir / mol2_src.name
    frc_dst  = ligand_dir / frc_src.name

    shutil.copy2(mol2_src, mol2_dst)
    shutil.copy2(frc_src, frc_dst)

    # --------------------------------------------------------
    # LEaP input
    # --------------------------------------------------------

    leap = topology_dir / "build_corrected.leap"

    leap.write_text(f"""
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

LIG = loadmol2 ../Ligand/{system}_gaff2_neutral.mol2
loadamberparams ../Ligand/{system}_correct.frcmod

PRO = loadpdb ../Protein/MmSCD_dockingframe_noH.pdb

COM = combine {{PRO LIG}}

check LIG
check PRO
check COM

savepdb COM {system}_corrected_unsolvated.pdb
saveamberparm COM {system}_corrected_unsolvated.prmtop {system}_corrected_unsolvated.inpcrd

solvatebox COM TIP3PBOX 12.0

addions COM Na+ 0
addions COM Cl- 0

savepdb COM {system}_corrected_solvated.pdb
saveamberparm COM {system}_corrected_solvated.prmtop {system}_corrected_solvated.inpcrd

quit
""")

    # --------------------------------------------------------
    # Run tleap
    # --------------------------------------------------------

    logfile = topology_dir / "leap_corrected.log"

    with logfile.open("w") as log:

        result = subprocess.run(
            ["tleap","-f","build_corrected.leap"],
            cwd=topology_dir,
            stdout=log,
            stderr=subprocess.STDOUT
        )

    if result.returncode != 0:

        print("FAIL: tleap")
        summary.append((system,"FAIL_TLEAP","",""))
        continue

    unsolv = topology_dir / f"{system}_corrected_unsolvated.pdb"
    solv   = topology_dir / f"{system}_corrected_solvated.pdb"

    if not unsolv.exists() or not solv.exists():

        print("FAIL: output PDB missing")
        summary.append((system,"FAIL_OUTPUT","",""))
        continue

    # --------------------------------------------------------
    # Geometry QC
    # --------------------------------------------------------

    try:

        d1, a1, l1 = active_site_distance(unsolv)
        d2, a2, l2 = active_site_distance(solv)

    except Exception as e:

        print("FAIL QC:",e)
        summary.append((system,"FAIL_QC","",""))
        continue

    status = (
        "PASS"
        if d1 <= 5.0 and d2 <= 5.0
        else "FAIL_DISTANCE"
    )

    print(f"Unsolvated distance : {d1:.2f} Å")
    print(f"Solvated distance   : {d2:.2f} Å")
    print("Active-site atoms   :",a2)
    print("Ligand atoms        :",l2)
    print("QC STATUS           :",status)

    summary.append(
        (
            system,
            status,
            f"{d1:.3f}",
            f"{d2:.3f}"
        )
    )


# ------------------------------------------------------------
# Write summary
# ------------------------------------------------------------

summary_file = NEW / "Corrected_MD11_PreMD_QC.tsv"

with summary_file.open("w") as f:

    f.write(
        "System\tQC_Status\t"
        "Unsolvated_Min_Distance_A\t"
        "Solvated_Min_Distance_A\n"
    )

    for row in summary:
        f.write("\t".join(row) + "\n")


print()
print("=" * 70)
print("FINAL CORRECTED MD11 QC")
print("=" * 70)

for row in summary:

    print(
        row[0],
        "Status =",row[1],
        "Unsolvated =",row[2],
        "Solvated =",row[3]
    )

passed = sum(r[1] == "PASS" for r in summary)

print()
print("PASS:",passed,"/",len(SYSTEMS))
print("Summary:",summary_file)

if passed != len(SYSTEMS):

    print()
    print("WARNING: NOT ALL SYSTEMS PASSED.")
    print("DO NOT START MD.")

    sys.exit(1)

print()
print("ALL 11 SYSTEMS PASSED PRE-MD GEOMETRY QC.")


