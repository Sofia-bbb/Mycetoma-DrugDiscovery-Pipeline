from pathlib import Path
import subprocess

ROOT = Path("04_interaction_profiling")
LIGDIR = ROOT / "Ligands"
PDBDIR = ROOT / "Ligands_PDB"
COMPLEXDIR = ROOT / "Complexes"

RECEPTOR = Path("02_target_preparation/receptor/MmSCD_AF_H.pdb")

PDBDIR.mkdir(parents=True, exist_ok=True)
COMPLEXDIR.mkdir(parents=True, exist_ok=True)

ligands = sorted(LIGDIR.glob("*_model1.pdbqt"))

success = 0
failed = []

# Read protein atoms only once
protein_lines = []
with open(RECEPTOR) as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            protein_lines.append(line)

for src in ligands:

    name = src.name.replace("_model1.pdbqt", "")

    ligand_pdb = PDBDIR / f"{name}_model1.pdb"
    complex_pdb = COMPLEXDIR / f"{name}_complex.pdb"

    try:
        # Proper ligand conversion
        subprocess.run(
            [
                "obabel",
                str(src),
                "-O",
                str(ligand_pdb)
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        ligand_lines = []

        with open(ligand_pdb) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):

                    # make ligand explicit HETATM
                    line = "HETATM" + line[6:]

                    # assign chain Z
                    if len(line) >= 22:
                        line = line[:21] + "Z" + line[22:]

                    ligand_lines.append(line)

        if not ligand_lines:
            raise RuntimeError("No ligand atoms after conversion")

        with open(complex_pdb, "w") as f:
            f.writelines(protein_lines)
            f.write("TER\n")
            f.writelines(ligand_lines)
            f.write("END\n")

        success += 1

    except Exception as e:
        failed.append((name, str(e)))

print("===== COMPLEX BUILD =====")
print("Expected:", len(ligands))
print("Built:", success)
print("Failed:", len(failed))

if failed:
    print("\nFailures:")
    for x in failed:
        print(*x)


