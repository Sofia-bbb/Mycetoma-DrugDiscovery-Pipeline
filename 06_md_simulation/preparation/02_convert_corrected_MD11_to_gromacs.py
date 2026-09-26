#!/usr/bin/env python3

from pathlib import Path
import parmed as pmd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "06_md_simulation/systems"

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

passed = 0
failed = 0

print("===== CORRECTED MD11 AMBER -> GROMACS =====")

for system in SYSTEMS:

    print("\n=====", system, "=====")

    try:

        topdir = BASE / system / "Topology"
        outdir = BASE / system / "Solvation"

        outdir.mkdir(parents=True, exist_ok=True)

        prmtop = topdir / f"{system}_corrected_solvated.prmtop"
        inpcrd = topdir / f"{system}_corrected_solvated.inpcrd"

        gro = outdir / f"{system}_solvated.gro"
        top = outdir / "topol.top"

        if not prmtop.exists():
            raise FileNotFoundError(prmtop)

        if not inpcrd.exists():
            raise FileNotFoundError(inpcrd)

        structure = pmd.load_file(
            str(prmtop),
            xyz=str(inpcrd)
        )

        charge = sum(atom.charge for atom in structure.atoms)

        print("Atoms:", len(structure.atoms))
        print("Residues:", len(structure.residues))
        print("Total charge:", charge)

        structure.save(
            str(gro),
            format="gro",
            overwrite=True
        )

        structure.save(
            str(top),
            format="gromacs",
            overwrite=True
        )

        if not gro.exists() or gro.stat().st_size == 0:
            raise RuntimeError("GRO missing")

        if not top.exists() or top.stat().st_size == 0:
            raise RuntimeError("Topology missing")

        print("PASS")
        passed += 1

    except Exception as exc:

        print("FAILED:", exc)
        failed += 1

print("\n===== CONVERSION SUMMARY =====")
print("Passed:", passed)
print("Failed:", failed)

if failed:
    raise SystemExit(1)



