#!/usr/bin/env python3

"""
Structurally align the MmSCD AlphaFold model with chain A of the
experimental 2STD scytalone-dehydratase structure using CEAligner.
"""

from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.cealign import CEAligner


PROJECT_ROOT = Path(__file__).resolve().parents[1]

MM_PDB = (
    PROJECT_ROOT
    / "02_target_preparation"
    / "receptor"
    / "MmSCD_AF.pdb"
)

REFERENCE_PDB = (
    PROJECT_ROOT
    / "02_target_preparation"
    / "reference"
    / "2STD_carpropamid_complex.pdb"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "02_target_preparation"
    / "reference"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SUPERPOSED_FILE = PROJECT_ROOT / "02_target_preparation" / "receptor" / "MmSCD_CE_superposed_on_2STD.pdb"
RESULTS_FILE = OUTPUT_DIR / "CE_Structural_Superposition.tsv"


def get_chain(model, preferred_chain="A"):
    """Return the preferred protein chain, or the first available chain."""

    if preferred_chain in model:
        return model[preferred_chain]

    chains = list(model.get_chains())

    if not chains:
        raise RuntimeError("No chains were found in the structure.")

    return chains[0]


def main() -> None:
    """Run CE structural alignment and save the transformed MmSCD model."""

    if not MM_PDB.exists():
        raise FileNotFoundError(f"MmSCD structure not found: {MM_PDB}")

    if not REFERENCE_PDB.exists():
        raise FileNotFoundError(
            f"Reference structure not found: {REFERENCE_PDB}"
        )

    parser = PDBParser(QUIET=True)

    mm_structure = parser.get_structure("MmSCD", str(MM_PDB))
    reference_structure = parser.get_structure(
        "2STD",
        str(REFERENCE_PDB),
    )

    mm_model = mm_structure[0]
    reference_model = reference_structure[0]

    mm_chain = get_chain(mm_model, "A")
    reference_chain = get_chain(reference_model, "A")

    aligner = CEAligner()

    aligner.set_reference(reference_model)
    aligner.align(mm_model, transform=True)

    io = PDBIO()
    io.set_structure(mm_structure)
    io.save(str(SUPERPOSED_FILE))

    with RESULTS_FILE.open("w", encoding="utf-8") as handle:
        handle.write("Metric\tValue\n")
        handle.write(f"Reference_structure\t{REFERENCE_PDB.name}\n")
        handle.write("Reference_model\t0\n")
        handle.write("Moving_model\t0\n")
        handle.write(f"CE_RMSD_Angstrom\t{aligner.rms:.4f}\n")
        handle.write(f"Superposed_output\t{SUPERPOSED_FILE.name}\n")

    print("CE structural alignment completed.")
    print("Reference model: 0")
    print("Moving model: 0")
    print(f"CE RMSD: {aligner.rms:.4f} Å")
    print(f"Saved structure: {SUPERPOSED_FILE}")
    print(f"Saved results: {RESULTS_FILE}")


if __name__ == "__main__":
    main()

