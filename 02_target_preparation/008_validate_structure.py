#!/usr/bin/env python3

"""
Basic structural validation of the AlphaFold MmSCD model.

Outputs:
- Structure_Validation.tsv
- Structure_Validation.md
"""

from pathlib import Path
from statistics import mean
from Bio.PDB import PDBParser
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "project.yaml"

STRUCTURE_FILE = (
    PROJECT_ROOT
    / "02_target_preparation"
    / "receptor"
    / "MmSCD_AF.pdb"
)

VALIDATION_DIR = (
    PROJECT_ROOT
    / "02_target_preparation"
    / "validation"
)

VALIDATION_DIR.mkdir(parents=True, exist_ok=True)


def load_configuration() -> dict:
    """Load project configuration."""

    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def confidence_category(plddt: float) -> str:
    """Return the AlphaFold confidence category."""

    if plddt >= 90:
        return "Very high confidence"
    if plddt >= 70:
        return "Confident"
    if plddt >= 50:
        return "Low confidence"
    return "Very low confidence"


def main() -> None:
    """Parse the PDB structure and calculate basic validation statistics."""

    if not STRUCTURE_FILE.exists():
        raise FileNotFoundError(
            f"Structure file not found: {STRUCTURE_FILE}"
        )

    config = load_configuration()

    target_name = config["target"]["protein_name"]
    accession = config["target"]["uniprot_accession"]
    gene_id = config["target"]["gene_id"]
    organism = config["organism"]["scientific_name"]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("MmSCD", STRUCTURE_FILE)

    models = list(structure.get_models())
    chains = list(structure.get_chains())
    residues = [
        residue
        for residue in structure.get_residues()
        if residue.id[0] == " "
    ]
    atoms = list(structure.get_atoms())

    chain_ids = sorted({chain.id for chain in chains})

    residue_numbers = [
        residue.id[1]
        for residue in residues
    ]

    # AlphaFold stores residue-level pLDDT values in the B-factor column.
    residue_plddt = []

    for residue in residues:
        atom_confidences = [
            atom.get_bfactor()
            for atom in residue.get_atoms()
        ]

        if atom_confidences:
            residue_plddt.append(mean(atom_confidences))

    if not residue_plddt:
        raise ValueError("No pLDDT/B-factor values were found.")

    average_plddt = mean(residue_plddt)
    minimum_plddt = min(residue_plddt)
    maximum_plddt = max(residue_plddt)

    very_high_count = sum(value >= 90 for value in residue_plddt)
    confident_count = sum(70 <= value < 90 for value in residue_plddt)
    low_count = sum(50 <= value < 70 for value in residue_plddt)
    very_low_count = sum(value < 50 for value in residue_plddt)

    low_confidence_residues = [
        residue.id[1]
        for residue, value in zip(residues, residue_plddt)
        if value < 70
    ]

    sequence_start = min(residue_numbers)
    sequence_end = max(residue_numbers)

    tsv_file = VALIDATION_DIR / "Structure_Validation.tsv"

    with tsv_file.open("w", encoding="utf-8") as handle:
        handle.write("Metric\tValue\n")
        handle.write(f"Target\t{target_name}\n")
        handle.write(f"Organism\t{organism}\n")
        handle.write(f"Gene_ID\t{gene_id}\n")
        handle.write(f"UniProt_accession\t{accession}\n")
        handle.write(f"Structure_file\t{STRUCTURE_FILE.name}\n")
        handle.write(f"Models\t{len(models)}\n")
        handle.write(f"Chains\t{','.join(chain_ids)}\n")
        handle.write(f"Residues\t{len(residues)}\n")
        handle.write(f"Residue_range\t{sequence_start}-{sequence_end}\n")
        handle.write(f"Atoms\t{len(atoms)}\n")
        handle.write(f"Average_pLDDT\t{average_plddt:.2f}\n")
        handle.write(f"Minimum_pLDDT\t{minimum_plddt:.2f}\n")
        handle.write(f"Maximum_pLDDT\t{maximum_plddt:.2f}\n")
        handle.write(f"Confidence_category\t{confidence_category(average_plddt)}\n")
        handle.write(f"Residues_pLDDT_90_or_more\t{very_high_count}\n")
        handle.write(f"Residues_pLDDT_70_to_89\t{confident_count}\n")
        handle.write(f"Residues_pLDDT_50_to_69\t{low_count}\n")
        handle.write(f"Residues_pLDDT_below_50\t{very_low_count}\n")
        handle.write(
            "Residues_below_pLDDT_70\t"
            + (
                ",".join(map(str, low_confidence_residues))
                if low_confidence_residues
                else "None"
            )
            + "\n"
        )

    markdown_file = VALIDATION_DIR / "Structure_Validation.md"

    with markdown_file.open("w", encoding="utf-8") as handle:
        handle.write("# MmSCD Structure Validation\n\n")
        handle.write(f"- **Target:** {target_name}\n")
        handle.write(f"- **Organism:** {organism}\n")
        handle.write(f"- **Gene:** {gene_id}\n")
        handle.write(f"- **UniProt:** {accession}\n")
        handle.write(f"- **Structure:** {STRUCTURE_FILE.name}\n")
        handle.write(f"- **Chains:** {', '.join(chain_ids)}\n")
        handle.write(f"- **Residues:** {len(residues)}\n")
        handle.write(f"- **Residue range:** {sequence_start}–{sequence_end}\n")
        handle.write(f"- **Atoms:** {len(atoms)}\n")
        handle.write(f"- **Average pLDDT:** {average_plddt:.2f}\n")
        handle.write(f"- **Minimum pLDDT:** {minimum_plddt:.2f}\n")
        handle.write(f"- **Maximum pLDDT:** {maximum_plddt:.2f}\n")
        handle.write(
            f"- **Overall confidence:** "
            f"{confidence_category(average_plddt)}\n\n"
        )

        handle.write("## Confidence distribution\n\n")
        handle.write(
            f"- pLDDT ≥ 90: {very_high_count} residues\n"
        )
        handle.write(
            f"- pLDDT 70–89: {confident_count} residues\n"
        )
        handle.write(
            f"- pLDDT 50–69: {low_count} residues\n"
        )
        handle.write(
            f"- pLDDT < 50: {very_low_count} residues\n\n"
        )

        handle.write("## Preliminary interpretation\n\n")

        if average_plddt >= 90:
            handle.write(
                "The AlphaFold model has very high average confidence. "
                "This supports further active-site validation, structural "
                "comparison, pocket analysis, and protein preparation.\n"
            )
        elif average_plddt >= 70:
            handle.write(
                "The AlphaFold model is generally confident, but "
                "lower-confidence regions should be inspected before docking.\n"
            )
        else:
            handle.write(
                "The model contains substantial uncertainty and should not "
                "be used for docking without additional validation or "
                "structural refinement.\n"
            )

    print("Structure validation completed.")
    print(f"Models: {len(models)}")
    print(f"Chains: {', '.join(chain_ids)}")
    print(f"Residues: {len(residues)}")
    print(f"Atoms: {len(atoms)}")
    print(f"Average pLDDT: {average_plddt:.2f}")
    print(f"Minimum pLDDT: {minimum_plddt:.2f}")
    print(f"Maximum pLDDT: {maximum_plddt:.2f}")
    print(f"Saved: {tsv_file}")
    print(f"Saved: {markdown_file}")


if __name__ == "__main__":
    main()

