#!/usr/bin/env python3
from pathlib import Path

from Bio.PDB import PDBParser, NeighborSearch

pdb_file = Path(__file__).resolve().parents[1] / "02_target_preparation" / "reference" / "2STD_carpropamid_complex.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("2STD", pdb_file)

atoms = list(structure.get_atoms())
ns = NeighborSearch(atoms)

ligand_atoms = []

for residue in structure.get_residues():
    if residue.get_resname() == "CRP":
        ligand_atoms.extend(list(residue.get_atoms()))

binding_residues = set()

for atom in ligand_atoms:
    neighbors = ns.search(atom.coord, 4.5)

    for n in neighbors:
        res = n.get_parent()

        if res.get_resname() != "CRP":
            if res.id[0] == " ":
                binding_residues.add(
                    (
                        res.get_parent().id,
                        res.get_resname(),
                        res.id[1],
                    )
                )

print("\nResidues within 4.5 Å of carpropamid\n")

for chain, resname, number in sorted(binding_residues,
                                     key=lambda x: x[2]):
    print(f"{chain} {resname} {number}")

