#!/usr/bin/env python3

from Bio import AlignIO
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

alignment_file = (
    PROJECT_ROOT /
    "02_target_preparation" /
    "reference" /
    "MmSCD_P56221_aligned.fasta"
)

alignment = AlignIO.read(alignment_file, "fasta")

mm = alignment[0]
ref = alignment[1]

reference_residues = {
    30: "Y",
    31: "D",
    50: "Y",
    85: "H",
    110: "H"
}

ref_counter = 0
mm_counter = 0

print("\nCatalytic residue mapping\n")
print("Reference\tMmSCD\tReferenceAA\tMmSCDAA\tConserved")

for i in range(len(ref.seq)):

    if ref.seq[i] != "-":
        ref_counter += 1

    if mm.seq[i] != "-":
        mm_counter += 1

    if ref_counter in reference_residues:

        ref_aa = ref.seq[i]
        mm_aa = mm.seq[i]

        conserved = "YES" if ref_aa == mm_aa else "NO"

        print(
            f"{ref_counter}\t"
            f"{mm_counter}\t"
            f"{ref_aa}\t"
            f"{mm_aa}\t"
            f"{conserved}"
        )

