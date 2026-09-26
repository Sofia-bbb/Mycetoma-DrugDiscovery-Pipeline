#!/usr/bin/env python3

from pathlib import Path
from Bio import AlignIO

alignment = AlignIO.read(
    str(Path(__file__).resolve().parents[1] / "02_target_preparation" / "reference" / "MmSCD_P56221_aligned.fasta"),
    "fasta",
)

mm = alignment[0].seq
ref = alignment[1].seq

reference_binding = [
26,30,50,53,54,69,70,75,76,
85,106,108,110,127,129,
131,147,149,151,153,
158,162,165,169
]

mapping = {}

mm_pos = 0
ref_pos = 0

for i in range(len(mm)):

    if mm[i] != "-":
        mm_pos += 1

    if ref[i] != "-":
        ref_pos += 1

    if ref_pos in reference_binding and ref[i] != "-":

        if mm[i] != "-":

            mapping[ref_pos] = (
                mm_pos,
                ref[i],
                mm[i]
            )

print("\nReference\tMmSCD\tReferenceAA\tMmSCDAA")

for refpos in reference_binding:

    if refpos in mapping:

        mmpos, aa1, aa2 = mapping[refpos]

        print(
            f"{refpos}\t{mmpos}\t{aa1}\t{aa2}"
        )

