#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]

fasta = ROOT / "02_target_preparation/reference/MmSCD_P56221.fasta"

records = list(SeqIO.parse(fasta, "fasta"))

if len(records) != 2:
    raise RuntimeError("Expected exactly two protein sequences")

mm = records[0]
ref = records[1]

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

alignment = aligner.align(mm.seq, ref.seq)[0]

print("===== ALIGNMENT SCORE =====")
print(alignment.score)

print()
print("===== ALIGNMENT =====")
print(alignment)

# residue sets used in the project
catalytic_reference = [30, 31, 50, 85, 110]

binding_reference = [
    26,30,50,53,54,69,70,75,76,
    85,106,108,110,127,129,131,
    147,149,151,153,158,162,165,169
]

def build_mapping(aln):
    """
    Map residue numbers from reference sequence (target 1)
    to MmSCD sequence (target 0), based on aligned blocks.
    """
    mapping = {}

    mm_blocks, ref_blocks = aln.aligned

    for (mm_start, mm_end), (ref_start, ref_end) in zip(mm_blocks, ref_blocks):
        block_len = min(mm_end-mm_start, ref_end-ref_start)

        for i in range(block_len):
            # convert Python 0-based coordinates to residue numbering
            mm_pos = mm_start + i + 1
            ref_pos = ref_start + i + 1
            mapping[ref_pos] = mm_pos

    return mapping

mapping = build_mapping(alignment)

print()
print("===== CATALYTIC RESIDUE MAPPING =====")

for r in catalytic_reference:
    print(f"{r}\t{mapping.get(r, 'GAP')}")

print()
print("===== 24-RESIDUE POCKET MAPPING =====")

for r in binding_reference:
    print(f"{r}\t{mapping.get(r, 'GAP')}")


