#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


ROOT = Path("05_prioritization")

INPUT = ROOT / "Results" / "PLIP121_Brenk.csv"

OUTDIR = ROOT / "Results"

LOGDIR = ROOT / "Logs"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "PLIP121_Quinone_Chemistry.csv"
SUMMARY = OUTDIR / "PLIP121_Quinone_Chemistry_Summary.txt"
FAILURES = LOGDIR / "PLIP121_quinone_chemistry_failures.tsv"


# ------------------------------------------------------------
# Same structural patterns as original Phase 4.4
# ------------------------------------------------------------

SMARTS = {
    "Benzoquinone_1_4": "O=C1C=CC(=O)C=C1",
    "Naphthoquinone_1_4": "O=C1C=CC2=CC=CC=C2C1=O",
    "Anthraquinone": "O=C1c2ccccc2C(=O)c2ccccc12",
    "Hydroquinone": "Oc1ccc(O)cc1",
    "Michael_Acceptor": "[C,c]=[C,c]-[CX3](=O)[#6,#7,#8]",
    "Catechol": "Oc1ccccc1O",
}

PATTERNS = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in SMARTS.items()
}


def count_matches(molecule, pattern):
    if pattern is None:
        return 0

    return len(
        molecule.GetSubstructMatches(
            pattern,
            uniquify=True
        )
    )


def classify_quinone(feature_counts):

    if feature_counts["Anthraquinone"] > 0:
        return "Anthraquinone"

    if feature_counts["Naphthoquinone_1_4"] > 0:
        return "Naphthoquinone"

    if feature_counts["Benzoquinone_1_4"] > 0:
        return "Benzoquinone"

    if feature_counts["Hydroquinone"] > 0:
        return "Hydroquinone-like"

    return "Other_Quinone_or_Unclassified"


with INPUT.open(
    newline="",
    encoding="utf-8-sig"
) as f:

    rows = list(csv.DictReader(f))


print("Input compounds:", len(rows))


results = []
failures = []

candidate_classes = Counter()
reference_classes = Counter()

candidate_michael = 0
candidate_hydroquinone = 0
candidate_catechol = 0

reference_michael = 0
reference_hydroquinone = 0
reference_catechol = 0


for source_row in rows:

    ligand = source_row["Ligand"].strip()
    smiles = source_row["Canonical_SMILES"].strip()
    ctype = source_row["Compound_Type"].strip()

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        failures.append(
            (ligand, "Invalid SMILES", smiles)
        )
        continue

    feature_counts = {
        name: count_matches(
            molecule,
            pattern
        )
        for name, pattern in PATTERNS.items()
    }

    quinone_class = classify_quinone(
        feature_counts
    )

    michael_count = feature_counts[
        "Michael_Acceptor"
    ]

    hydroquinone_count = feature_counts[
        "Hydroquinone"
    ]

    catechol_count = feature_counts[
        "Catechol"
    ]

    if ctype == "Candidate":

        candidate_classes[quinone_class] += 1

        if michael_count > 0:
            candidate_michael += 1

        if hydroquinone_count > 0:
            candidate_hydroquinone += 1

        if catechol_count > 0:
            candidate_catechol += 1

    else:

        reference_classes[quinone_class] += 1

        if michael_count > 0:
            reference_michael += 1

        if hydroquinone_count > 0:
            reference_hydroquinone += 1

        if catechol_count > 0:
            reference_catechol += 1


    aromatic_rings = (
        rdMolDescriptors.CalcNumAromaticRings(
            molecule
        )
    )

    total_rings = rdMolDescriptors.CalcNumRings(
        molecule
    )


    result = dict(source_row)

    result.update(
        {
            "Quinone_Class": quinone_class,

            "Benzoquinone_Count":
                feature_counts["Benzoquinone_1_4"],

            "Naphthoquinone_Count":
                feature_counts["Naphthoquinone_1_4"],

            "Anthraquinone_Count":
                feature_counts["Anthraquinone"],

            "Hydroquinone_Count":
                hydroquinone_count,

            "Catechol_Count":
                catechol_count,

            "Michael_Acceptor_Count":
                michael_count,

            "Michael_Acceptor_Flag":
                int(michael_count > 0),

            "Hydroquinone_Flag":
                int(hydroquinone_count > 0),

            "Catechol_Flag":
                int(catechol_count > 0),

            "Aromatic_Ring_Count_404":
                aromatic_rings,

            "Total_Ring_Count_404":
                total_rings,

            # Candidate quinone membership is expected.
            # References are retained regardless of class.
            "Quinone_Scaffold_Expected":
                int(ctype == "Candidate"),

            "Reactivity_Review_Flag":
                int(
                    michael_count > 0
                    or catechol_count > 0
                ),

            # Do not eliminate automatically here.
            "Retain_For_Further_Analysis": 1,
        }
    )

    results.append(result)


if not results:
    raise RuntimeError(
        "No compounds were successfully processed."
    )


fields = list(results[0].keys())

with OUTPUT.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(results)


with FAILURES.open(
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "Ligand\tReason\tStructure\n"
    )

    for failure in failures:
        f.write(
            "\t".join(failure) + "\n"
        )


references = [
    r for r in results
    if r["Compound_Type"] == "Reference"
]


summary = [
    "===== PLIP121 QUINONE CHEMISTRY =====",
    f"Total input: {len(rows)}",
    f"Successfully processed: {len(results)}",
    f"Failures: {len(failures)}",
    "",
    "CANDIDATES",
]

for qclass, count in candidate_classes.most_common():
    summary.append(
        f"{qclass}: {count}"
    )

summary.extend(
    [
        "",
        f"Michael-acceptor positive: {candidate_michael}",
        f"Hydroquinone positive: {candidate_hydroquinone}",
        f"Catechol positive: {candidate_catechol}",
        "",
        "REFERENCES",
    ]
)

for r in references:

    summary.append(
        f'{r["Ligand"]}: '
        f'Class={r["Quinone_Class"]}; '
        f'Michael={r["Michael_Acceptor_Flag"]}; '
        f'Hydroquinone={r["Hydroquinone_Flag"]}; '
        f'Catechol={r["Catechol_Flag"]}; '
        f'ReactivityReview={r["Reactivity_Review_Flag"]}'
    )


summary.extend(
    [
        "",
        "Interpretation:",
        (
            "Quinone scaffold membership is descriptive and "
            "is not an automatic exclusion criterion."
        ),
        (
            "Michael-acceptor and catechol motifs are retained "
            "as reactivity-review flags for later integration."
        ),
    ]
)


SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)


print("\n".join(summary))
print()
print("Saved:", OUTPUT)


