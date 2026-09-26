#!/usr/bin/env python3

from pathlib import Path
import csv

ROOT = Path("05_prioritization")

PRED = ROOT / "ADMET" / "Results" / "PLIP121_ADMET_AI_Predictions.csv"
EVIDENCE = ROOT / "Results" / "PLIP121_Developability_Profile.csv"
OUTDIR = ROOT / "ADMET" / "Results"

OUTDIR.mkdir(parents=True, exist_ok=True)

MASTER = OUTDIR / "PLIP121_ADMET_Integrated_Master.csv"
INTERPRET = OUTDIR / "PLIP121_ADMET_Interpretation.csv"
SUMMARY = OUTDIR / "PLIP121_ADMET_Interpretation_Summary.txt"


# ------------------------------------------------------------
# Directional ADMET endpoints
# ------------------------------------------------------------

LOWER_BETTER = [
    "AMES",
    "Carcinogens_Lagunin",
    "ClinTox",
    "DILI",
    "hERG",
    "Skin_Reaction",
    "CYP1A2_Veith",
    "CYP2C19_Veith",
    "CYP2C9_Veith",
    "CYP2D6_Veith",
    "CYP3A4_Veith",
    "Pgp_Broccatelli",
]

HIGHER_BETTER = [
    "HIA_Hou",
    "Bioavailability_Ma",
    "PAMPA_NCATS",
]

CONTEXTUAL = [
    "BBB_Martins",
    "Caco2_Wang",
    "PPBR_AZ",
    "VDss_Lombardo",
    "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ",
    "Half_Life_Obach",
    "LD50_Zhu",
    "Solubility_AqSolDB",
    "Lipophilicity_AstraZeneca",
    "HydrationFreeEnergy_FreeSolv",
]


def num(x):
    try:
        return float(x)
    except:
        return None


def status(value, lower_better=True):

    if value is None:
        return "MISSING"

    if lower_better:

        if value < 0.30:
            return "FAVORABLE"

        if value <= 0.70:
            return "INTERMEDIATE"

        return "UNFAVORABLE"

    else:

        if value > 0.70:
            return "FAVORABLE"

        if value >= 0.30:
            return "INTERMEDIATE"

        return "UNFAVORABLE"


# ------------------------------------------------------------
# Load platform evidence
# ------------------------------------------------------------

with EVIDENCE.open(
    newline="",
    encoding="utf-8-sig"
) as f:

    evidence_rows = list(csv.DictReader(f))


evidence = {
    r["Ligand"].strip(): r
    for r in evidence_rows
}


# ------------------------------------------------------------
# Load ADMET predictions
# ------------------------------------------------------------

with PRED.open(
    newline="",
    encoding="utf-8-sig"
) as f:

    reader = csv.DictReader(f)
    predictions = list(reader)
    prediction_fields = reader.fieldnames or []


print("ADMET predictions:", len(predictions))
print("Evidence rows:", len(evidence))


# ------------------------------------------------------------
# Full integrated master
# ------------------------------------------------------------

master_rows = []

for pred in predictions:

    compound = pred["Compound"].strip()

    ev = evidence.get(compound)

    if ev is None:
        raise RuntimeError(
            f"Missing platform evidence for {compound}"
        )

    result = dict(ev)

    for key, value in pred.items():

        if key in {"Compound", "Type", "smiles"}:
            continue

        result[f"ADMET_{key}"] = value

    master_rows.append(result)


fields = list(master_rows[0].keys())

with MASTER.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(master_rows)


# ------------------------------------------------------------
# ADMET interpretation
# ------------------------------------------------------------

interpret_rows = []

for pred in predictions:

    compound = pred["Compound"].strip()
    ev = evidence[compound]

    favorable = 0
    intermediate = 0
    unfavorable = 0
    evaluated = 0

    strengths = []
    concerns = []

    result = {
        "Compound": compound,
        "Type": pred["Type"],

        "Docking_Score":
            ev.get("Docking_Score", ""),

        "PLIP_Key_Residue_Count":
            ev.get("PLIP_Key_Residue_Count", ""),

        "PLIP_Key_Residues":
            ev.get("PLIP_Key_Residues", ""),

        "Hydrogen_Bonds":
            ev.get("Hydrogen_Bonds", ""),

        "Total_PLIP_Interactions":
            ev.get("Total_PLIP_Interactions", ""),

        "QED":
            ev.get("QED", ""),

        "ESOL_LogS":
            ev.get("ESOL_LogS", ""),

        "PAINS_Count":
            ev.get("PAINS_Count", ""),

        "Brenk_Count":
            ev.get("Brenk_Count", ""),

        "Michael_Acceptor_Flag":
            ev.get("Michael_Acceptor_Flag", ""),

        "Catechol_Flag":
            ev.get("Catechol_Flag", ""),
    }


    for endpoint in LOWER_BETTER:

        if endpoint not in prediction_fields:
            continue

        value = num(pred.get(endpoint))

        s = status(
            value,
            lower_better=True
        )

        result[endpoint] = pred.get(
            endpoint,
            ""
        )

        result[f"{endpoint}_Status"] = s

        if s == "MISSING":
            continue

        evaluated += 1

        if s == "FAVORABLE":
            favorable += 1
            strengths.append(endpoint)

        elif s == "INTERMEDIATE":
            intermediate += 1

        else:
            unfavorable += 1
            concerns.append(endpoint)


    for endpoint in HIGHER_BETTER:

        if endpoint not in prediction_fields:
            continue

        value = num(pred.get(endpoint))

        s = status(
            value,
            lower_better=False
        )

        result[endpoint] = pred.get(
            endpoint,
            ""
        )

        result[f"{endpoint}_Status"] = s

        if s == "MISSING":
            continue

        evaluated += 1

        if s == "FAVORABLE":
            favorable += 1
            strengths.append(endpoint)

        elif s == "INTERMEDIATE":
            intermediate += 1

        else:
            unfavorable += 1
            concerns.append(endpoint)


    for endpoint in CONTEXTUAL:

        if endpoint in prediction_fields:
            result[endpoint] = pred.get(
                endpoint,
                ""
            )


    if evaluated:

        score = (
            favorable
            + 0.5 * intermediate
        ) / evaluated * 100

    else:

        score = None


    critical = sum(

        result.get(
            f"{endpoint}_Status"
        ) == "UNFAVORABLE"

        for endpoint in [
            "AMES",
            "ClinTox",
            "DILI",
            "hERG",
        ]

    )


    if critical == 0:

        decision = "ADVANCE"

    elif critical == 1:

        decision = "ADVANCE_WITH_CAUTION"

    else:

        decision = "DEPRIORITIZE"


    result.update({

        "Directional_Endpoints_Evaluated":
            evaluated,

        "Favorable_Count":
            favorable,

        "Intermediate_Count":
            intermediate,

        "Unfavorable_Count":
            unfavorable,

        "Directional_ADMET_Score":
            round(score, 3)
            if score is not None
            else "",

        "Critical_Safety_Concerns":
            critical,

        "ADMET_Decision":
            decision,

        "ADMET_Strengths":
            ";".join(strengths),

        "ADMET_Concerns":
            ";".join(concerns),
    })


    interpret_rows.append(result)


# ------------------------------------------------------------
# Rank candidates descriptively
# ADMET first, then target engagement
# ------------------------------------------------------------

candidates = [
    r for r in interpret_rows
    if r["Type"] == "Candidate"
]

references = [
    r for r in interpret_rows
    if r["Type"] != "Candidate"
]


def fnum(x, default=0.0):
    try:
        return float(x)
    except:
        return default


candidates.sort(
    key=lambda r: (
        r["Critical_Safety_Concerns"],
        -fnum(r["Directional_ADMET_Score"]),
        -fnum(r["PLIP_Key_Residue_Count"]),
        fnum(r["Docking_Score"], 999),
        -fnum(r["Hydrogen_Bonds"]),
    )
)


rank = 1

for r in candidates:
    r["Integrated_Priority_Rank"] = rank
    rank += 1

for r in references:
    r["Integrated_Priority_Rank"] = "REFERENCE"


final_rows = candidates + references

fields = list(final_rows[0].keys())

with INTERPRET.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(final_rows)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

advance = sum(
    r["ADMET_Decision"] == "ADVANCE"
    for r in candidates
)

caution = sum(
    r["ADMET_Decision"] == "ADVANCE_WITH_CAUTION"
    for r in candidates
)

deprioritize = sum(
    r["ADMET_Decision"] == "DEPRIORITIZE"
    for r in candidates
)


summary = [
    "===== PLIP121 INTEGRATED ADMET =====",
    f"Candidates: {len(candidates)}",
    f"References: {len(references)}",
    "",
    f"ADVANCE: {advance}",
    f"ADVANCE_WITH_CAUTION: {caution}",
    f"DEPRIORITIZE: {deprioritize}",
    "",
    "REFERENCES",
]


for r in references:

    summary.append(
        f'{r["Compound"]}: '
        f'ADMET={r["ADMET_Decision"]}; '
        f'Score={r["Directional_ADMET_Score"]}; '
        f'Critical={r["Critical_Safety_Concerns"]}; '
        f'Docking={r["Docking_Score"]}; '
        f'PLIP={r["PLIP_Key_Residue_Count"]}/6'
    )


summary.extend([
    "",
    "TOP 20 CANDIDATES",
])


for r in candidates[:20]:

    summary.append(
        f'{r["Integrated_Priority_Rank"]:2d}. '
        f'{r["Compound"]:14s} '
        f'ADMET={r["ADMET_Decision"]:20s} '
        f'Score={fnum(r["Directional_ADMET_Score"]):6.1f} '
        f'Critical={r["Critical_Safety_Concerns"]} '
        f'PLIP={r["PLIP_Key_Residue_Count"]}/6 '
        f'Dock={fnum(r["Docking_Score"]):7.3f} '
        f'HB={r["Hydrogen_Bonds"]}'
    )


SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)


print("\n".join(summary))
print()
print("Master:", MASTER)
print("Interpretation:", INTERPRET)


