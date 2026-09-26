from pathlib import Path
import xml.etree.ElementTree as ET
import csv
from collections import Counter

ROOT = Path("04_interaction_profiling")

PLIP = ROOT / "PLIP_Results"
OUT  = ROOT / "Analysis"

OUT.mkdir(parents=True, exist_ok=True)

# Key MmSCD residues used throughout the corrected workflow
KEY_RESIDUES = {
    28:  "TYR28",
    48:  "TYR48",
    83:  "HIS83",
    109: "HIS109",
    128: "SER128",
    130: "ASN130",
}

# PLIP XML containers and individual interaction tags
INTERACTION_TYPES = {
    "hydrophobic": (
        "hydrophobic_interactions",
        "hydrophobic_interaction"
    ),
    "hbond": (
        "hydrogen_bonds",
        "hydrogen_bond"
    ),
    "water_bridge": (
        "water_bridges",
        "water_bridge"
    ),
    "salt_bridge": (
        "salt_bridges",
        "salt_bridge"
    ),
    "pi_stack": (
        "pi_stacks",
        "pi_stack"
    ),
    "pi_cation": (
        "pi_cation_interactions",
        "pi_cation_interaction"
    ),
    "halogen_bond": (
        "halogen_bonds",
        "halogen_bond"
    ),
    "metal_complex": (
        "metal_complexes",
        "metal_complex"
    ),
}


def text_of(parent, tag):
    x = parent.find(tag)
    if x is None or x.text is None:
        return ""
    return x.text.strip()


def parse_residue(interaction):
    """
    Read the protein residue from a PLIP interaction record.
    Standard PLIP XML uses restype/resnr/reschain.
    """

    restype = text_of(interaction, "restype")
    resnr   = text_of(interaction, "resnr")
    chain   = text_of(interaction, "reschain")

    try:
        resnr_int = int(resnr)
    except:
        resnr_int = None

    if restype and resnr:
        residue = f"{restype}{resnr}"
    elif resnr:
        residue = resnr
    else:
        residue = ""

    return residue, resnr_int, chain


summary_rows = []
detail_rows = []

xml_files = sorted(PLIP.glob("*/*report.xml"))

print("XML reports found:", len(xml_files))

for xml_file in xml_files:

    ligand_id = xml_file.parent.name

    tree = ET.parse(xml_file)
    root = tree.getroot()

    counts = Counter()

    all_residues = set()
    key_residues = set()

    binding_sites = root.findall(".//bindingsite")

    has_interactions = False

    for bs in binding_sites:

        if bs.attrib.get("has_interactions", "").lower() == "true":
            has_interactions = True

        inter_root = bs.find("interactions")

        if inter_root is None:
            continue

        for interaction_type, (container_tag, item_tag) in INTERACTION_TYPES.items():

            container = inter_root.find(container_tag)

            if container is None:
                continue

            entries = container.findall(item_tag)

            counts[interaction_type] += len(entries)

            for entry in entries:

                residue, resnr, chain = parse_residue(entry)

                if residue:
                    all_residues.add(residue)

                is_key = (
                    resnr in KEY_RESIDUES
                    if resnr is not None
                    else False
                )

                if is_key:
                    key_residues.add(KEY_RESIDUES[resnr])

                # Retain useful PLIP geometry values when present
                distance = ""

                for possible in [
                    "dist",
                    "dist_h-a",
                    "dist_d-a",
                    "centdist",
                    "distance",
                ]:
                    value = text_of(entry, possible)
                    if value:
                        distance = value
                        break

                detail_rows.append({
                    "Ligand_ID": ligand_id,
                    "Interaction_Type": interaction_type,
                    "Residue": residue,
                    "Residue_Number": resnr if resnr is not None else "",
                    "Chain": chain,
                    "Key_MmSCD_Residue": "YES" if is_key else "NO",
                    "Distance_A": distance,
                })

    total_interactions = sum(counts.values())

    summary_rows.append({
        "Ligand_ID": ligand_id,

        "Hydrophobic": counts["hydrophobic"],
        "Hydrogen_Bonds": counts["hbond"],
        "Water_Bridges": counts["water_bridge"],
        "Salt_Bridges": counts["salt_bridge"],
        "Pi_Stacks": counts["pi_stack"],
        "Pi_Cation": counts["pi_cation"],
        "Halogen_Bonds": counts["halogen_bond"],
        "Metal_Complexes": counts["metal_complex"],

        "Total_PLIP_Interactions": total_interactions,

        "All_Interacting_Residues":
            ";".join(sorted(all_residues)),

        "Key_MmSCD_Residues":
            ";".join(
                sorted(
                    key_residues,
                    key=lambda x: int(
                        "".join(c for c in x if c.isdigit())
                    )
                )
            ),

        "Key_Residue_Count": len(key_residues),

        "Has_PLIP_Interactions":
            "YES" if has_interactions else "NO",

        "Reference_Status":
            (
                "REFERENCE"
                if ligand_id in {"Carpropamid", "Fenoxanil"}
                else "CANDIDATE"
            ),
    })


# ============================================================
# Write Supplementary Table S9
# ============================================================

summary_out = OUT / "Supplementary_Table_S9_PLIP_Master_Interactions.tsv"

summary_fields = [
    "Ligand_ID",
    "Reference_Status",
    "Hydrophobic",
    "Hydrogen_Bonds",
    "Water_Bridges",
    "Salt_Bridges",
    "Pi_Stacks",
    "Pi_Cation",
    "Halogen_Bonds",
    "Metal_Complexes",
    "Total_PLIP_Interactions",
    "All_Interacting_Residues",
    "Key_MmSCD_Residues",
    "Key_Residue_Count",
    "Has_PLIP_Interactions",
]

with open(summary_out, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=summary_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary_rows)


# ============================================================
# Detailed interaction-level table
# ============================================================

detail_out = OUT / "PLIP_Interaction_Detail.tsv"

detail_fields = [
    "Ligand_ID",
    "Interaction_Type",
    "Residue",
    "Residue_Number",
    "Chain",
    "Key_MmSCD_Residue",
    "Distance_A",
]

with open(detail_out, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=detail_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(detail_rows)


# ============================================================
# Ranked summary
# ============================================================

ranked = sorted(
    summary_rows,
    key=lambda r: (
        -r["Key_Residue_Count"],
        -r["Hydrogen_Bonds"],
        -r["Total_PLIP_Interactions"],
    )
)

rank_out = OUT / "PLIP_Ranked_By_Key_Interactions.tsv"

with open(rank_out, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=summary_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(ranked)


# ============================================================
# Console summary
# ============================================================

print()
print("===== PLIP MASTER ANALYSIS =====")
print("XML reports analyzed:", len(xml_files))
print("Summary rows:", len(summary_rows))
print("Detailed interaction rows:", len(detail_rows))

print()
print("REFERENCE COMPOUNDS")

for r in summary_rows:
    if r["Reference_Status"] == "REFERENCE":
        print(
            f'{r["Ligand_ID"]:12s} '
            f'Hbond={r["Hydrogen_Bonds"]:2d} '
            f'Hydrophobic={r["Hydrophobic"]:2d} '
            f'Total={r["Total_PLIP_Interactions"]:2d} '
            f'Key={r["Key_Residue_Count"]}/6 '
            f'{r["Key_MmSCD_Residues"]}'
        )

print()
print("TOP 20 BY KEY-RESIDUE / INTERACTION PROFILE")

for r in ranked[:20]:
    print(
        f'{r["Ligand_ID"]:14s} '
        f'Key={r["Key_Residue_Count"]}/6 '
        f'HB={r["Hydrogen_Bonds"]:2d} '
        f'Hyd={r["Hydrophobic"]:2d} '
        f'Total={r["Total_PLIP_Interactions"]:2d} '
        f'{r["Key_MmSCD_Residues"]}'
    )

print()
print("Saved:")
print(summary_out)
print(detail_out)
print(rank_out)


