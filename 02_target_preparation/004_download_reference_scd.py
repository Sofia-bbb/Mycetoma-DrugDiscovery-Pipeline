#!/usr/bin/env python3

from pathlib import Path
import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]

homology_dir = PROJECT_ROOT / "02_target_preparation" / "reference"
homology_dir.mkdir(parents=True, exist_ok=True)

downloads = {
    "P56221_reference_SCD.fasta":
        "https://rest.uniprot.org/uniprotkb/P56221.fasta",
    "2STD_carpropamid_complex.pdb":
        "https://files.rcsb.org/download/2STD.pdb",
}

for filename, url in downloads.items():
    print(f"Downloading {filename}...")
    response = requests.get(url, timeout=120)
    response.raise_for_status()

    output = homology_dir / filename
    output.write_bytes(response.content)

    print(f"Saved: {output}")

print("Reference retrieval completed.")

