#!/usr/bin/env python3

"""
001_download_uniprot.py

Mycetoma Drug Discovery Platform

Purpose
-------
Download the UniProt sequence and annotation for the target protein.

Input
-----
config/project.yaml

Outputs
-------
02_target_preparation/reference/MmSCD.fasta
02_target_preparation/reference/MmSCD.json
"""

from pathlib import Path
import requests
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "project.yaml"

with CONFIG_FILE.open("r", encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

accession = config["target"]["uniprot_accession"]

sequence_dir = PROJECT_ROOT / "02_target_preparation" / "reference"
sequence_dir.mkdir(parents=True, exist_ok=True)

print(f"Target accession: {accession}")

downloads = {
    "MmSCD.fasta": f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
    "MmSCD.json": f"https://rest.uniprot.org/uniprotkb/{accession}.json",
}

downloaded = 0
skipped = 0

for filename, url in downloads.items():

    output_file = sequence_dir / filename

    if output_file.exists():
        print(f"Skipping {filename} (already exists)")
        skipped += 1
        continue

    print(f"Downloading {filename}...")

    response = requests.get(url, timeout=60)
    response.raise_for_status()

    output_file.write_text(response.text, encoding="utf-8")

    print(
        f"Saved: {output_file.name} "
        f"({len(response.text):,} characters)"
    )

    downloaded += 1

print("=" * 60)
print("UniProt retrieval completed successfully.")
print(f"Target accession : {accession}")
print(f"Downloaded       : {downloaded}")
print(f"Skipped          : {skipped}")
print("=" * 60)


