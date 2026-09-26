#!/usr/bin/env python3

from pathlib import Path
import requests
import yaml

# Project root
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Read configuration
CONFIG = PROJECT_ROOT / "config" / "project.yaml"

with open(CONFIG, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

accession = cfg["target"]["uniprot_accession"]

print(f"Checking AlphaFold model for {accession}...")

api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"

response = requests.get(api_url, timeout=60)

if response.status_code == 404:
    print("No AlphaFold model found for this UniProt accession.")
    raise SystemExit()

response.raise_for_status()

data = response.json()

pdb_url = data[0]["pdbUrl"]

print("AlphaFold model found!")
print(f"Downloading structure from:\n{pdb_url}")

structure_dir = PROJECT_ROOT / "02_target_preparation" / "receptor"
structure_dir.mkdir(parents=True, exist_ok=True)

outfile = structure_dir / "MmSCD_AF.pdb"

if outfile.exists():
    print("AlphaFold structure already exists. Skipping download.")
    raise SystemExit(0)

pdb_response = requests.get(pdb_url, timeout=120)
pdb_response.raise_for_status()

outfile.write_text(
    pdb_response.text,
    encoding="utf-8",
)

print(f"\nStructure saved to:\n{outfile}")


