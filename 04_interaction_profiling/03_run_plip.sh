#!/bin/bash
set -e

ROOT="04_interaction_profiling"
COMPLEXDIR="$ROOT/Complexes"
OUTROOT="$ROOT/PLIP_Results"

mkdir -p "$OUTROOT"

for complex in "$COMPLEXDIR"/*_complex.pdb
do
    name=$(basename "$complex" _complex.pdb)
    out="$OUTROOT/$name"

    mkdir -p "$out"

    echo "Running PLIP: $name"

    plip \
      -f "$complex" \
      -o "$out" \
      -x \
      -t \
      --silent
done

echo "PLIP complete."

