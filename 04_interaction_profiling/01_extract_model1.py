from pathlib import Path

ROOT = Path("04_interaction_profiling")

IDS = Path("03_docking/Advanced_Candidates.tsv")
OUT = ROOT / "Ligands"

POSES = Path(
    "03_docking/poses"
)

BENCH = Path(
    "02_target_preparation/benchmark"
)

OUT.mkdir(parents=True, exist_ok=True)


def extract_model1(src, dst):
    """Extract MODEL 1 only from a Vina multi-model PDBQT."""

    if not src.exists():
        return False, "SOURCE_MISSING"

    lines = []
    inside = False
    found_model = False

    with open(src) as f:
        for line in f:

            if line.startswith("MODEL"):
                parts = line.split()

                if len(parts) >= 2 and parts[1] == "1":
                    inside = True
                    found_model = True
                    continue

                elif inside:
                    break

            if inside and line.startswith("ENDMDL"):
                break

            if inside:
                lines.append(line)

    if not found_model:
        return False, "MODEL1_NOT_FOUND"

    atoms = sum(
        1 for line in lines
        if line.startswith(("ATOM", "HETATM"))
    )

    if atoms == 0:
        return False, "NO_ATOMS"

    with open(dst, "w") as f:
        f.writelines(lines)

    return True, atoms


# -------------------------
# 1188 advanced candidates
# -------------------------

ids = [
    line.split("\t")[1]
    for line in IDS.read_text().splitlines()[1:]
    if line.strip()
]

success = 0
failed = []

for name in ids:

    src = POSES / f"{name}_docked.pdbqt"
    dst = OUT / f"{name}_model1.pdbqt"

    ok, info = extract_model1(src, dst)

    if ok:
        success += 1
    else:
        failed.append((name, info))


# -------------------------
# Corrected benchmarks
# -------------------------

benchmarks = ["Carpropamid", "Fenoxanil"]

benchmark_success = 0
benchmark_failed = []

for name in benchmarks:

    src = BENCH / f"{name}_docked.pdbqt"
    dst = OUT / f"{name}_model1.pdbqt"

    ok, info = extract_model1(src, dst)

    if ok:
        benchmark_success += 1
    else:
        benchmark_failed.append((name, info))


print("===== MODEL 1 EXTRACTION =====")
print("Candidates expected:", len(ids))
print("Candidates extracted:", success)
print("Candidate failures:", len(failed))

print()
print("Benchmarks expected:", len(benchmarks))
print("Benchmarks extracted:", benchmark_success)
print("Benchmark failures:", len(benchmark_failed))

if failed:
    print("\nFAILED CANDIDATES")
    for x in failed:
        print(*x)

if benchmark_failed:
    print("\nFAILED BENCHMARKS")
    for x in benchmark_failed:
        print(*x)

print("\nTotal MODEL 1 files:",
      len(list(OUT.glob("*_model1.pdbqt"))))


