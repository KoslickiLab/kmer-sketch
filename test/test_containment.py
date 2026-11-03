"""
test_containment.py
===================
Verify that `filter --metric containment` returns the correct directional
containment estimate for FracMinHash and BottomK sketches.

Before this fix, both sketch types silently fell back to Jaccard similarity
whenever `--metric containment` was requested.

Test strategy
-------------
We construct two synthetic genomes A and B with known k-mer overlap:

  Symmetric case  (k=31):
    |A| = 10 000 unique k-mers
    |B| = 10 000 unique k-mers, 5 000 shared with A
    Expected: C(A⊆B) = C(B⊆A) = 0.5

  Asymmetric case:
    |A| =  5 000 unique k-mers (all are in B)
    |B| = 10 000 unique k-mers (5 000 shared with A)
    Expected: C(A⊆B) ≈ 1.0,  C(B⊆A) ≈ 0.5

Each k-mer is written as its own 31-bp FASTA record (one sequence = one k-mer),
so the genome has exactly the k-mers we specify and the ground truth is exact.

Run from the repo root (scripts/kmer-sketch/):
  python3 test/test_containment.py
"""

import os
import random
import subprocess
import tempfile

# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
KMER_SIZE   = 31
BIN_SKETCH  = "bin/sketch"
BIN_FILTER  = "bin/filter"

FMH_SCALE   = 0.1   # ~10% sampling → ~500-1000 hashes; enough for stable estimates
BK_K        = 1000  # keep 1000 bottom hashes; plenty for 5 000-10 000 k-mer genomes

TOL_SYM     = 0.05  # ±5 percentage points for symmetric test
TOL_ASYM_HI = 0.05  # ±5 pp for the near-1.0 direction
TOL_ASYM_LO = 0.05  # ±5 pp for the near-0.5 direction

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _rand_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _unique_kmers(n: int, k: int, seed: int) -> list[str]:
    """Return n unique, random, canonical-free DNA k-mers of length k."""
    rng = random.Random(seed)
    seen: set[str] = set()
    result: list[str] = []
    while len(result) < n:
        kmer = _rand_dna(k, rng)
        if kmer not in seen:
            seen.add(kmer)
            result.append(kmer)
    return result


def _write_fasta(path: str, kmers: list[str]) -> None:
    """Write one FASTA record per k-mer (header = >kmer_<i>)."""
    with open(path, "w") as fh:
        for i, kmer in enumerate(kmers):
            fh.write(f">kmer_{i}\n{kmer}\n")


def _sketch(fasta: str, sketch: str, algo: str, **params) -> None:
    """Run bin/sketch to produce a sketch file."""
    cmd = [
        BIN_SKETCH,
        "--input",     fasta,
        "--kmer",      str(KMER_SIZE),
        "--algo",      algo,
        "--canonical",
        "--output",    sketch,
    ]
    for key, val in params.items():
        cmd += [f"--{key}", str(val)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"sketch failed:\n{r.stderr}")


def _filter_containment(query_sketch: str, ref_sketch: str) -> float:
    """
    Run `filter --metric containment` and return the containment score.
    The binary computes C(reference ⊆ query), i.e. the fraction of the
    *reference* sketch found in the *query* sketch.
    """
    r = subprocess.run(
        [
            BIN_FILTER,
            "--query",     query_sketch,
            "--refs",      ref_sketch,
            "--metric",    "containment",
            "--threshold", "0.0",
        ],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"filter failed:\n{r.stderr}")
    # Output: header line + one data line
    lines = [l for l in r.stdout.splitlines() if l and not l.startswith("#")]
    # Skip the header row
    data_lines = lines[1:] if lines and lines[0].startswith("reference") else lines
    if not data_lines:
        raise RuntimeError(f"No output from filter:\nstdout={r.stdout!r}\nstderr={r.stderr!r}")
    parts = data_lines[0].split("\t")
    return float(parts[1])


def _filter_jaccard(query_sketch: str, ref_sketch: str) -> float:
    """Run `filter --metric jaccard` and return the score."""
    r = subprocess.run(
        [
            BIN_FILTER,
            "--query",     query_sketch,
            "--refs",      ref_sketch,
            "--metric",    "jaccard",
            "--threshold", "0.0",
        ],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"filter failed:\n{r.stderr}")
    lines = [l for l in r.stdout.splitlines() if l and not l.startswith("#")]
    data_lines = lines[1:] if lines and lines[0].startswith("reference") else lines
    if not data_lines:
        raise RuntimeError(f"No output from filter:\nstdout={r.stdout!r}\nstderr={r.stderr!r}")
    parts = data_lines[0].split("\t")
    return float(parts[1])


# --------------------------------------------------------------------------- #
# Test cases
# --------------------------------------------------------------------------- #

def run_tests(algo: str, sketch_fn, tmpdir: str) -> None:
    """
    Run symmetric and asymmetric containment tests for a given algorithm.

    algo      – display name (e.g. "FracMinHash", "BottomK")
    sketch_fn – callable(fasta_path, sketch_path) that writes the sketch
    tmpdir    – temporary directory for intermediate files
    """
    print(f"\n{'='*60}")
    print(f"  Algorithm: {algo}")
    print(f"{'='*60}")

    # ---- Build k-mer pools ------------------------------------------------- #
    # 15 000 unique k-mers; first 5 000 = shared, next 5 000 = A-only,
    # next 5 000 = B-only.
    all_kmers  = _unique_kmers(15_000, KMER_SIZE, seed=2024)
    shared     = all_kmers[:5_000]
    a_only     = all_kmers[5_000:10_000]
    b_only     = all_kmers[10_000:15_000]

    # Symmetric case: A = shared + a_only (10 000), B = shared + b_only (10 000)
    kmers_A_sym = shared + a_only   # 10 000 k-mers
    kmers_B_sym = shared + b_only   # 10 000 k-mers

    # Asymmetric case: A_asym = shared (5 000), B_asym = shared + b_only (10 000)
    kmers_A_asym = shared            # 5 000 k-mers (all in B_asym)
    kmers_B_asym = shared + b_only   # 10 000 k-mers

    fasta_A_sym  = os.path.join(tmpdir, f"{algo}_A_sym.fasta")
    fasta_B_sym  = os.path.join(tmpdir, f"{algo}_B_sym.fasta")
    fasta_A_asym = os.path.join(tmpdir, f"{algo}_A_asym.fasta")
    fasta_B_asym = os.path.join(tmpdir, f"{algo}_B_asym.fasta")

    _write_fasta(fasta_A_sym,  kmers_A_sym)
    _write_fasta(fasta_B_sym,  kmers_B_sym)
    _write_fasta(fasta_A_asym, kmers_A_asym)
    _write_fasta(fasta_B_asym, kmers_B_asym)

    sk_A_sym  = os.path.join(tmpdir, f"{algo}_A_sym.sketch")
    sk_B_sym  = os.path.join(tmpdir, f"{algo}_B_sym.sketch")
    sk_A_asym = os.path.join(tmpdir, f"{algo}_A_asym.sketch")
    sk_B_asym = os.path.join(tmpdir, f"{algo}_B_asym.sketch")

    sketch_fn(fasta_A_sym,  sk_A_sym)
    sketch_fn(fasta_B_sym,  sk_B_sym)
    sketch_fn(fasta_A_asym, sk_A_asym)
    sketch_fn(fasta_B_asym, sk_B_asym)

    # ---- Symmetric test ----------------------------------------------------- #
    # filter --query A --refs B computes C(B ⊆ A) = |B ∩ A| / |B|
    # filter --query B --refs A computes C(A ⊆ B) = |A ∩ B| / |A|
    # Both should be ≈ 0.5 in the symmetric case.
    c_b_in_a_sym = _filter_containment(sk_A_sym, sk_B_sym)
    c_a_in_b_sym = _filter_containment(sk_B_sym, sk_A_sym)
    j_sym        = _filter_jaccard(sk_A_sym, sk_B_sym)

    expected_cont_sym = 0.5
    expected_jacc_sym = 5_000 / (10_000 + 10_000 - 5_000)   # 5000/15000 ≈ 0.333

    print(f"\n  [Symmetric]  |A|=|B|=10000, |A∩B|=5000")
    print(f"    C(B⊆A) = {c_b_in_a_sym:.4f}  (expected ≈ {expected_cont_sym:.4f})")
    print(f"    C(A⊆B) = {c_a_in_b_sym:.4f}  (expected ≈ {expected_cont_sym:.4f})")
    print(f"    Jaccard = {j_sym:.4f}  (expected ≈ {expected_jacc_sym:.4f})")

    assert abs(c_b_in_a_sym - expected_cont_sym) < TOL_SYM, (
        f"{algo} symmetric: C(B⊆A) = {c_b_in_a_sym:.4f}, expected {expected_cont_sym:.4f} ± {TOL_SYM}"
    )
    assert abs(c_a_in_b_sym - expected_cont_sym) < TOL_SYM, (
        f"{algo} symmetric: C(A⊆B) = {c_a_in_b_sym:.4f}, expected {expected_cont_sym:.4f} ± {TOL_SYM}"
    )
    # Sanity: containment must differ from Jaccard in the symmetric case
    assert abs(c_b_in_a_sym - j_sym) > 0.05, (
        f"{algo} symmetric: containment ({c_b_in_a_sym:.4f}) ≈ jaccard ({j_sym:.4f}); "
        "looks like the fix did not take effect"
    )
    print("    ✓ symmetric containment correct")

    # ---- Asymmetric test ---------------------------------------------------- #
    # A_asym ⊆ B_asym entirely → C(A_asym ⊆ B_asym) ≈ 1.0
    # filter --query B_asym --refs A_asym  →  C(A_asym ⊆ B_asym)
    c_a_in_b_asym = _filter_containment(sk_B_asym, sk_A_asym)
    # C(B_asym ⊆ A_asym) = |shared| / |B_asym| = 5000/10000 = 0.5
    # filter --query A_asym --refs B_asym  →  C(B_asym ⊆ A_asym)
    c_b_in_a_asym = _filter_containment(sk_A_asym, sk_B_asym)

    expected_a_in_b = 1.0
    expected_b_in_a = 0.5

    print(f"\n  [Asymmetric] |A|=5000 (all in B), |B|=10000")
    print(f"    C(A⊆B) = {c_a_in_b_asym:.4f}  (expected ≈ {expected_a_in_b:.4f})")
    print(f"    C(B⊆A) = {c_b_in_a_asym:.4f}  (expected ≈ {expected_b_in_a:.4f})")

    assert abs(c_a_in_b_asym - expected_a_in_b) < TOL_ASYM_HI, (
        f"{algo} asymmetric: C(A⊆B) = {c_a_in_b_asym:.4f}, expected {expected_a_in_b:.4f} ± {TOL_ASYM_HI}"
    )
    assert abs(c_b_in_a_asym - expected_b_in_a) < TOL_ASYM_LO, (
        f"{algo} asymmetric: C(B⊆A) = {c_b_in_a_asym:.4f}, expected {expected_b_in_a:.4f} ± {TOL_ASYM_LO}"
    )
    # Sanity: the two containment directions must differ significantly
    assert abs(c_a_in_b_asym - c_b_in_a_asym) > 0.3, (
        f"{algo} asymmetric: containment values are symmetric "
        f"({c_a_in_b_asym:.4f} vs {c_b_in_a_asym:.4f}); expected asymmetry > 0.3"
    )
    print("    ✓ asymmetric containment correct")

    print(f"\n  {algo}: ALL TESTS PASSED")


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    os.makedirs("outputs", exist_ok=True)

    with tempfile.TemporaryDirectory(dir="outputs", prefix="test_containment_") as tmpdir:

        # ---- FracMinHash ---------------------------------------------------- #
        def sketch_fmh(fasta: str, out: str) -> None:
            _sketch(fasta, out, "fracminhash", scale=FMH_SCALE, seed=42)

        run_tests("FracMinHash", sketch_fmh, tmpdir)

        # ---- BottomK -------------------------------------------------------- #
        def sketch_bk(fasta: str, out: str) -> None:
            _sketch(fasta, out, "bottomk", k=BK_K, seed=42)

        run_tests("BottomK", sketch_bk, tmpdir)

    print("\n" + "="*60)
    print("  ALL CONTAINMENT TESTS PASSED")
    print("="*60)
