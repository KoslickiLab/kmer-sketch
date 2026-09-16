# MaxGeomHash & $\alpha$-MaxGeomHash  
A C++ implementation of MaxGeomHash and $\alpha$-MaxGeomHash for 
- k-mer sketching of FASTA/FASTQ files, and 
- pairwise similarity estimation from k-mer sketches.

## Table of Contents
1. [Overview](#overview)
2. [Build Instructions](#build-instructions)
   - [Building](#building)
   - [Verifying Installation](#verifying-installation)
3. [Quick start and demo](#quick-start-and-demo)
4. [Detailed usage](#detailed-usage)
   - [`sketch`](#sketch)
   - [`pwsimilarity`](#pwsimilarity)
   - [`expt_growth`](#expt_growth)
5. [Citing](#citing)

## Overview
This repository provides efficient C++ implementations of **MaxGeomHash** and **$\alpha$-MaxGeomHash**, two hashing–based sketching algorithms. The tools support several other sketching methods -- including FracMinHash, MinHash (bottom-k). From FASTA or FASTQ files, the tools can compute k-mer sketches using any of these sketching methods. The tools also allow for rapid pairwise similarity estimation using a number of sketches.

All dependencies are either header-only or included with the repository. Compilation requires only `make`.

---

## Build Instructions

### Building

Clone the repository and build. The executables are generated in the bin/ directory. Kindly make sure to add the bin/ directory in your PATH variable.

```bash
git clone https://github.com/mahmudhera/kmer-sketch
cd kmer-sketch
make
export PATH=$(pwd)/bin:$PATH
```

### Verifying Installation

Run `sketch -h` to verify installation.

## Quick start and demo
Following is a simple example showing how to compute MaxGeomhash sketches of three plant genomes, and how to compute their pairwise k-mer based Jaccard similarity scores.

```bash
# computing sketches
sketch --input data/ecoli-k12.fasta --kmer 31 --algo maxgeom --b 90 --w 64 --seed 42 --canonical --output data/ecoli-k12.maxgeom.sketch
sketch --input data/ecoli-nctc.fasta --kmer 31 --algo maxgeom --b 90 --w 64 --seed 42 --canonical --output data/ecoli-nctc.maxgeom.sketch
sketch --input data/ecoli-o157.fasta --kmer 31 --algo maxgeom --b 90 --w 64 --seed 42 --canonical --output data/ecoli-o157.maxgeom.sketch

# computing pairwise similarity scores
pwsimilarity --metric jaccard --output data/ecoli_pw_scores.tsv data/*.maxgeom.sketch

# view pairwise similarity scores
cat data/ecoli_pw_scores.tsv
```

## Detailed usage

Two main programs are included:
- `sketch`
- `pwsimilarity`

### `sketch`
This program creates the sketches. Following are the arguments.

| Argument            | Argument Type | Default Value | What It Means |
|---------------------|---------------|----------------|----------------|
| `--input FILE`      | string (path) | **required**   | Input sequence file in FASTA or FASTQ format. |
| `--kmer N`          | integer       | 31             | k-mer size to break sequences into. |
| `--algo ALGO`       | string        | **required**   | Sketching algorithm to use: `maxgeom`, `alphamaxgeom`, `fracminhash`, `minhash`, `bottomk`. |
| `--k K`             | integer       | 1000           | (bottom-k) Sketch size *K* for bottom-k hashing. |
| `--b B`             | integer       | 90             | (maxgeom) Bucket capacity *B* used in MaxGeomHash. |
| `--w W`             | integer       | 64             | (maxgeom, alphamaxgeom) Maximum number of buckets *W*. |
| `--alpha A`         | float         | 0.45           | (alphamaxgeom) Alpha parameter for α-MaxGeomHash. |
| `--scale S`         | float         | 0.001          | (fracminhash) Scale parameter controlling sampling probability. |
| `--num-perm K`      | integer       | 1000           | (minhash) Number of permutations for classical MinHash. |
| `--seed SEED`       | integer       | 42             | Random seed for reproducibility. |
| `--canonical`       | flag          | false          | Treat each k-mer as canonical (min of forward/reverse complement). |
| `--keep-ambiguous`  | flag          | false          | Keep k-mers containing ambiguous bases instead of skipping. |
| `--output OUT`      | string (path) | **required**   | Output sketch file path. |

## `pwsimilarity`

This program computes pairwise similarity from a list of sketch files. The arguments are as follows:

| Argument               | Argument Type      | Default Value | What It Means |
|------------------------|--------------------|----------------|----------------|
| `--metric METRIC`      | string             | jaccard        | Similarity metric to compute between sketches. Options: `jaccard`, `cosine`. |
| `--output OUT.tsv`     | string (path)      | pairs.tsv      | Output TSV file containing pairwise similarity results. |
| `SKETCH1 SKETCH2 ...`  | list of file paths | **required**   | Input sketch files to compare pairwise. |


## `expt_growth`

Run repeated-seed experiments at growing set sizes. Jaccard and cosine retain
the existing equal-size set synthesis. Containment is directional:
`C(A,B) = |A intersect B| / |A|`, not overlap divided by the smaller set size.

```bash
bin/expt_growth --t 0.5 --metric containment --seeds 500 --steps 10 \
  --growth x2 --out results/fixed_containment_expt_amgh_t0.5_a0.45 \
  --algo alphamaxgeom --alpha 0.45 --base_n 100000 --size_multiplier 40
```

For step `s` (starting at zero), `|A| = base_n * growth_factor^s`,
`|B| = round(|A| * size_multiplier)`, and the intersection contains
`round(t * |A|)` distinct elements. Thus the example starts with 100,000 and
4,000,000 elements, then 200,000 and 8,000,000, and so on. `true_sim` records
the realized containment after integer rounding; MSE is measured against this
value, not the unrounded target. `--growth` accepts `x2` or `x10`.

`--size_multiplier` defaults to 1 and accepts positive finite numbers. Values
below 1 are allowed when the target fits inside B (`t <= size_multiplier`).
Empty rounded B sizes and overflowing growth schedules are rejected. Non-unit
multipliers require `--metric containment`.

Containment supports `alphamaxgeom`, `maxgeom`, `bottomk`, and `fracminhash`.
Use `--k` for MaxGeom bucket capacity or bottom-k sketch size, `--alpha` for
AlphaMaxGeom, and `--scale` for FracMinHash. For example, replace the algorithm
options in the command above with `--algo maxgeom --k 50`,
`--algo bottomk --k 2000`, or `--algo fracminhash --scale 0.001`.
Classical `--algo minhash` is not supported for containment; it fails explicitly
rather than silently returning a different metric. This section describes
`expt_growth`; the `pwsimilarity` CLI remains unchanged.

Containment output keeps the existing TSV columns and appends `valid_seeds`.
All requested trials run, but a trial with no A elements in the coordinated
sample has an undefined containment estimate. Such trials are **not replaced
with 0 or 1**: `mean_est` and `mse` use only the defined estimates, and a warning
reports the excluded count. Both statistics are `nan` if there are no valid
trials. For an MSE over all 500 trials, ensure `valid_seeds` is 500. In particular,
the default bottom-k `--k 50` can be too sparse at a 40-fold size imbalance;
increase `--k` before comparing MSEs across algorithms. The two sample-size
columns average the stored sketch sizes over all trials, not just valid trials.

Containment synthesis streams distinct integer IDs through a seeded 64-bit
permutation directly into the actual sketch implementations; shared elements
receive identical hashes. It does not allocate a universal string pool or
materialize A and B. The seed controls reproducibility. Runtime still includes
all sketch updates in every trial, and FracMinHash memory still grows with its
sampling rate and input size; a 500-seed, ten-step, 40-fold experiment is large.
Existing Jaccard/cosine synthesis and TSV columns are unchanged.

Run the dependency-free containment unit and CLI smoke tests with:

```bash
python3 test/test_containment.py
```

## Citing

*This work was accepted and presented at RECOMB 2025 and is now under review at Genome Research.*

Please cite the following preprint if you use MaxGeomHash or α-MaxGeomHash.
```
Hera, M. Rahman, Koslicki, D., & Martínez, C. (2025). MaxGeomHash: An algorithm for variable-size random sampling of distinct elements [Preprint]. bioRxiv. https://doi.org/10.1101/2025.11.11.687920
```