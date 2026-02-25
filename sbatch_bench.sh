#!/bin/bash
#SBATCH --job-name=barcode_bench
#SBATCH --cpus-per-task=16
#SBATCH -o barcode_bench.%j.out
#SBATCH -e barcode_bench.%j.err
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --qos=giustirodriguezp-b

module load java
module load rust/1.87.0

PROJ=/blue/giustirodriguezp/PROJECTS/Sprite/Guttman_Ismagilov_Labs/KK_scSPRITE/scSPRITE
cd "$PROJ"

# Build the Rust release binary for this node's architecture (Linux x86_64).
cd barcode_id && cargo build --release && cd ..

# Run the benchmark: generates 5M synthetic read pairs, times Java and Rust
# (1 thread + 16 threads), checks correctness, appends results to benchmarks/results.md.
bash tests/bench_vs_java.sh --n-reads 5000000 --threads 16
