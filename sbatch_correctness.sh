#!/bin/bash
#SBATCH --job-name=barcode_correctness
#SBATCH --cpus-per-task=16
#SBATCH -o barcode_correctness.%j.out
#SBATCH -e barcode_correctness.%j.err
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --qos=giustirodriguezp-b

module load java
module load rust/1.87.0

PROJ=/blue/giustirodriguezp/PROJECTS/Sprite/Guttman_Ismagilov_Labs/KK_scSPRITE/scSPRITE
cd "$PROJ"

# Build the Rust release binary for this node's architecture (Linux x86_64).
# Must rebuild here — the macOS binary from your laptop won't run on the cluster.
cd barcode_id && cargo build --release && cd ..

# ---------------------------------------------------------------------------
# Pick a FASTQ pair.
# Edit these two lines to point at a specific sample if needed.
# ---------------------------------------------------------------------------
FASTQ_DIR=/blue/giustirodriguezp/PROJECTS/Sprite/Guttman_Ismagilov_Labs/scSPRITE/raw_fastq
INPUT1=$(ls "$FASTQ_DIR"/*_R1*.fastq.gz "$FASTQ_DIR"/*_R1*.fq.gz 2>/dev/null | head -1)
INPUT2=$(ls "$FASTQ_DIR"/*_R2*.fastq.gz "$FASTQ_DIR"/*_R2*.fq.gz 2>/dev/null | head -1)

if [[ -z "$INPUT1" || -z "$INPUT2" ]]; then
    echo "ERROR: Could not find R1/R2 FASTQ files in $FASTQ_DIR" >&2
    exit 1
fi

echo "Using R1: $INPUT1"
echo "Using R2: $INPUT2"

bash barcode_id/tests/compare_java.sh \
    --jar    misc/sprite-pipeline/java/BarcodeIdentification_v1.2.0.jar \
    --config misc/config_dpm6_y-stag_scSPRITE2.txt \
    --input1 "$INPUT1" \
    --input2 "$INPUT2"
