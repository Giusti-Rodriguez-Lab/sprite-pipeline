#!/usr/bin/env bash
# compare_java.sh
#
# Validates that the Rust barcode_id binary produces bit-identical output
# (post-decompression) to the Java BarcodeIdentification_v1.2.0.jar on the
# same input FASTQ files.
#
# Usage:
#   ./tests/compare_java.sh \
#       --jar     /path/to/BarcodeIdentification_v1.2.0.jar \
#       --config  /path/to/barcode_config.txt \
#       --input1  sample_R1.fq.gz \
#       --input2  sample_R2.fq.gz
#
# Must be run after `cargo build --release` so the binary is at
# target/release/barcode_id.
#
# Exit codes: 0 = pass, 1 = diff found or error.

set -euo pipefail

JAR=""
CONFIG=""
INPUT1=""
INPUT2=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --jar)     JAR="$2";    shift 2 ;;
        --config)  CONFIG="$2"; shift 2 ;;
        --input1)  INPUT1="$2"; shift 2 ;;
        --input2)  INPUT2="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$JAR" || -z "$CONFIG" || -z "$INPUT1" || -z "$INPUT2" ]]; then
    echo "Usage: $0 --jar JAR --config CONFIG --input1 R1.fq.gz --input2 R2.fq.gz"
    exit 1
fi

BINARY="$(dirname "$0")/../target/release/barcode_id"
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: Binary not found at $BINARY. Run 'cargo build --release' first."
    exit 1
fi

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== Running Java ==="
time java -jar "$JAR" \
    --input1  "$INPUT1" \
    --input2  "$INPUT2" \
    --output1 "$TMPDIR/java_R1.fq.gz" \
    --output2 "$TMPDIR/java_R2.fq.gz" \
    --config  "$CONFIG"

echo ""
echo "=== Running Rust ==="
time "$BINARY" \
    --input1  "$INPUT1" \
    --input2  "$INPUT2" \
    --output1 "$TMPDIR/rust_R1.fq.gz" \
    --output2 "$TMPDIR/rust_R2.fq.gz" \
    --config  "$CONFIG"

echo ""
echo "=== Comparing R1 ==="
diff <(zcat "$TMPDIR/java_R1.fq.gz") <(zcat "$TMPDIR/rust_R1.fq.gz") \
    && echo "R1: IDENTICAL" || { echo "R1: DIFF FOUND"; exit 1; }

echo "=== Comparing R2 ==="
diff <(zcat "$TMPDIR/java_R2.fq.gz") <(zcat "$TMPDIR/rust_R2.fq.gz") \
    && echo "R2: IDENTICAL" || { echo "R2: DIFF FOUND"; exit 1; }

# Read-count conservation
INPUT_COUNT=$(zcat "$INPUT1" | awk 'NR%4==1' | wc -l | tr -d ' ')
RUST_COUNT=$(zcat "$TMPDIR/rust_R1.fq.gz" | awk 'NR%4==1' | wc -l | tr -d ' ')
if [[ "$INPUT_COUNT" -ne "$RUST_COUNT" ]]; then
    echo "ERROR: Read count mismatch! Input=$INPUT_COUNT Rust=$RUST_COUNT"
    exit 1
fi
echo "Read count conserved: $RUST_COUNT pairs"

echo ""
echo "=== ALL CHECKS PASSED ==="
