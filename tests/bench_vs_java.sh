#!/usr/bin/env bash
# bench_vs_java.sh — wall-clock comparison: Rust barcode_id vs Java JAR
#
# Usage (from project root):
#   chmod +x tests/bench_vs_java.sh
#   ./tests/bench_vs_java.sh [--n-reads N] [--threads T]
#
# Outputs a Markdown results table appended to benchmarks/results.md.

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
N=500000    # read pairs
T=0         # threads for Rust multi (0 = all CPUs)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --n-reads) N="$2"; shift 2 ;;
        --threads) T="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths (relative to project root; script lives in tests/)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

JAR="$ROOT/misc/sprite-pipeline/java/BarcodeIdentification_v1.2.0.jar"
BIN="$ROOT/barcode_id/target/release/barcode_id"
CFG="$ROOT/misc/config_dpm6_y-stag_scSPRITE2.txt"
RESULTS="$ROOT/benchmarks/results.md"
GENERATOR="$SCRIPT_DIR/gen_test_fastq.py"

# ---------------------------------------------------------------------------
# Prerequisite checks
# ---------------------------------------------------------------------------
command -v java    >/dev/null 2>&1 || { echo "ERROR: 'java' not found in PATH" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' not found in PATH" >&2; exit 1; }
command -v gzip    >/dev/null 2>&1 || { echo "ERROR: 'gzip' not found in PATH" >&2; exit 1; }

[[ -f "$JAR"       ]] || { echo "ERROR: JAR not found: $JAR" >&2; exit 1; }
[[ -f "$CFG"       ]] || { echo "ERROR: config not found: $CFG" >&2; exit 1; }
[[ -f "$GENERATOR" ]] || { echo "ERROR: generator not found: $GENERATOR" >&2; exit 1; }

if [[ ! -f "$BIN" ]]; then
    echo "ERROR: Rust binary not found: $BIN" >&2
    echo "  Build it first:  cd barcode_id && cargo build --release" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Temporary directory (cleaned up on exit)
# ---------------------------------------------------------------------------
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

R1="$TMP/R1.fq.gz"
R2="$TMP/R2.fq.gz"
JAVA_O1="$TMP/java_out1.fq.gz";  JAVA_O2="$TMP/java_out2.fq.gz"
RUST1_O1="$TMP/rust1_out1.fq.gz"; RUST1_O2="$TMP/rust1_out2.fq.gz"
RUSTN_O1="$TMP/rustn_out1.fq.gz"; RUSTN_O2="$TMP/rustn_out2.fq.gz"

# ---------------------------------------------------------------------------
# Generate synthetic FASTQ data
# ---------------------------------------------------------------------------
echo "=== Generating $N synthetic read pairs ==="
python3 "$GENERATOR" \
    --config  "$CFG" \
    --out-r1  "$R1" \
    --out-r2  "$R2" \
    --n-reads "$N"

# ---------------------------------------------------------------------------
# Timing helper
#
# TIMEFORMAT='%R' makes bash's time builtin print only the real elapsed
# seconds to the group's stderr.  Redirecting that stderr to stdout with
# } 2>&1, inside $(), captures the time string.  The command's own
# stdout/stderr are suppressed by >/dev/null 2>&1 inside the group.
# ---------------------------------------------------------------------------
TIMEFORMAT='%R'

# ---------------------------------------------------------------------------
# Java JAR
# ---------------------------------------------------------------------------
echo ""
echo "=== Java JAR ==="
JAVA_TIME=$({ time java -jar "$JAR" \
    --input1 "$R1"    --input2 "$R2" \
    --output1 "$JAVA_O1" --output2 "$JAVA_O2" \
    --config "$CFG" >/dev/null 2>&1; } 2>&1)
echo "  wall: ${JAVA_TIME}s"

# ---------------------------------------------------------------------------
# Rust — 1 thread
# ---------------------------------------------------------------------------
echo ""
echo "=== Rust (1 thread) ==="
RUST1_TIME=$({ time "$BIN" \
    --input1 "$R1"     --input2 "$R2" \
    --output1 "$RUST1_O1" --output2 "$RUST1_O2" \
    --config "$CFG" --threads 1 >/dev/null 2>&1; } 2>&1)
echo "  wall: ${RUST1_TIME}s"

# ---------------------------------------------------------------------------
# Resolve thread count for multi-threaded run
# ---------------------------------------------------------------------------
if [[ "$T" -eq 0 ]]; then
    T=$(python3 -c "import os; print(os.cpu_count())")
fi

# ---------------------------------------------------------------------------
# Rust — N threads
# ---------------------------------------------------------------------------
echo ""
echo "=== Rust ($T threads) ==="
RUSTN_TIME=$({ time "$BIN" \
    --input1 "$R1"     --input2 "$R2" \
    --output1 "$RUSTN_O1" --output2 "$RUSTN_O2" \
    --config "$CFG" --threads "$T" >/dev/null 2>&1; } 2>&1)
echo "  wall: ${RUSTN_TIME}s"

# ---------------------------------------------------------------------------
# Correctness check: Java output vs Rust single-thread output (R1 file)
# ---------------------------------------------------------------------------
echo ""
echo "=== Correctness check (Java vs Rust 1-thread, R1 output) ==="
if diff <(gzip -dc "$JAVA_O1") <(gzip -dc "$RUST1_O1") >/dev/null; then
    CORRECTNESS="PASS"
    echo "  PASS: outputs are identical"
else
    CORRECTNESS="FAIL"
    echo "  FAIL: outputs differ — run diff manually to inspect" >&2
fi

# ---------------------------------------------------------------------------
# Throughput and speedup (computed in Python to handle float arithmetic)
# ---------------------------------------------------------------------------
JAVA_RPS=$( python3 -c "print(int($N / $JAVA_TIME))")
RUST1_RPS=$(python3 -c "print(int($N / $RUST1_TIME))")
RUSTN_RPS=$(python3 -c "print(int($N / $RUSTN_TIME))")

JAVA_KPS=$( python3 -c "print(f'{$JAVA_RPS  / 1000:.1f}')")
RUST1_KPS=$(python3 -c "print(f'{$RUST1_RPS / 1000:.1f}')")
RUSTN_KPS=$(python3 -c "print(f'{$RUSTN_RPS / 1000:.1f}')")

RUST1_X=$(python3 -c "print(f'{$JAVA_TIME / $RUST1_TIME:.2f}')")
RUSTN_X=$(python3 -c "print(f'{$JAVA_TIME / $RUSTN_TIME:.2f}')")

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
DATE=$(date +%Y-%m-%d)
OS=$(uname -s)
ARCH=$(uname -m)

echo ""
echo "=== Summary ==="
printf "  %-20s  %8ss  %8sK reads/sec\n" "Java JAR"          "$JAVA_TIME"  "$JAVA_KPS"
printf "  %-20s  %8ss  %8sK reads/sec  %s× vs Java\n" \
    "Rust (1 thread)"  "$RUST1_TIME" "$RUST1_KPS" "$RUST1_X"
printf "  %-20s  %8ss  %8sK reads/sec  %s× vs Java\n" \
    "Rust ($T threads)" "$RUSTN_TIME" "$RUSTN_KPS" "$RUSTN_X"
echo "  Correctness: $CORRECTNESS"

# ---------------------------------------------------------------------------
# Append Markdown table to benchmarks/results.md
# ---------------------------------------------------------------------------
mkdir -p "$(dirname "$RESULTS")"

{
echo ""
echo "## Run: $DATE — $OS/$ARCH / ${N} reads / correctness: $CORRECTNESS"
echo ""
echo "| Runner              | Wall time    | Reads/sec     | Speedup vs Java |"
echo "|---------------------|--------------|---------------|-----------------|"
printf "| %-19s | %-12s | %-13s | %-15s |\n" \
    "Java JAR"           "${JAVA_TIME}s"  "${JAVA_KPS}K/s"  "1.00×"
printf "| %-19s | %-12s | %-13s | %-15s |\n" \
    "Rust (1 thread)"    "${RUST1_TIME}s" "${RUST1_KPS}K/s" "${RUST1_X}×"
printf "| %-19s | %-12s | %-13s | %-15s |\n" \
    "Rust ($T threads)"  "${RUSTN_TIME}s" "${RUSTN_KPS}K/s" "${RUSTN_X}×"
} >> "$RESULTS"

echo ""
echo "Results appended to: $RESULTS"
