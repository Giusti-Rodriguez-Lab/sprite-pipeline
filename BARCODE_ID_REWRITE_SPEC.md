# BarcodeIdentification Rewrite Specification
## Rust replacement for `BarcodeIdentification_v1.2.0.jar`

**Version:** 0.1-draft
**Status:** Draft
**Context:** scSPRITE pipeline, `run_barcode` rule in Snakefile

---

## 1. Background and Motivation

### 1.1 What BarcodeIdentification Does in the Pipeline

The SPRITE (Split-Pool Recognition of Interactions by Tag Extension) protocol physically tags RNA and DNA molecules with combinatorial barcodes before sequencing. Each sequenced read encodes a barcode "address" that identifies which spatial cluster the molecule came from.

`BarcodeIdentification` is **Step 1** of the pipeline: it reads raw paired-end FASTQ files and annotates each read's name with the barcodes it found. Every downstream step — ligation efficiency, alignment, BAM filtering, clustering, contact maps — depends entirely on the output of this step. It is the single most I/O- and CPU-intensive step in the pipeline.

### 1.2 Why the Java Implementation Is a Bottleneck

The current implementation (`BarcodeIdentification_v1.2.0.jar`, Guttman Lab v1.2.0) has several fundamental performance limitations:

**Single-threaded processing.** The Java code processes one read pair at a time in a single thread. A typical scSPRITE experiment produces 500M–1B read pairs. On a 32-core HPC node, 31 cores sit idle during this step.

**String allocation in the hot inner loop.** For every tag position on every read, the code calls `String.substring()` — which allocates a new heap object — up to `(maxAdvance × numTagLengths)` times. At 6 laxity positions and multiple tag lengths, this is hundreds of allocations per read pair, or tens of billions of short-lived objects over a full run. This drives continuous minor GC.

**JVM GC pauses.** All those short-lived String allocations trigger frequent stop-the-world minor collections. On long runs this adds up to minutes of cumulative pause time and unpredictable latency spikes.

**Sequential GZIP I/O.** Java's `GZIPInputStream`/`GZIPOutputStream` decompress and compress on a single thread, even though the input FASTQs are independent and could be streamed in parallel.

**JVM startup overhead.** Each sample in the Snakefile spawns a new JVM process. With many samples this is 200–500 ms of startup overhead per sample, plus JIT warm-up time before the code reaches full throughput.

**Quantified impact.** In a 16-sample scSPRITE run (~800M total read pairs), the `run_barcode` rule consumes roughly 60–80% of total wall-clock time. Reducing this step from 8 hours to 1–2 hours would make iterating on pipeline parameters (barcode configs, laxity settings) practical rather than prohibitive.

### 1.3 Why Rust

The workload is byte-string matching on streaming data: compute-bound, latency-sensitive, memory-allocation-dominated. Rust provides:

- **Zero-copy byte slices** (`&[u8]`): the hot-path hash lookup can operate directly on a slice of the input buffer with no allocation, vs. Java's mandatory `String` heap object.
- **`rayon` data parallelism**: process batches of read pairs across all cores with a single `.par_iter()` — near-linear scaling to core count.
- **`AHashMap`**: non-cryptographic hash map with ~2x throughput over Java's `HashMap` for short byte-key lookups.
- **No GC**: deterministic memory behavior, no pauses.
- **`flate2`/`bgzf`**: GZIP with optional multi-threaded compression via `bgzf` (blocked gzip, compatible with standard `gzip -d`).
- **Single static binary**: no JVM, no classpath, no module conflicts. Drop in and run.

Expected speedup: **5–15x** wall-clock on a 16-core node (based on similar bioinformatics tool rewrites: fastp vs. Trimmomatic, minimap2 vs. BWA-MEM).

---

## 2. Scope

### In Scope
- Exact functional replication of `BarcodeIdentification_v1.2.0.jar`
- Paired-end FASTQ input (gzipped)
- Gzipped FASTQ output
- Same config file format (tab-delimited, as documented in `example_config.txt`)
- Same output format: barcode string appended to read name as `readname::[Tag1][Tag2]...[NOT_FOUND]...`
- Multi-core parallelism
- Drop-in replacement in the existing Snakefile (`run_barcode` rule)

### Out of Scope
- Changes to downstream pipeline steps
- New barcode matching algorithms (e.g., edit-distance beyond Hamming)
- Support for single-end reads (not used in SPRITE)
- SAM/BAM input (handled by `filter_all_tags.py` later in the pipeline)
- GUI or web interface

---

## 3. Functional Requirements

### FR-1: Config File Parsing

The tool must parse the same config format as the Java version:

```
READ1 = <CATEGORY>[|<CATEGORY>...]
READ2 = <CATEGORY>[|<CATEGORY>...]
SPACER = <integer>          # optional, default 6
LAXITY = <integer>          # optional, default 6
<CATEGORY>\t<NAME>\t<SEQUENCE>\t<NUM_MISMATCHES>
```

- Valid categories: `EVEN`, `ODD`, `Y`, `RPM`, `DPM`, `LIGTAG`, `SPACER`
- Lines beginning with `#` or empty lines are ignored
- If `SPACER` appears twice, warn and use last value (match Java behavior)
- If `LAXITY` appears twice, warn and use last value
- Exit with a clear error message if no layouts or no tags are found
- Warn (but do not exit) if two tags within the same category are within `mismatches` Hamming distance of each other (ambiguity warning)

### FR-2: Hamming Neighbor Pre-computation

At startup, for each tag with `NUM_MISMATCHES = k > 0`, generate all sequences within Hamming distance `k` and insert them into the lookup map pointing back to the canonical tag name. The alphabet is `{A, C, G, T, N}`.

- This matches the Java `Tag.generateTagsWithinHammingOf()` logic exactly
- For `k=0`, only the exact sequence is inserted
- Maps are per-category: `HashMap<category, HashMap<sequence, canonical_tag_name>>`

### FR-3: Read Processing

For each read pair (R1, R2):

1. Build a single barcode label string by processing R1 then R2 according to their respective layouts.
2. For each element in the layout, left to right:
   - If `SPACER`: skip exactly `spacerLength` bases. If fewer bases remain, consume all remaining.
   - Otherwise (a tag category): call the tag-matching procedure.
3. Tag-matching procedure:
   - Starting at `advance = 0`, up to `advance = maxAdvance` (inclusive):
     - For each tag length known for this category (ascending order):
       - Extract `bases[advance..advance+len]`
       - Hash-lookup in the category's map
       - If found: append `[canonical_name]` to the label, advance the cursor by `advance + len`, stop searching
     - If not found at this advance, increment advance and retry
   - If no match at any advance position: append `[NOT_FOUND]`, advance cursor by `advance + lengths[last]` (i.e., using the last/longest tag length at the final advance position)
4. Append the label to **both** R1 and R2 read names: `<original_name>::<label>`
   - The `<original_name>` is the part before the first space in the FASTQ header line
5. Write the modified R1 and R2 records to their respective output FASTQ files

### FR-4: FASTQ I/O

- Input: gzipped FASTQ (`.fastq.gz` / `.fq.gz`)
- Output: gzipped FASTQ, standard gzip format (readable by `zcat`, `gzip -d`, Python `gzip`, etc.)
- Preserve all 4 lines of each FASTQ record (header, sequence, `+`, quality)
- Only the header line (line 1) is modified
- Quality scores and sequence are passed through unchanged

### FR-5: CLI Interface

The binary must accept (at minimum) these flags, matching the Java JAR interface:

```
--input1   <path>   Input FASTQ R1 (gzipped)
--input2   <path>   Input FASTQ R2 (gzipped)
--output1  <path>   Output FASTQ R1 (gzipped)
--output2  <path>   Output FASTQ R2 (gzipped)
--config   <path>   Config file
--threads  <int>    Number of worker threads (default: all available cores)
--version          Print version and exit
--help             Print usage and exit
```

The Snakefile invocation must work unchanged after substituting the binary path for the JAR:

```bash
# Before (Java):
java -jar BarcodeIdentification_v1.2.0.jar --input1 R1.fq.gz --input2 R2.fq.gz \
    --output1 R1.barcoded.fastq.gz --output2 R2.barcoded.fastq.gz --config config.txt

# After (Rust):
barcode_id --input1 R1.fq.gz --input2 R2.fq.gz \
    --output1 R1.barcoded.fastq.gz --output2 R2.barcoded.fastq.gz --config config.txt
```

### FR-6: Progress Logging

Log to stderr (not stdout):
- Startup: config file path, number of tags loaded per category, spacer length, laxity
- Progress: every 500,000 read pairs processed (matching Java's log cadence)
- Completion: total read pairs processed, wall-clock elapsed time

---

## 4. Non-Functional Requirements

### NFR-1: Correctness (Highest Priority)

The output FASTQs must be **bit-for-bit identical** to the Java version on the same input, modulo GZIP metadata (timestamps, OS byte). The decompressed content must be identical. This is the gate condition for all other requirements.

### NFR-2: Performance

| Metric | Target | Rationale |
|--------|--------|-----------|
| Throughput (single-core) | ≥ 2× Java single-core | Eliminate allocation overhead alone |
| Throughput (16-core) | ≥ 8× Java single-core | Near-linear scaling to ~50% efficiency |
| Throughput (32-core) | ≥ 14× Java single-core | Target for HPC nodes |
| Memory (peak) | < 4 GB | Fit within typical HPC job allocations |
| Startup time | < 50 ms | vs. ~300 ms JVM cold start |
| Binary size | < 20 MB | Static binary, no runtime dep |

### NFR-3: Compatibility

- Output must be accepted by all downstream pipeline steps without modification:
  - `get_ligation_efficiency.py` (parses `[TagName]` regex from FASTQ headers)
  - `filter_all_tags.py` (searches for `NOT_FOUND` string in read names)
  - STAR/Parabricks alignment (read names must be valid SAM query names after `::` appended)
- Config files must be parsed without modification to existing `.txt` config files

### NFR-4: Reliability

- Must not silently drop or corrupt read pairs
- Read pair count in output must equal read pair count in input
- Must exit with non-zero status on I/O errors, malformed config, or unreadable input

### NFR-5: Build and Deployment

- Single static binary (no shared library dependencies)
- Builds with `cargo build --release` on Linux x86_64 (HPC target) and macOS arm64/x86_64 (dev)
- No external runtime dependencies (no JVM, no Python, no conda)
- Can be placed anywhere on `$PATH`; no install step required

### NFR-6: Maintainability

- Code structured into clear modules: `config`, `tag`, `fastq`, `match`, `io`
- Inline doc comments on public functions
- All logic unit-testable without disk I/O

---

## 5. Validation Strategy

### 5.1 Validation Philosophy

The Rust implementation is a **correctness-first rewrite**. The validation strategy follows three phases:

1. **Unit correctness**: Each component produces correct output independently
2. **Integration correctness**: Full pipeline output matches Java on known inputs
3. **Performance validation**: Throughput meets NFR-2 targets

### 5.2 Test Data Sets

Three test datasets should be used:

| Dataset | Description | Read Pairs | Purpose |
|---------|-------------|-----------|---------|
| **Tiny** | 1,000 synthetic read pairs, all barcodes perfect matches | ~1000 | Fast iteration, debugging |
| **Small** | 100,000 real reads from a published scSPRITE run | ~100K | Real-world distribution |
| **Full** | One complete sample from the scSPRITE run | ~50–100M | Performance benchmarking |

The **Tiny** and **Small** datasets must have Java-generated ground truth (run the JAR, save output) committed alongside the test inputs.

### 5.3 Correctness Test Suite

#### T-1: Exact sequence match against Java output

```bash
# Run Java on small dataset
java -jar BarcodeIdentification_v1.2.0.jar \
  --input1 test/small_R1.fq.gz --input2 test/small_R2.fq.gz \
  --output1 java_out_R1.fq.gz --output2 java_out_R2.fq.gz \
  --config test/config.txt

# Run Rust on same input
barcode_id \
  --input1 test/small_R1.fq.gz --input2 test/small_R2.fq.gz \
  --output1 rust_out_R1.fq.gz --output2 rust_out_R2.fq.gz \
  --config test/config.txt

# Decompress and diff (ignoring GZIP metadata)
diff <(zcat java_out_R1.fq.gz) <(zcat rust_out_R1.fq.gz)
diff <(zcat java_out_R2.fq.gz) <(zcat rust_out_R2.fq.gz)
# Both diffs must produce zero output
```

Pass criteria: **zero diff** on both R1 and R2 across all datasets.

#### T-2: Read pair count conservation

```bash
# Input count
input_count=$(zcat test/small_R1.fq.gz | wc -l | awk '{print $1/4}')
# Output count
rust_count=$(zcat rust_out_R1.fq.gz | wc -l | awk '{print $1/4}')
[ "$input_count" -eq "$rust_count" ] || echo "FAIL: read count mismatch"
```

Pass criteria: output read count equals input read count, for all datasets.

#### T-3: Header format validation

Every output read name must match the pattern:
```
^@[^ ]+::[(\[NOT_FOUND\]|\[[a-zA-Z0-9_\-]+\])+$
```
i.e., original name (up to first space), `::`, then one or more barcode tokens.

```bash
zcat rust_out_R1.fq.gz | awk 'NR%4==1' | \
  grep -vP '^@[^ ]+::(\[NOT_FOUND\]|\[[a-zA-Z0-9_\-]+\])+$' | wc -l
# Must be 0
```

#### T-4: NOT_FOUND rate within expected range

On real data, the NOT_FOUND rate for each tag position should fall within a biologically expected range. If >95% of reads show NOT_FOUND for all positions, the matching logic is likely broken.

```bash
# Use existing get_ligation_efficiency.py on Rust output:
python get_ligation_efficiency.py rust_out_R1.fq.gz
# Compare distribution to Java output — must be statistically identical
```

#### T-5: Ambiguity warning parity

On a config with known ambiguous barcodes (two tags within Hamming distance of each other), both Java and Rust must emit a warning. Construct a synthetic config for this.

#### T-6: Edge cases

| Test | Input | Expected Output |
|------|-------|----------------|
| Read shorter than first tag | Short read (< tag length) | `[NOT_FOUND]` for all positions |
| Read exhausted mid-layout | Read shorter than full layout | `[NOT_FOUND]` for remaining positions |
| Perfect match at advance=0 | Exact barcode at start | `[TagName]` found, cursor advanced |
| Match only at advance=maxAdvance | Tag shifted by 6 bp | `[TagName]` found at max laxity |
| Tag with k=0 mismatches | Exact match required | Only exact sequence matches |
| Tag with k=2 mismatches | 2 substitutions in barcode | Still matches to canonical name |
| Config with only READ2 layout | No READ1 layout | R1 gets empty `::` label, R2 gets full label |
| LAXITY=0 | No laxity | Only checks advance=0 |
| Empty reads (length 0) | Zero-length sequence | All `[NOT_FOUND]` |

### 5.4 Unit Tests (in Rust `#[cfg(test)]`)

#### Unit Test Suite

| Test | What It Tests |
|------|--------------|
| `test_hamming_k0` | `generate_neighbors(seq, 0)` returns exactly `{seq}` |
| `test_hamming_k1_len5` | A 5-mer with k=1 returns exactly 25 neighbors (5 positions × 5 bases) |
| `test_hamming_k2_spot_check` | Known sequence at Hamming-2 is present in neighbors |
| `test_hamming_n_included` | `N` is one of the substitution bases (matches Java behavior) |
| `test_config_parse_basic` | Round-trip parse of `example_config.txt` produces correct tag counts |
| `test_config_spacer_default` | Config without `SPACER` line defaults to 6 |
| `test_config_laxity_default` | Config without `LAXITY` line defaults to 6 |
| `test_match_exact` | Exact-match tag at advance=0 found immediately |
| `test_match_with_laxity` | Tag at advance=3 found when maxAdvance=6 |
| `test_match_not_found` | Unknown sequence returns NOT_FOUND and advances by longest tag |
| `test_spacer_skip` | SPACER skips exactly N bases |
| `test_spacer_short_read` | SPACER on read shorter than spacer length returns empty |
| `test_full_read2_layout` | Full `Y|SPACER|EVEN|SPACER|ODD|SPACER|EVEN|SPACER|ODD|SPACER|DPM` produces correct label |
| `test_header_append` | `readname::barcodes` format is correct, space-trimmed name |
| `test_fastq_roundtrip` | Parse 4-line FASTQ block, modify name, serialize back to identical format |

### 5.5 Performance Benchmarks

Run on the **Full** dataset (one complete sample, ~50M read pairs) on the HPC node (32-core, 256 GB RAM):

```bash
# Java baseline (single-threaded, as it currently runs)
time java -jar BarcodeIdentification_v1.2.0.jar \
  --input1 full_R1.fq.gz --input2 full_R2.fq.gz \
  --output1 /dev/null --output2 /dev/null \
  --config config.txt
# Record: wall time, CPU time, peak RSS

# Rust single-threaded
time barcode_id --threads 1 \
  --input1 full_R1.fq.gz --input2 full_R2.fq.gz \
  --output1 /dev/null --output2 /dev/null \
  --config config.txt

# Rust multi-threaded
for t in 4 8 16 32; do
  time barcode_id --threads $t \
    --input1 full_R1.fq.gz --input2 full_R2.fq.gz \
    --output1 /dev/null --output2 /dev/null \
    --config config.txt
done
```

Record: reads/second, wall time, CPU efficiency (wall × threads / CPU time).

Pass criteria:
- Rust 1-thread ≥ 2× Java wall time
- Rust 16-thread ≥ 8× Java wall time
- Rust 32-thread ≥ 14× Java wall time

---

## 6. Implementation Tasks and Dependencies

```
Task dependency graph:

T1 ─────────────────────────────────────────────────────┐
T2 ────────────────────────────────────────┐            │
T3 ──────────┐                             │            │
             ▼                             ▼            ▼
T4 ──────── T5 ──────── T6 ──────── T7 ──────── T8 ──── T9
             │                      │
             ▼                      ▼
           T10                    T11
```

### T1 — Project scaffold
**Depends on:** nothing
**Description:** Initialize a Rust binary crate (`cargo new barcode_id`). Add dependencies to `Cargo.toml`: `flate2`, `ahash`, `rayon`, `clap` (v4), `thiserror`, `log`, `env_logger`. Set up module structure: `main.rs`, `config.rs`, `tag.rs`, `fastq.rs`, `matcher.rs`, `io.rs`. Add CI (GitHub Actions) with `cargo test` + `cargo clippy`.

### T2 — Test fixtures
**Depends on:** nothing (parallel with T1)
**Description:** Create `tests/fixtures/`. Generate the **Tiny** synthetic dataset (1,000 reads with known barcodes at known positions, mix of exact matches, 1-mismatch, 2-mismatch, NOT_FOUND). Run the Java JAR on this dataset and on `example_config.txt`. Commit: inputs + Java ground-truth outputs + a mini config file. These fixtures are the acceptance gate for all integration tests.

### T3 — Config parser (`config.rs`)
**Depends on:** T1
**Description:** Implement parsing of the tab-delimited config format. Parse `READ1`/`READ2` layouts into `Vec<TagCategory>`. Parse `SPACER`, `LAXITY`. Parse tag lines into `Vec<TagDef>` (category, name, sequence, mismatches). Expose `Config` struct. Unit tests: T-`test_config_*` from §5.4. No I/O in tests (pass config as string, not path).

### T4 — Hamming neighbor generation (`tag.rs`)
**Depends on:** T3
**Description:** Implement `generate_neighbors(seq: &[u8], k: u8) -> HashSet<Vec<u8>>` matching Java's `generateTagsWithinHammingOf()`. Alphabet = `[A, C, G, T, N]`. Build the per-category lookup map: `HashMap<TagCategory, AHashMap<Vec<u8>, &str>>`. Unit tests: `test_hamming_*` from §5.4. Verify against Java: a known sequence at Hamming-2 resolves to the canonical name.

### T5 — Core matching logic (`matcher.rs`)
**Depends on:** T3, T4
**Description:** Implement `check_read(bases: &[u8], layout: &[TagCategory], tag_maps: &TagMaps, spacer_len: usize, max_advance: usize) -> String`. This is the direct port of Java's `checkRead()` + `checkReadForTagAndReturnRemainder()`. Must use zero-copy `&bases[start..end]` slices for hash lookups — no `String` allocation in the inner loop. Unit tests: all `test_match_*` and `test_spacer_*` from §5.4.

### T6 — FASTQ parser/serializer (`fastq.rs`)
**Depends on:** T1
**Description:** Implement a streaming FASTQ record parser. A `FastqRecord` holds: `name: Vec<u8>`, `seq: Vec<u8>`, `plus: Vec<u8>`, `qual: Vec<u8>`. Implement parsing from a `BufRead` source and serialization back to 4-line text. Handle gzip decompression/compression via `flate2`. Unit test: `test_fastq_roundtrip` from §5.4. Do not load the full file into memory; process record by record.

### T7 — Single-threaded integration (`main.rs`)
**Depends on:** T3, T4, T5, T6
**Description:** Wire together: parse args (clap), load config, build tag maps, open input FASTQs, iterate record pairs, call `check_read` for R1 and R2, append label to both headers, write output FASTQs. Log progress every 500K pairs. Integration test: run against Tiny fixtures, compare decompressed output to Java ground truth with `diff`. This is the correctness gate — no parallelism yet.

### T8 — Parallel processing (`main.rs` + `io.rs`)
**Depends on:** T7
**Description:** Add `--threads` flag. Implement a producer/consumer pipeline: one reader thread that batches read pairs into chunks (e.g., 10,000 pairs per chunk), a rayon thread pool that processes chunks in parallel (each chunk calls `check_read` on each pair), one writer thread that collects results in order and writes to output. Order preservation is mandatory — output record order must match input order. Benchmark against Java baseline.

### T9 — Integration test harness
**Depends on:** T7, T2
**Description:** Write `tests/integration_test.rs` (or a shell script `tests/compare_java.sh`) that: (1) runs the Rust binary on Tiny and Small fixtures, (2) decompresses both Rust and Java outputs, (3) `diff`s them, (4) checks read counts, (5) checks header format with regex. This test runs in CI on every PR. Also implement T-4 (NOT_FOUND rate check) using a Python helper that calls `get_ligation_efficiency.py` on both outputs and asserts distributions are identical.

### T10 — Performance benchmarks
**Depends on:** T8
**Description:** Add a `benches/` directory with `criterion` benchmarks for: (a) `generate_neighbors` for a single 15-mer at k=2, (b) `check_read` for a single read with the full scSPRITE2 layout, (c) end-to-end throughput on the Tiny dataset. Run full benchmarks (T-§5.5) on the HPC node against the Java baseline. Record and commit results to `benchmarks/results.md`.

### T11 — Snakefile integration and documentation
**Depends on:** T9
**Description:** Update `Snakefile` `run_barcode` rule to accept either the JAR (for backward compat) or the Rust binary based on a `config.yaml` flag (e.g., `RUNBC_RUST: true`). Update `README.md` with build instructions (`cargo build --release`), placement of binary, and performance notes. Add a `--dry-run` / diff mode to the binary itself that processes N reads and prints what labels would have been assigned (useful for debugging config files).

---

## 7. Task Summary Table

| Task | Description | Depends on | Estimated Complexity |
|------|-------------|-----------|---------------------|
| T1 | Project scaffold + CI | — | Small |
| T2 | Test fixtures (Java ground truth) | — | Small |
| T3 | Config parser | T1 | Medium |
| T4 | Hamming neighbor generation | T3 | Medium |
| T5 | Core matching logic | T3, T4 | Medium |
| T6 | FASTQ parser/serializer | T1 | Medium |
| T7 | Single-threaded integration | T3, T4, T5, T6 | Medium |
| T8 | Parallel processing | T7 | Large |
| T9 | Integration test harness | T7, T2 | Medium |
| T10 | Performance benchmarks | T8 | Small |
| T11 | Snakefile integration + docs | T9 | Small |

**Critical path:** T1 → T3 → T4 → T5 → T7 → T9 (correctness gate) → T8 → T10

T2 and T6 can be developed in parallel with T3/T4/T5.
T10 and T11 can be developed in parallel after T9 passes.

---

## 8. Acceptance Criteria

The rewrite is considered complete when all of the following pass:

- [ ] `diff <(zcat java_out_R1.fq.gz) <(zcat rust_out_R1.fq.gz)` produces zero output on **all three** datasets
- [ ] `diff <(zcat java_out_R2.fq.gz) <(zcat rust_out_R2.fq.gz)` produces zero output on **all three** datasets
- [ ] Read count is conserved across all datasets
- [ ] Header format regex passes across all datasets
- [ ] All unit tests pass (`cargo test`)
- [ ] All integration tests pass in CI
- [ ] Rust 16-thread throughput ≥ 8× Java single-thread on the Full dataset
- [ ] `cargo clippy` produces zero warnings
- [ ] Binary builds successfully on Linux x86_64 (HPC) and macOS arm64 (dev)
- [ ] Snakefile `run_barcode` rule works with new binary on a real sample end-to-end
