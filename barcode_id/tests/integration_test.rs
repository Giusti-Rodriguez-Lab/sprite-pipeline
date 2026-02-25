//! Integration tests: run the barcode_id binary end-to-end on synthetic FASTQ
//! data and verify the output matches manually-computed expected values.
//!
//! These tests are the "correctness gate" described in the spec (T9).  They
//! do not require the Java JAR.  A separate shell script
//! `tests/compare_java.sh` handles the Java-vs-Rust diff test (run on the
//! HPC where the JAR is available).

use std::io::{Read, Write};
use std::path::Path;
use std::process::Command;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Write bytes to a gzipped file.
fn write_gz(path: &Path, data: &[u8]) {
    let file    = std::fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(file, Compression::default());
    enc.write_all(data).unwrap();
    enc.finish().unwrap();
}

/// Read a gzipped file and return its decompressed bytes.
fn read_gz(path: &Path) -> Vec<u8> {
    let file = std::fs::File::open(path).unwrap();
    let mut dec = GzDecoder::new(file);
    let mut buf = Vec::new();
    dec.read_to_end(&mut buf).unwrap();
    buf
}

/// Run the barcode_id binary and return exit status.
fn run_binary(
    config:  &Path,
    input1:  &Path,
    input2:  &Path,
    output1: &Path,
    output2: &Path,
) -> std::process::ExitStatus {
    run_binary_threaded(config, input1, input2, output1, output2, 1).0
}

/// Run with an explicit --threads count; returns (exit status, captured stderr).
fn run_binary_threaded(
    config:  &Path,
    input1:  &Path,
    input2:  &Path,
    output1: &Path,
    output2: &Path,
    threads: usize,
) -> (std::process::ExitStatus, String) {
    let binary = env!("CARGO_BIN_EXE_barcode_id");
    let out = Command::new(binary)
        .env("RUST_LOG", "info")
        .arg("--config")  .arg(config)
        .arg("--input1")  .arg(input1)
        .arg("--input2")  .arg(input2)
        .arg("--output1") .arg(output1)
        .arg("--output2") .arg(output2)
        .arg("--threads") .arg(threads.to_string())
        .output()
        .expect("Failed to launch barcode_id binary");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (out.status, stderr)
}

/// FASTQ fixture: build a minimal 4-line FASTQ string from name + sequence.
/// The quality string is the same length as the sequence, filled with 'I'.
fn make_fastq(pairs: &[(&str, &str)]) -> Vec<u8> {
    let mut out = Vec::new();
    for &(name, seq) in pairs {
        let qual: String = "I".repeat(seq.len());
        let line = format!("@{}\n{}\n+\n{}\n", name, seq, qual);
        out.extend_from_slice(line.as_bytes());
    }
    out
}

/// Parse decompressed FASTQ bytes into (header, seq) pairs.
fn parse_headers(data: &[u8]) -> Vec<String> {
    let text = std::str::from_utf8(data).unwrap();
    text.lines()
        .enumerate()
        .filter(|(i, _)| i % 4 == 0)
        .map(|(_, line)| line.to_string())
        .collect()
}

fn config_path() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/mini_config.txt")
}

// ---------------------------------------------------------------------------
// Barcode identification correctness tests
// ---------------------------------------------------------------------------

/// Pair 1 — all perfect matches.
///
/// Layout:
///   READ1 = DPM  (start-of-read, zero laxity)
///   READ2 = Y | SPACER(3) | ODD  (ODD has laxity=2)
///
/// R1: DPM_A(8) + padding
/// R2: Y_1(8) + NNN(spacer) + Odd_1(8)
///
/// Expected label: [DPM_A][Y_1][Odd_1]
#[test]
fn test_pair_all_perfect_matches() {
    let tmp = TempDir::new().unwrap();

    // R1: AAAACCCC (DPM_A, 8bp) + 10 bp padding
    let r1_data = make_fastq(&[("read1", "AAAACCCCXXXXXXXXXX")]);
    // R2: TTTTAAAA (Y_1, 8bp) + NNN (spacer=3) + ATCGATCG (Odd_1, 8bp)
    let r2_data = make_fastq(&[("read1", "TTTTAAAANNNNATCGATCG")]);
    // Note: spacer=3, so after consuming Y_1 (8bp), we skip 3bp → remaining is NATCGATCG
    // ODD advance=0: NATCGATC no match; advance=1: ATCGATCG → Odd_1 ✓

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    let h2 = parse_headers(&read_gz(&out2));
    // Both reads get the same combined label
    assert_eq!(h1[0], "@read1::[DPM_A][Y_1][Odd_1]");
    assert_eq!(h2[0], "@read1::[DPM_A][Y_1][Odd_1]");
}

/// Pair 2 — NOT_FOUND in R1, Hamming-1 match in ODD.
///
/// R1: unknown DPM → [NOT_FOUND]
/// R2: Y_1 + NNN(spacer) + ATCGATCC (Odd_1 with 1 mismatch at position 7)
///
/// Expected label: [NOT_FOUND][Y_1][Odd_1]
#[test]
fn test_pair_not_found_r1_mismatch_odd() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[("read2", "XXXXXXXXXXXXXXXXXX")]);
    // ATCGATCC is Hamming-1 from Odd_1=ATCGATCG (G→C at pos 7)
    let r2_data = make_fastq(&[("read2", "TTTTAAAANNNNATCGATCC")]);

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    let h2 = parse_headers(&read_gz(&out2));
    assert_eq!(h1[0], "@read2::[NOT_FOUND][Y_1][Odd_1]");
    assert_eq!(h2[0], "@read2::[NOT_FOUND][Y_1][Odd_1]");
}

/// Pair 3 — ODD found via laxity (shifted 2 bp after the spacer).
///
/// After Y_1 and spacer(3), the ODD tag sits 2 bases late — needs advance=2.
/// laxity=2 in the config so it should still be found.
///
/// Expected label: [DPM_B][Y_1][Odd_1]
#[test]
fn test_pair_odd_found_via_laxity() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[("read3", "GGGGTTTTXXXXXXXXXX")]);
    // Y_1(8) + spacer(3) + 2-byte offset + Odd_1(8) = 8+3+2+8 = 21 bytes
    let r2_data = make_fastq(&[("read3", "TTTTAAAANNNXXATCGATCG")]);
    //                                    Y_1     spc off Odd_1

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    assert_eq!(h1[0], "@read3::[DPM_B][Y_1][Odd_1]");
}

/// Pair 4 — ODD exceeds laxity limit: shift=3 with laxity=2 → NOT_FOUND.
///
/// Expected label: [DPM_A][Y_1][NOT_FOUND]
#[test]
fn test_pair_odd_exceeds_laxity() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[("read4", "AAAACCCCXXXXXXXXXX")]);
    // Y_1(8) + spacer(3) + 3-byte offset (exceeds laxity=2) + Odd_1(8) = 22 bytes
    let r2_data = make_fastq(&[("read4", "TTTTAAAANNNXXXATCGATCG")]);

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    assert_eq!(h1[0], "@read4::[DPM_A][Y_1][NOT_FOUND]");
}

/// Pair 5 — DPM has zero laxity (start-of-read): a shifted DPM must NOT match.
///
/// DPM_A starts at position 2 in R1, but since it is the first tag (zero
/// laxity), it must not be found.
///
/// Expected label: [NOT_FOUND][Y_1][Odd_1]
#[test]
fn test_pair_dpm_shift_rejected_at_start_of_read() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[("read5", "XXAAAACCCCXXXXXXXX")]);
    //                                    ^^ 2-byte offset — must not match (zero laxity)
    let r2_data = make_fastq(&[("read5", "TTTTAAAANNNNATCGATCG")]);

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    assert_eq!(h1[0], "@read5::[NOT_FOUND][Y_1][Odd_1]");
}

/// Pair 6 — read shorter than the first tag: graceful NOT_FOUND.
#[test]
fn test_pair_read_too_short() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[("read6", "AAA")]);  // only 3 bytes; DPM needs 8
    let r2_data = make_fastq(&[("read6", "TTT")]);  // only 3 bytes; Y needs 7+

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    // Both tags NOT_FOUND, all positions
    assert!(h1[0].contains("[NOT_FOUND]"));
    assert!(!h1[0].contains("[DPM_A]"));
}

/// Pair 7 — multiple read pairs in the same file; count must be conserved.
#[test]
fn test_multiple_pairs_count_conserved() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[
        ("p1", "AAAACCCCXXXXXXXXXX"),
        ("p2", "GGGGTTTTXXXXXXXXXX"),
        ("p3", "XXXXXXXXXXXXXXXXXX"),
    ]);
    let r2_data = make_fastq(&[
        ("p1", "TTTTAAAANNNNATCGATCG"),
        ("p2", "CCCCGGGNNNATCGATCG"),
        ("p3", "TTTTAAAANNNNATCGATCG"),
    ]);

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    let status = run_binary(&config_path(), &in1, &in2, &out1, &out2);
    assert!(status.success());

    let h1 = parse_headers(&read_gz(&out1));
    let h2 = parse_headers(&read_gz(&out2));
    assert_eq!(h1.len(), 3, "output1 should have 3 records");
    assert_eq!(h2.len(), 3, "output2 should have 3 records");
}

/// Pair 8 — sequence is passed through unchanged (only header is modified).
#[test]
fn test_sequence_and_qual_passthrough() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[("r", "AAAACCCCXXXXXXXXXX")]);
    let r2_data = make_fastq(&[("r", "TTTTAAAANNNNATCGATCG")]);

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    run_binary(&config_path(), &in1, &in2, &out1, &out2);

    let out_text = String::from_utf8(read_gz(&out1)).unwrap();
    let lines: Vec<&str> = out_text.lines().collect();
    // Line 1 (index 1) is sequence
    assert_eq!(lines[1], "AAAACCCCXXXXXXXXXX");
    // Line 3 (index 3) is quality
    assert_eq!(lines[3], "I".repeat(18));
}

/// Pair 9 — header format: every output header must match the expected pattern.
#[test]
fn test_header_format_regex() {
    let tmp = TempDir::new().unwrap();

    let r1_data = make_fastq(&[
        ("SRR123.1 lane=1", "AAAACCCCXXXXXXXXXX"),  // has comment after space
        ("SRR123.2",         "GGGGTTTTXXXXXXXXXX"),
    ]);
    let r2_data = make_fastq(&[
        ("SRR123.1 lane=1", "TTTTAAAANNNNATCGATCG"),
        ("SRR123.2",         "CCCCGGGNNNATCGATCG"),
    ]);

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &r1_data);
    write_gz(&in2, &r2_data);

    run_binary(&config_path(), &in1, &in2, &out1, &out2);

    let headers = parse_headers(&read_gz(&out1));

    // Comment ("lane=1") must be stripped — only name before first space is kept.
    assert!(headers[0].starts_with("@SRR123.1::"),
        "Expected @SRR123.1::..., got {}", headers[0]);
    assert!(!headers[0].contains("lane"),
        "Comment should be stripped: {}", headers[0]);

    // All headers must match pattern: @name::[tags]
    for h in &headers {
        assert!(h.starts_with('@'), "Header must start with @: {}", h);
        assert!(h.contains("::"),  "Header must contain '::': {}", h);
        let after = h.split("::").nth(1).unwrap();
        // Must contain at least one barcode token
        assert!(
            after.contains('[') && after.contains(']'),
            "Header barcode section must have [...] tokens: {}", h,
        );
    }
}

// ---------------------------------------------------------------------------
// Parallelization tests
// ---------------------------------------------------------------------------

/// Pair 10 — output order is preserved under multi-threaded execution.
///
/// Generate 500 pairs with sequential names (read_0000..read_0499).  Run with
/// --threads 4.  If rayon's par_iter+collect ordering guarantee is broken the
/// headers will come back in an arbitrary order; this test catches that.
#[test]
fn test_parallel_output_order_preserved() {
    let tmp = TempDir::new().unwrap();
    let n: usize = 500;

    let names: Vec<String> = (0..n).map(|i| format!("read_{:04}", i)).collect();

    let r1_pairs: Vec<(&str, &str)> = names
        .iter()
        .map(|name| (name.as_str(), "AAAACCCCXXXXXXXXXX"))
        .collect();
    let r2_pairs: Vec<(&str, &str)> = names
        .iter()
        .map(|name| (name.as_str(), "TTTTAAAANNNNATCGATCG"))
        .collect();

    let in1  = tmp.path().join("r1.fq.gz");
    let in2  = tmp.path().join("r2.fq.gz");
    let out1 = tmp.path().join("out1.fq.gz");
    let out2 = tmp.path().join("out2.fq.gz");
    write_gz(&in1, &make_fastq(&r1_pairs));
    write_gz(&in2, &make_fastq(&r2_pairs));

    let (status, stderr) = run_binary_threaded(&config_path(), &in1, &in2, &out1, &out2, 4);
    assert!(status.success());

    // Verify rayon actually configured 4 threads (not silently falling back to 1)
    assert!(
        stderr.contains("Using 4 rayon threads."),
        "Expected 'Using 4 rayon threads.' in stderr, got:\n{}", stderr,
    );

    let headers1 = parse_headers(&read_gz(&out1));
    let headers2 = parse_headers(&read_gz(&out2));

    assert_eq!(headers1.len(), n, "output1 record count mismatch");
    assert_eq!(headers2.len(), n, "output2 record count mismatch");

    for (i, (h1, h2)) in headers1.iter().zip(headers2.iter()).enumerate() {
        let expected_name = format!("@read_{:04}::", i);
        assert!(
            h1.starts_with(&expected_name),
            "output1 record {} out of order: expected prefix '{}', got '{}'",
            i, expected_name, h1,
        );
        assert_eq!(h1, h2, "output1 and output2 headers differ at record {}", i);
    }
}

/// Pair 11 — results are identical regardless of thread count.
///
/// Run the same input through --threads 1 and --threads 4, then compare every
/// output header.  Verifies that parallelism does not change correctness.
#[test]
fn test_parallel_matches_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let n: usize = 200;

    let names: Vec<String> = (0..n).map(|i| format!("r{}", i)).collect();
    let r1_pairs: Vec<(&str, &str)> = names
        .iter()
        .map(|name| (name.as_str(), "AAAACCCCXXXXXXXXXX"))
        .collect();
    let r2_pairs: Vec<(&str, &str)> = names
        .iter()
        .map(|name| (name.as_str(), "TTTTAAAANNNNATCGATCG"))
        .collect();

    let in1 = tmp.path().join("r1.fq.gz");
    let in2 = tmp.path().join("r2.fq.gz");
    write_gz(&in1, &make_fastq(&r1_pairs));
    write_gz(&in2, &make_fastq(&r2_pairs));

    // Single-threaded run
    let out1_single = tmp.path().join("out1_single.fq.gz");
    let out2_single = tmp.path().join("out2_single.fq.gz");
    let (s1, stderr1) = run_binary_threaded(&config_path(), &in1, &in2, &out1_single, &out2_single, 1);
    assert!(s1.success());
    assert!(stderr1.contains("Using 1 rayon threads."),
        "Expected 'Using 1 rayon threads.' in stderr, got:\n{}", stderr1);

    // Multi-threaded run
    let out1_multi = tmp.path().join("out1_multi.fq.gz");
    let out2_multi = tmp.path().join("out2_multi.fq.gz");
    let (s4, stderr4) = run_binary_threaded(&config_path(), &in1, &in2, &out1_multi, &out2_multi, 4);
    assert!(s4.success());
    assert!(stderr4.contains("Using 4 rayon threads."),
        "Expected 'Using 4 rayon threads.' in stderr, got:\n{}", stderr4);

    let h_single = parse_headers(&read_gz(&out1_single));
    let h_multi  = parse_headers(&read_gz(&out1_multi));

    assert_eq!(h_single.len(), h_multi.len(), "record count differs between 1-thread and 4-thread runs");
    for (i, (hs, hm)) in h_single.iter().zip(h_multi.iter()).enumerate() {
        assert_eq!(hs, hm, "header mismatch at record {} between 1-thread and 4-thread runs", i);
    }
}
