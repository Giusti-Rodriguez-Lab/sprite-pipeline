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
    let binary = env!("CARGO_BIN_EXE_barcode_id");
    Command::new(binary)
        .arg("--config")  .arg(config)
        .arg("--input1")  .arg(input1)
        .arg("--input2")  .arg(input2)
        .arg("--output1") .arg(output1)
        .arg("--output2") .arg(output2)
        .status()
        .expect("Failed to launch barcode_id binary")
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
