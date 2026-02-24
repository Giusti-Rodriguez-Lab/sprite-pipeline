use crate::config::{Config, TagCategory};
use crate::fastq::Record;
use crate::tag::TagMaps;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Process a read pair against the config layouts.
///
/// The combined barcode label is built by processing R1 (layout1) then R2
/// (layout2) and appending each tag result — exactly as the Java version does
/// by sharing a single `StringBuilder sb` across both reads.  Both output
/// records get the same combined label appended to their names.
pub fn process_pair(
    r1:       &Record,
    r2:       &Record,
    config:   &Config,
    tag_maps: &TagMaps,
) -> (Record, Record) {
    let mut label: Vec<u8> = Vec::with_capacity(128);

    check_read(&r1.seq, &config.layout1, tag_maps, config.spacer_len, config.laxity, &mut label);
    check_read(&r2.seq, &config.layout2, tag_maps, config.spacer_len, config.laxity, &mut label);

    let out_r1 = r1.with_barcode_label(&label);
    let out_r2 = r2.with_barcode_label(&label);
    (out_r1, out_r2)
}

// ---------------------------------------------------------------------------
// Per-read processing
// ---------------------------------------------------------------------------

/// Scan `seq` according to `layout`, appending barcode results to `label`.
///
/// Mirrors Java's `checkRead()`:
/// - SPACER elements skip `spacer_len` bytes.
/// - Tag elements run the sliding-window hash lookup.
/// - The very first layout element (whether spacer or tag) is the "start of
///   read": if it is a tag position, zero laxity is used (maxAdvance = 0).
///   All subsequent tag positions use `max_advance = laxity`.
pub fn check_read(
    seq:        &[u8],
    layout:     &[TagCategory],
    tag_maps:   &TagMaps,
    spacer_len: usize,
    laxity:     usize,
    label:      &mut Vec<u8>,
) {
    let mut bases: &[u8] = seq;
    let mut start_of_read = true;

    for &category in layout {
        if category == TagCategory::Spacer {
            bases = skip_spacer(bases, spacer_len);
        } else {
            let effective_advance = if start_of_read { 0 } else { laxity };
            bases = match_tag_and_advance(bases, category, tag_maps, effective_advance, label);
        }
        start_of_read = false;
    }
}

// ---------------------------------------------------------------------------
// Inner sliding-window matching
// ---------------------------------------------------------------------------

/// Try to find a tag from `category` in `bases`, sliding up to `max_advance`
/// positions from the current cursor.
///
/// Appends `[TagName]` on success or `[NOT_FOUND]` on failure to `label`.
/// Returns the remaining bytes after consuming the tag (or the fallback
/// advancement on NOT_FOUND).
///
/// Mirrors Java's `checkReadForTagAndReturnRemainder`:
/// ```text
/// int advance = -1;
/// while (match == null && advance < maxAdvance) {
///     advance++;                    // 0..=maxAdvance
///     for each length: hash lookup
/// }
/// ```
fn match_tag_and_advance<'a>(
    bases:       &'a [u8],
    category:    TagCategory,
    tag_maps:    &TagMaps,
    max_advance: usize,
    label:       &mut Vec<u8>,
) -> &'a [u8] {
    let lengths = tag_maps.lengths(category);

    let mut advance: usize = 0;
    loop {
        for &len in lengths {
            let end = advance + len;
            if end <= bases.len() {
                let candidate = &bases[advance..end];
                if let Some(name) = tag_maps.lookup(category, candidate) {
                    label.push(b'[');
                    label.extend_from_slice(name.as_bytes());
                    label.push(b']');
                    return &bases[end..];
                }
            }
        }

        if advance >= max_advance {
            break;
        }
        advance += 1;
    }

    // --- NOT_FOUND ---
    // Advance the cursor by (advance + last_length), using the largest known
    // length for this category.  Matches Java's fallback behaviour.
    label.extend_from_slice(b"[NOT_FOUND]");
    let last_len = lengths.last().copied().unwrap_or(0);
    let skip     = advance + last_len;
    if skip >= bases.len() {
        &bases[bases.len()..]
    } else {
        &bases[skip..]
    }
}

/// Skip exactly `spacer_len` bases.  If fewer bases remain, consume all.
fn skip_spacer(bases: &[u8], spacer_len: usize) -> &[u8] {
    if spacer_len >= bases.len() {
        &bases[bases.len()..]
    } else {
        &bases[spacer_len..]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::tag::TagMaps;

    const MINI_CONFIG: &str = "\
READ1 = DPM
READ2 = Y|SPACER|ODD
SPACER = 3
LAXITY = 2
DPM\tDPM_A\tAAAACCCC\t0
DPM\tDPM_B\tGGGGTTTT\t0
Y\tY_1\tTTTTAAAA\t0
Y\tY_2\tCCCCGGGG\t0
ODD\tOdd_1\tATCGATCG\t1
ODD\tOdd_2\tGCTAGCTA\t1
";

    fn setup() -> (Config, TagMaps) {
        let cfg  = Config::from_str_content(MINI_CONFIG).unwrap();
        let maps = TagMaps::build(&cfg);
        (cfg, maps)
    }

    // ---- skip_spacer ----

    #[test]
    fn test_spacer_skip_normal() {
        let bases: &[u8] = b"NNNNACGT";
        assert_eq!(skip_spacer(bases, 4), b"ACGT");
    }

    #[test]
    fn test_spacer_skip_exact_length() {
        let bases: &[u8] = b"NNNN";
        assert_eq!(skip_spacer(bases, 4), b"");
    }

    #[test]
    fn test_spacer_skip_longer_than_read() {
        let bases: &[u8] = b"NN";
        assert_eq!(skip_spacer(bases, 6), b"");
    }

    #[test]
    fn test_spacer_skip_zero() {
        let bases: &[u8] = b"ACGT";
        assert_eq!(skip_spacer(bases, 0), b"ACGT");
    }

    // ---- match_tag_and_advance ----

    #[test]
    fn test_exact_match_at_advance0() {
        let (_, maps) = setup();
        let mut label = Vec::new();
        let rest = match_tag_and_advance(
            b"AAAACCCCXXXXXXXX",
            TagCategory::Dpm,
            &maps,
            0,   // max_advance = 0 (start of read)
            &mut label,
        );
        assert_eq!(label, b"[DPM_A]");
        assert_eq!(rest,  b"XXXXXXXX");
    }

    #[test]
    fn test_not_found_zero_laxity() {
        // At start-of-read (maxAdvance=0), a shifted tag is NOT found.
        let (_, maps) = setup();
        let mut label = Vec::new();
        // DPM_A starts at position 2 — should NOT be found with maxAdvance=0
        let rest = match_tag_and_advance(
            b"XXAAAACCCCYYY",
            TagCategory::Dpm,
            &maps,
            0,
            &mut label,
        );
        assert_eq!(label, b"[NOT_FOUND]");
        // cursor = advance(0) + last_len(8) = 8 → &bases[8..] = "CCCYYY"... wait
        // bases = XXAAAACCCCYYY (13 chars)
        // skip = 0 + 8 = 8 → &bases[8..] = "CCYYY" ... let me recount
        // X X A A A A C C C C Y Y Y
        // 0 1 2 3 4 5 6 7 8 9 ...
        // skip=8 → &bases[8..] = "CCYYY"
        assert_eq!(rest, b"CCYYY");
    }

    #[test]
    fn test_match_found_with_laxity() {
        // Tag shifted by 2 bases; found with maxAdvance=2
        let (_, maps) = setup();
        let mut label = Vec::new();
        let rest = match_tag_and_advance(
            b"XXAAAACCCC",
            TagCategory::Dpm,
            &maps,
            2,   // max_advance = 2 (laxity)
            &mut label,
        );
        assert_eq!(label, b"[DPM_A]");
        assert_eq!(rest,  b"");  // consumed XX + AAAACCCC = all 10 bytes
    }

    #[test]
    fn test_hamming1_neighbor_found() {
        // ATCGATCC is Hamming-1 from Odd_1=ATCGATCG → should resolve to Odd_1
        let (_, maps) = setup();
        let mut label = Vec::new();
        let rest = match_tag_and_advance(
            b"ATCGATCC",
            TagCategory::Odd,
            &maps,
            0,
            &mut label,
        );
        assert_eq!(label, b"[Odd_1]");
        assert_eq!(rest,  b"");
    }

    #[test]
    fn test_not_found_advances_by_longest_len() {
        // When no match at any advance position, cursor = advance + longest_len
        let (_, maps) = setup();
        let mut label = Vec::new();
        // maxAdvance=2, no match anywhere → advance ends at 2
        // longest DPM length = 8 → skip = 2 + 8 = 10
        let rest = match_tag_and_advance(
            b"XXXXXXXXXXREST",
            TagCategory::Dpm,
            &maps,
            2,
            &mut label,
        );
        assert_eq!(label, b"[NOT_FOUND]");
        assert_eq!(rest,  b"REST");
    }

    #[test]
    fn test_empty_bases_returns_not_found() {
        let (_, maps) = setup();
        let mut label = Vec::new();
        let rest = match_tag_and_advance(
            b"",
            TagCategory::Dpm,
            &maps,
            0,
            &mut label,
        );
        assert_eq!(label, b"[NOT_FOUND]");
        assert_eq!(rest,  b"");
    }

    // ---- check_read (full layout) ----

    #[test]
    fn test_check_read_layout1_perfect_match() {
        // READ1 = DPM; start-of-read → zero laxity
        let (cfg, maps) = setup();
        let mut label = Vec::new();
        check_read(
            b"AAAACCCCXXXX",
            &cfg.layout1,
            &maps,
            cfg.spacer_len,
            cfg.laxity,
            &mut label,
        );
        assert_eq!(label, b"[DPM_A]");
    }

    #[test]
    fn test_check_read_layout2_full() {
        // READ2 = Y|SPACER|ODD
        // R2 = TTTTAAAA (Y_1, 8bp) + NNN (spacer=3) + ATCGATCG (Odd_1, exact)
        let (cfg, maps) = setup();
        let mut label = Vec::new();
        check_read(
            b"TTTTAAAANNNATCGATCG",
            &cfg.layout2,
            &maps,
            cfg.spacer_len,
            cfg.laxity,
            &mut label,
        );
        assert_eq!(label, b"[Y_1][Odd_1]");
    }

    #[test]
    fn test_check_read_layout2_odd_with_laxity() {
        // ODD is not start-of-read → can use laxity
        // After Y+SPACER, ODD_1 is shifted by 2 within the available laxity
        // R2 = TTTTAAAA (Y_1) + NNN (spacer=3) + XX (offset 2) + ATCGATCG (Odd_1)
        let (cfg, maps) = setup();
        let mut label = Vec::new();
        check_read(
            b"TTTTAAAAXXXATCGATCG",  // spacer=NNN(3)+offset=XX = 5 "wasted" bytes before ODD
            // Actually: after Y(8) + spacer(3) = 11 bytes consumed → remaining = "XXATCGATCG"
            // Wait: TTTTAAAA(8) + NNN(3) consumed from "TTTTAAAAXXXATCGATCG"... let me recount
            // Total: T T T T A A A A X X X A T C G A T C G = 19 chars
            //        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
            // After Y(8 consumed): bases = XXXATCGATCG (11 chars)
            // After SPACER(3): bases = ATCGATCG... wait: XXX = 3 → ATCGATCG (8 chars)
            // After SPACER(3): remaining = ATCGATCG — exact match! No laxity needed here.
            &cfg.layout2,
            &maps,
            cfg.spacer_len,
            cfg.laxity,
            &mut label,
        );
        assert_eq!(label, b"[Y_1][Odd_1]");
    }

    #[test]
    fn test_check_read_layout2_odd_laxity_needed() {
        // After Y + spacer, ODD shifted by 2 extra bytes (needs laxity=2)
        // R2 = TTTTAAAA + NNN + XX + ATCGATCG  (spacer=3, then 2-byte offset before ODD)
        let (cfg, maps) = setup();
        let mut label = Vec::new();
        // Build: Y_1(8) + spacer_bytes(3) + offset(2) + Odd_1(8) = 21 bytes
        let seq = b"TTTTAAAANNNXXATCGATCG";
        //          Y_1        spc off Odd_1
        // After Y(8): bases = NNNXXATCGATCG (13)
        // After spacer(3): bases = XXATCGATCG (10)
        // ODD: advance=0 → XXATCGAT (no match), advance=1 → XATCGATC (no match),
        //      advance=2 → ATCGATCG (Odd_1 exact match!)
        check_read(seq, &cfg.layout2, &maps, cfg.spacer_len, cfg.laxity, &mut label);
        assert_eq!(label, b"[Y_1][Odd_1]");
    }

    // ---- process_pair ----

    #[test]
    fn test_process_pair_both_reads_get_combined_label() {
        let (cfg, maps) = setup();
        // R1: DPM_A perfect match
        // R2: Y_1 + NNN(spacer) + Odd_1 perfect match
        let r1 = crate::fastq::Record {
            header: b"@read1".to_vec(),
            seq:    b"AAAACCCCXXXXXX".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIIIIIIIIIIIII".to_vec(),
        };
        let r2 = crate::fastq::Record {
            header: b"@read1".to_vec(),
            seq:    b"TTTTAAAANNNATCGATCG".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIIIIIIIIIIIIIIIIII".to_vec(),
        };
        let (out1, out2) = process_pair(&r1, &r2, &cfg, &maps);
        assert_eq!(out1.header, b"@read1::[DPM_A][Y_1][Odd_1]");
        assert_eq!(out2.header, b"@read1::[DPM_A][Y_1][Odd_1]");
        // Sequence and quality pass through unchanged
        assert_eq!(out1.seq, b"AAAACCCCXXXXXX");
        assert_eq!(out2.seq, b"TTTTAAAANNNATCGATCG");
    }

    #[test]
    fn test_process_pair_not_found_in_r1() {
        let (cfg, maps) = setup();
        let r1 = crate::fastq::Record {
            header: b"@read2".to_vec(),
            seq:    b"XXXXXXXXXXXXXXXXXX".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIIIIIIIIIIIIIIIII".to_vec(),
        };
        let r2 = crate::fastq::Record {
            header: b"@read2".to_vec(),
            seq:    b"TTTTAAAANNNATCGATCG".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIIIIIIIIIIIIIIIIII".to_vec(),
        };
        let (out1, out2) = process_pair(&r1, &r2, &cfg, &maps);
        let expected = b"@read2::[NOT_FOUND][Y_1][Odd_1]";
        assert_eq!(out1.header, expected);
        assert_eq!(out2.header, expected);
    }

    #[test]
    fn test_process_pair_empty_layout1() {
        // Config with only READ2 — R1 contributes nothing to the label
        let cfg_str = "READ2 = Y\nY\tY_1\tTTTTAAAA\t0\n";
        let cfg      = Config::from_str_content(cfg_str).unwrap();
        let maps     = TagMaps::build(&cfg);

        let r1 = crate::fastq::Record {
            header: b"@r".to_vec(),
            seq:    b"ANYTHING".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIIIIIII".to_vec(),
        };
        let r2 = crate::fastq::Record {
            header: b"@r".to_vec(),
            seq:    b"TTTTAAAA".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIIIIIII".to_vec(),
        };
        let (out1, out2) = process_pair(&r1, &r2, &cfg, &maps);
        assert_eq!(out1.header, b"@r::[Y_1]");
        assert_eq!(out2.header, b"@r::[Y_1]");
    }
}
