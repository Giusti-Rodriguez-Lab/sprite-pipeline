use std::collections::{BTreeSet, HashMap, HashSet};

use ahash::AHashMap;

use crate::config::{Config, TagCategory};

// ---------------------------------------------------------------------------
// Hamming neighbor generation (mirrors Java Tag.generateTagsWithinHammingOf)
// ---------------------------------------------------------------------------

/// Alphabet used for substitution — must match Java: {A, C, G, T, N}.
const BASES: &[u8] = b"ACGTN";

/// Generate all sequences within Hamming distance `k` from `seq`.
///
/// The alphabet is {A, C, G, T, N}. This mirrors the Java implementation
/// which replaces every position with each of the five bases (including the
/// original), producing a set of all sequences at Hamming distance 0..=k.
pub fn generate_neighbors(seq: &[u8], k: u8) -> HashSet<Vec<u8>> {
    let mut current: HashSet<Vec<u8>> = HashSet::new();
    current.insert(seq.to_vec());

    for _ in 0..k {
        current = expand_by_one(current);
    }

    current
}

/// Given a set of sequences, return all sequences at Hamming distance 0 or 1
/// from any sequence in the set.  This is the inner loop of the Java
/// `generateTagsWithinHammingOfOne` method.
fn expand_by_one(seqs: HashSet<Vec<u8>>) -> HashSet<Vec<u8>> {
    let mut result: HashSet<Vec<u8>> = seqs.clone();
    for seq in &seqs {
        for i in 0..seq.len() {
            for &base in BASES {
                let mut variant = seq.clone();
                variant[i] = base;
                result.insert(variant);
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// 2-bit packed encoding for fast exact-match lookup
// ---------------------------------------------------------------------------

/// Encode `seq` as a packed `u64` using 2 bits per base (A=00, C=01, G=10, T=11),
/// most-significant pair = first base.  Returns `None` if any base is not ACGT
/// (e.g. N) or if the sequence is longer than 32 bp.
pub fn encode_packed(seq: &[u8]) -> Option<u64> {
    if seq.len() > 32 {
        return None;
    }
    let mut packed: u64 = 0;
    for &b in seq {
        let bits: u64 = match b {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            _    => return None, // N or any non-ACGT base
        };
        packed = (packed << 2) | bits;
    }
    Some(packed)
}

/// Count the number of mismatched positions between two packed sequences given
/// the XOR of their encodings and the sequence length.
///
/// Uses a 2-bit-per-position collapse: any position where either bit differs
/// contributes 1 to the count.
pub fn count_mismatches_packed(xor: u64, len: usize) -> u32 {
    // Mask off any bits above the encoded length (upper bits are 0 from encode_packed,
    // but guard anyway).
    let mask = if len >= 32 { u64::MAX } else { (1u64 << (2 * len)) - 1 };
    let xor  = xor & mask;
    // OR the two bits of each 2-bit position together, then count set bits.
    let collapsed = (xor | (xor >> 1)) & 0x5555_5555_5555_5555_u64;
    collapsed.count_ones()
}

// ---------------------------------------------------------------------------
// Packed tag tables (exact-match fast path)
// ---------------------------------------------------------------------------

struct PackedTag {
    packed: u64,
    len:    usize,  // length of the original sequence; guards against leading-A collisions
    name:   String,
}

struct PackedCategory {
    tags: Vec<PackedTag>,
}

// ---------------------------------------------------------------------------
// Tag lookup maps
// ---------------------------------------------------------------------------

/// Pre-computed per-category lookup tables built from the config file.
///
/// For each tag category:
///   - `maps[category]` maps a sequence (bytes) to the canonical tag name.
///     Hamming-neighbor sequences also map to the canonical name.
///   - `lengths[category]` is a sorted list of unique tag lengths for that
///     category (used by the sliding-window matcher to try each length).
///   - `packed[category]` holds canonical sequences encoded as u64 for an
///     O(n) exact-match fast path that avoids heap allocation.
pub struct TagMaps {
    maps:    HashMap<TagCategory, AHashMap<Vec<u8>, String>>,
    lengths: HashMap<TagCategory, Vec<usize>>,
    packed:  HashMap<TagCategory, PackedCategory>,
}

impl TagMaps {
    pub fn build(config: &Config) -> Self {
        let mut maps: HashMap<TagCategory, AHashMap<Vec<u8>, String>> = HashMap::new();
        let mut packed: HashMap<TagCategory, PackedCategory> = HashMap::new();
        // Use BTreeSet so the lengths come out sorted ascending — matches Java TreeSet.
        let mut length_sets: HashMap<TagCategory, BTreeSet<usize>> = HashMap::new();

        for tag_def in &config.tag_defs {
            if tag_def.category == TagCategory::Spacer {
                // SPACER is handled by the matcher directly; no lookup map needed.
                continue;
            }

            let category_map  = maps.entry(tag_def.category).or_default();
            let length_set    = length_sets.entry(tag_def.category).or_default();

            // Record the canonical length (Hamming expansion preserves length).
            length_set.insert(tag_def.seq.len());

            let neighbors = generate_neighbors(&tag_def.seq, tag_def.mismatches);
            for neighbor in neighbors {
                // Ambiguity check: two different canonical names map to the same seq.
                if let Some(existing) = category_map.get(&neighbor) {
                    if existing != &tag_def.name {
                        log::warn!(
                            "Ambiguity: {} is too close in sequence to {}",
                            tag_def.name, existing
                        );
                    }
                }
                category_map.insert(neighbor, tag_def.name.clone());
            }

            // Packed fast-path: canonical seq only (≤32 bp, no N).
            if let Some(p) = encode_packed(&tag_def.seq) {
                let pc = packed.entry(tag_def.category).or_insert_with(|| PackedCategory {
                    tags: Vec::new(),
                });
                pc.tags.push(PackedTag { packed: p, len: tag_def.seq.len(), name: tag_def.name.clone() });
            }
        }

        let lengths: HashMap<TagCategory, Vec<usize>> = length_sets
            .into_iter()
            .map(|(cat, set)| (cat, set.into_iter().collect()))
            .collect();

        Self { maps, lengths, packed }
    }

    /// Look up a sequence in a category's map. Returns the canonical tag name
    /// if found.
    pub fn lookup(&self, category: TagCategory, seq: &[u8]) -> Option<&str> {
        self.maps.get(&category)?.get(seq).map(String::as_str)
    }

    /// Hybrid lookup: tries a packed exact-match fast path first, then falls
    /// back to the AHashMap (which covers k=1/k=2 neighbors and N-containing
    /// sequences).
    ///
    /// Use this in the hot path (matcher) instead of `lookup`.
    pub fn lookup_packed(&self, category: TagCategory, seq: &[u8]) -> Option<&str> {
        // Fast path: distance-0 match via packed u64 comparison.
        if let (Some(pc), Some(query)) = (
            self.packed.get(&category),
            encode_packed(seq),
        ) {
            for pt in &pc.tags {
                if pt.packed == query && pt.len == seq.len() {
                    return Some(&pt.name);
                }
            }
            // No exact match found via packed — the sequence might be a Hamming
            // neighbor; fall through to AHashMap.
        }
        // Fallback: AHashMap covers neighbors (pre-expanded at build time) and
        // N-containing sequences.
        self.lookup(category, seq)
    }

    /// Sorted (ascending) list of unique sequence lengths for a category.
    /// The matcher iterates these lengths at each advance position.
    pub fn lengths(&self, category: TagCategory) -> &[usize] {
        self.lengths
            .get(&category)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{Config, TagCategory};

    // ---- Hamming expansion ----

    #[test]
    fn test_hamming_k0_returns_only_original() {
        let neighbors = generate_neighbors(b"ATCG", 0);
        assert_eq!(neighbors.len(), 1);
        assert!(neighbors.contains(&b"ATCG".to_vec()));
    }

    #[test]
    fn test_hamming_k1_contains_original() {
        let neighbors = generate_neighbors(b"ATCG", 1);
        assert!(neighbors.contains(&b"ATCG".to_vec()));
    }

    #[test]
    fn test_hamming_k1_contains_single_substitutions() {
        let neighbors = generate_neighbors(b"ATCG", 1);
        // pos 0: A → C
        assert!(neighbors.contains(&b"CTCG".to_vec()));
        // pos 1: T → A
        assert!(neighbors.contains(&b"AACG".to_vec()));
        // pos 2: C → N
        assert!(neighbors.contains(&b"ATNG".to_vec()));
        // pos 3: G → A
        assert!(neighbors.contains(&b"ATCA".to_vec()));
    }

    #[test]
    fn test_hamming_k1_does_not_contain_two_mismatches() {
        let neighbors = generate_neighbors(b"ATCG", 1);
        // Two substitutions: A→C and T→A
        assert!(!neighbors.contains(&b"CACG".to_vec()));
    }

    #[test]
    fn test_hamming_k2_contains_two_mismatch_sequence() {
        let neighbors = generate_neighbors(b"ATCG", 2);
        // Two mismatches
        assert!(neighbors.contains(&b"CACG".to_vec())); // pos0 A→C, pos1 T→A... wait
        // pos0 A→C gives CTCG; then pos1 T→A gives CACG — yes, 2 mismatches from original
        assert!(neighbors.contains(&b"CACG".to_vec()));
    }

    #[test]
    fn test_hamming_n_is_valid_substitution() {
        let neighbors = generate_neighbors(b"ATCG", 1);
        // N should be one of the substitution options at every position
        assert!(neighbors.contains(&b"NTCG".to_vec())); // pos 0
        assert!(neighbors.contains(&b"ANCG".to_vec())); // pos 1
        assert!(neighbors.contains(&b"ATNG".to_vec())); // pos 2
        assert!(neighbors.contains(&b"ATCN".to_vec())); // pos 3
    }

    // ---- encode_packed ----

    #[test]
    fn test_encode_packed_acgt() {
        // A=00 C=01 G=10 T=11 → ACGT = 0b_00_01_10_11 = 27
        assert_eq!(encode_packed(b"ACGT"), Some(0b_00_01_10_11));
    }

    #[test]
    fn test_encode_packed_all_a() {
        assert_eq!(encode_packed(b"AAAA"), Some(0));
    }

    #[test]
    fn test_encode_packed_all_t() {
        // TTTT = 11_11_11_11 = 0xFF
        assert_eq!(encode_packed(b"TTTT"), Some(0xFF));
    }

    #[test]
    fn test_encode_packed_empty() {
        assert_eq!(encode_packed(b""), Some(0));
    }

    #[test]
    fn test_encode_packed_n_returns_none() {
        assert_eq!(encode_packed(b"ACGN"), None);
        assert_eq!(encode_packed(b"NACGT"), None);
    }

    #[test]
    fn test_encode_packed_long_seq_returns_none() {
        let seq: Vec<u8> = b"ACGT".iter().cycle().take(33).cloned().collect();
        assert_eq!(encode_packed(&seq), None);
    }

    // ---- count_mismatches_packed ----

    #[test]
    fn test_count_mismatches_packed_identical() {
        assert_eq!(count_mismatches_packed(0, 8), 0);
    }

    #[test]
    fn test_count_mismatches_packed_one_mismatch() {
        // A=00, C=01 → XOR = 01 → 1 mismatch
        let a = encode_packed(b"A").unwrap();
        let c = encode_packed(b"C").unwrap();
        assert_eq!(count_mismatches_packed(a ^ c, 1), 1);
    }

    #[test]
    fn test_count_mismatches_packed_all_differ() {
        // ACGT vs TGCA — all 4 positions differ
        let seq1 = encode_packed(b"ACGT").unwrap();
        let seq2 = encode_packed(b"TGCA").unwrap();
        assert_eq!(count_mismatches_packed(seq1 ^ seq2, 4), 4);
    }

    #[test]
    fn test_count_mismatches_packed_two_mismatches() {
        // AAAA vs ACAT — positions 1 and 3 differ
        let seq1 = encode_packed(b"AAAA").unwrap();
        let seq2 = encode_packed(b"ACAT").unwrap();
        assert_eq!(count_mismatches_packed(seq1 ^ seq2, 4), 2);
    }

    // ---- TagMaps ----

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

    fn mini_maps() -> TagMaps {
        let cfg = Config::from_str_content(MINI_CONFIG).unwrap();
        TagMaps::build(&cfg)
    }

    #[test]
    fn test_exact_match_dpm() {
        let maps = mini_maps();
        assert_eq!(maps.lookup(TagCategory::Dpm, b"AAAACCCC"), Some("DPM_A"));
        assert_eq!(maps.lookup(TagCategory::Dpm, b"GGGGTTTT"), Some("DPM_B"));
    }

    #[test]
    fn test_no_match_wrong_category() {
        let maps = mini_maps();
        // DPM_A sequence should not match in the ODD map
        assert_eq!(maps.lookup(TagCategory::Odd, b"AAAACCCC"), None);
    }

    #[test]
    fn test_hamming1_neighbor_resolves_to_canonical() {
        let maps = mini_maps();
        // Odd_1 = ATCGATCG with k=1; ATCGATCC (last G→C) should resolve to Odd_1
        assert_eq!(maps.lookup(TagCategory::Odd, b"ATCGATCC"), Some("Odd_1"));
        // ATCGATCT (last G→T) should also resolve to Odd_1
        assert_eq!(maps.lookup(TagCategory::Odd, b"ATCGATCT"), Some("Odd_1"));
    }

    #[test]
    fn test_unknown_sequence_returns_none() {
        let maps = mini_maps();
        assert_eq!(maps.lookup(TagCategory::Dpm, b"XXXXXXXXX"), None);
    }

    #[test]
    fn test_lengths_dpm() {
        let maps = mini_maps();
        // Both DPM tags are 8 bp
        assert_eq!(maps.lengths(TagCategory::Dpm), &[8]);
    }

    #[test]
    fn test_lengths_y_variable() {
        // Build a config where Y tags have two different lengths
        let content = "READ2 = Y\nY\tY_short\tTTTT\t0\nY\tY_long\tTTTTAAAA\t0\n";
        let cfg = Config::from_str_content(content).unwrap();
        let maps = TagMaps::build(&cfg);
        // Lengths should be [4, 8] sorted ascending
        assert_eq!(maps.lengths(TagCategory::Y), &[4, 8]);
    }

    #[test]
    fn test_spacer_category_not_in_maps() {
        // SPACER entries in the config should not produce a lookup map
        let maps = mini_maps();
        assert_eq!(maps.lengths(TagCategory::Spacer), &[]);
        assert_eq!(maps.lookup(TagCategory::Spacer, b"XXXXXX"), None);
    }

    // ---- lookup_packed ----

    #[test]
    fn test_lookup_packed_exact_match() {
        let maps = mini_maps();
        assert_eq!(maps.lookup_packed(TagCategory::Dpm, b"AAAACCCC"), Some("DPM_A"));
        assert_eq!(maps.lookup_packed(TagCategory::Dpm, b"GGGGTTTT"), Some("DPM_B"));
        assert_eq!(maps.lookup_packed(TagCategory::Y,   b"TTTTAAAA"), Some("Y_1"));
    }

    #[test]
    fn test_lookup_packed_no_match() {
        let maps = mini_maps();
        assert_eq!(maps.lookup_packed(TagCategory::Dpm, b"XXXXXXXXX"), None);
    }

    #[test]
    fn test_lookup_packed_neighbor_falls_back_to_ahashmap() {
        let maps = mini_maps();
        // ATCGATCC is Hamming-1 from Odd_1=ATCGATCG → packed exact path misses it,
        // AHashMap fallback finds it.
        assert_eq!(maps.lookup_packed(TagCategory::Odd, b"ATCGATCC"), Some("Odd_1"));
        assert_eq!(maps.lookup_packed(TagCategory::Odd, b"ATCGATCT"), Some("Odd_1"));
    }

    #[test]
    fn test_lookup_packed_n_containing_falls_back() {
        let maps = mini_maps();
        // DPM has mismatches=0, so N-substituted DPM sequences are not in AHashMap.
        assert_eq!(maps.lookup_packed(TagCategory::Dpm, b"AAANCCCC"), None);
        // ODD has mismatches=1, so N-substituted ODD sequences ARE in AHashMap.
        assert_eq!(maps.lookup_packed(TagCategory::Odd, b"ATCGATCN"), Some("Odd_1"));
    }

    #[test]
    fn test_lookup_packed_no_length_collision() {
        // encode_packed strips leading A-bits, so sequences of different lengths
        // that share a common suffix encode identically:
        //   encode_packed(b"CC")   = Some(5)   (C=01: 1, C=01: (1<<2)|1 = 5)
        //   encode_packed(b"AACC") = Some(5)   (A=0,A=0,C=1,C=5 — same result)
        // A category with only a len-4 tag "AACC" must NOT match a len-2 query "CC".
        assert_eq!(encode_packed(b"CC"), encode_packed(b"AACC"),
            "Precondition: these two seqs must collide for the test to be meaningful");

        let content = "READ1 = EVEN\nEVEN\tEVEN_AACC\tAACC\t0\n";
        let cfg  = Config::from_str_content(content).unwrap();
        let maps = TagMaps::build(&cfg);

        // Exact match at the correct length still works.
        assert_eq!(maps.lookup_packed(TagCategory::Even, b"AACC"), Some("EVEN_AACC"));

        // A shorter query that collides in packed space must return None,
        // not the name of the longer tag.
        assert_eq!(maps.lookup_packed(TagCategory::Even, b"CC"), None);
    }

    #[test]
    fn test_lookup_packed_agrees_with_lookup() {
        let maps = mini_maps();
        let cases: &[(&[u8], TagCategory)] = &[
            (b"AAAACCCC", TagCategory::Dpm),
            (b"GGGGTTTT", TagCategory::Dpm),
            (b"TTTTAAAA", TagCategory::Y),
            (b"CCCCGGGG", TagCategory::Y),
            (b"ATCGATCG", TagCategory::Odd),
            (b"ATCGATCC", TagCategory::Odd),  // k=1 neighbor
            (b"ATCGATCN", TagCategory::Odd),  // N-containing k=1 neighbor
            (b"XXXXXXXX", TagCategory::Dpm),  // unknown
        ];
        for &(seq, cat) in cases {
            assert_eq!(
                maps.lookup_packed(cat, seq),
                maps.lookup(cat, seq),
                "Mismatch for {:?} in {:?}", seq, cat
            );
        }
    }
}
