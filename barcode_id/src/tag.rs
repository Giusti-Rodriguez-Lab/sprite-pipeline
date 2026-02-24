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
// Tag lookup maps
// ---------------------------------------------------------------------------

/// Pre-computed per-category lookup tables built from the config file.
///
/// For each tag category:
///   - `maps[category]` maps a sequence (bytes) to the canonical tag name.
///     Hamming-neighbor sequences also map to the canonical name.
///   - `lengths[category]` is a sorted list of unique tag lengths for that
///     category (used by the sliding-window matcher to try each length).
pub struct TagMaps {
    maps:    HashMap<TagCategory, AHashMap<Vec<u8>, String>>,
    lengths: HashMap<TagCategory, Vec<usize>>,
}

impl TagMaps {
    pub fn build(config: &Config) -> Self {
        let mut maps: HashMap<TagCategory, AHashMap<Vec<u8>, String>> = HashMap::new();
        // Use BTreeSet so the lengths come out sorted ascending — matches Java TreeSet.
        let mut length_sets: HashMap<TagCategory, BTreeSet<usize>> = HashMap::new();

        for tag_def in &config.tag_defs {
            if tag_def.category == TagCategory::Spacer {
                // SPACER is handled by the matcher directly; no lookup map needed.
                continue;
            }

            let category_map = maps.entry(tag_def.category).or_default();
            let length_set   = length_sets.entry(tag_def.category).or_default();

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
        }

        let lengths: HashMap<TagCategory, Vec<usize>> = length_sets
            .into_iter()
            .map(|(cat, set)| (cat, set.into_iter().collect()))
            .collect();

        Self { maps, lengths }
    }

    /// Look up a sequence in a category's map. Returns the canonical tag name
    /// if found.
    pub fn lookup(&self, category: TagCategory, seq: &[u8]) -> Option<&str> {
        self.maps.get(&category)?.get(seq).map(String::as_str)
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
}
