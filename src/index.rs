use crate::kmer::{KmerSimilarityTable, N_DIMERS};
use crate::model::Node;

pub const N_COMPOSITE_KEYS: usize = N_DIMERS * N_DIMERS;

/// Largest number of composite keys a single node can occupy.
pub const MAX_COMPOSITE_KEYS: usize = 9;

/// Which residue-column pairs inside a terminal 3-mer are indexed.
///
/// The constrained full alignment requires at least
/// `--minimum-terminal-match-length` residue-to-residue columns drawn from the
/// first three residues of both peptides, and the same at the C terminus. For
/// the default of two, those columns may sit at any ordered position pair
/// `(i1 < i2)` of one peptide against any `(j1 < j2)` of the other. Indexing
/// only the contiguous pairs therefore cannot retrieve a pair whose required
/// terminal columns are shifted, which is the common case when the two
/// peptides differ in length.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TerminalSeed {
    /// Historical geometry: contiguous dimers (1,2) and (2,3) only.
    Contiguous,
    /// All three ordered column pairs (1,2), (1,3) and (2,3).
    AllColumnPairs,
}

impl TerminalSeed {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "contiguous" => Ok(Self::Contiguous),
            "all_column_pairs" | "all-column-pairs" => Ok(Self::AllColumnPairs),
            _ => Err(format!("invalid --terminal-seed: {value}")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Contiguous => "contiguous",
            Self::AllColumnPairs => "all_column_pairs",
        }
    }

    /// Ordered position pairs inside one terminal 3-mer.
    #[inline]
    fn column_pairs(self) -> &'static [(usize, usize)] {
        match self {
            Self::Contiguous => &[(0, 1), (1, 2)],
            Self::AllColumnPairs => &[(0, 1), (0, 2), (1, 2)],
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct KeyRelation {
    pub first: u32,
    pub second: u32,
}

#[derive(Clone, Copy, Debug)]
pub struct IndexStats {
    pub occupied_keys: usize,
    pub largest_bucket: usize,
    pub similar_key_relations: usize,
    /// Upper bound before node-pair deduplication. Cross-bucket intersections
    /// can include a node paired to itself; actual generation removes these.
    pub candidate_occurrence_upper_bound: u128,
}

/// Ordered 2-mer codes for one terminus. `offset` is 0 for the N-terminal
/// 3-mer and 3 for the C-terminal 3-mer of `node.codes`.
#[inline]
pub fn terminal_dimers(node: &Node, offset: usize, seed: TerminalSeed) -> ([u32; 3], usize) {
    let pairs = seed.column_pairs();
    let mut result = [0u32; 3];
    for (slot, &(first, second)) in pairs.iter().enumerate() {
        result[slot] = node.codes[offset + first] as u32 * 20 + node.codes[offset + second] as u32;
    }
    (result, pairs.len())
}

/// Composite `front * N_DIMERS + end` keys occupied by one node. Retrieval
/// expands the neighbours of both components, so enumerating the full cross
/// product makes candidacy equivalent to "at least one neighbouring front
/// column pair AND at least one neighbouring end column pair".
#[inline]
pub fn composite_keys(node: &Node, seed: TerminalSeed) -> ([u32; MAX_COMPOSITE_KEYS], usize) {
    let (front, n_front) = terminal_dimers(node, 0, seed);
    let (back, n_back) = terminal_dimers(node, 3, seed);
    let mut keys = [0u32; MAX_COMPOSITE_KEYS];
    let mut count = 0usize;
    for f in 0..n_front {
        for b in 0..n_back {
            keys[count] = front[f] * N_DIMERS as u32 + back[b];
            count += 1;
        }
    }
    keys[..count].sort_unstable();
    let mut unique = 0usize;
    for i in 0..count {
        if i == 0 || keys[i] != keys[unique - 1] {
            keys[unique] = keys[i];
            unique += 1;
        }
    }
    (keys, unique)
}

pub fn build_exact_index(nodes: &[Node], seed: TerminalSeed) -> Vec<Vec<u32>> {
    build_exact_index_subset(nodes, 0..nodes.len(), seed)
}

pub fn build_exact_index_subset<I>(nodes: &[Node], node_ids: I, seed: TerminalSeed) -> Vec<Vec<u32>>
where
    I: IntoIterator<Item = usize>,
{
    let mut buckets = vec![Vec::<u32>::new(); N_COMPOSITE_KEYS];
    for node_id in node_ids {
        let (keys, count) = composite_keys(&nodes[node_id], seed);
        for &key in &keys[..count] {
            buckets[key as usize].push(node_id as u32);
        }
    }
    buckets
}

/// Retrieve distinct nodes sharing at least one neighbouring front column pair
/// and at least one neighbouring end column pair with `query`. The returned
/// identifiers are canonical and independent of posting-list traversal order.
pub fn retrieve_candidates(
    query_id: Option<u32>,
    query: &Node,
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
) -> Vec<u32> {
    let mut result = Vec::<u32>::new();
    let (keys, count) = composite_keys(query, seed);
    for &key in &keys[..count] {
        let front = key as usize / N_DIMERS;
        let back = key as usize % N_DIMERS;
        for &other_front in table.neighbours(front) {
            for &other_back in table.neighbours(back) {
                let other = other_front as usize * N_DIMERS + other_back as usize;
                result.extend_from_slice(&buckets[other]);
            }
        }
    }
    result.sort_unstable();
    result.dedup();
    if let Some(id) = query_id {
        if let Ok(position) = result.binary_search(&id) {
            result.remove(position);
        }
    }
    result
}

pub fn build_similar_key_relations(
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
) -> (Vec<KeyRelation>, IndexStats) {
    let mut relations = Vec::<KeyRelation>::new();
    let mut occupied_keys = 0usize;
    let mut largest_bucket = 0usize;
    for bucket in buckets {
        if !bucket.is_empty() {
            occupied_keys += 1;
            largest_bucket = largest_bucket.max(bucket.len());
        }
    }

    let mut candidate_occurrence_upper_bound = 0u128;
    for first in 0..N_COMPOSITE_KEYS {
        if buckets[first].is_empty() {
            continue;
        }
        let front = first / N_DIMERS;
        let back = first % N_DIMERS;
        for &other_front in table.neighbours(front) {
            for &other_back in table.neighbours(back) {
                let second = other_front as usize * N_DIMERS + other_back as usize;
                if second < first || buckets[second].is_empty() {
                    continue;
                }
                relations.push(KeyRelation {
                    first: first as u32,
                    second: second as u32,
                });
                let n_first = buckets[first].len() as u128;
                let n_second = buckets[second].len() as u128;
                if first == second {
                    candidate_occurrence_upper_bound += n_first * n_first.saturating_sub(1) / 2;
                } else {
                    candidate_occurrence_upper_bound += n_first * n_second;
                }
            }
        }
    }
    relations.sort_unstable_by_key(|relation| (relation.first, relation.second));
    let similar_key_relations = relations.len();
    (
        relations,
        IndexStats {
            occupied_keys,
            largest_bucket,
            similar_key_relations,
            candidate_occurrence_upper_bound,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::aa_code;

    fn node(anchor: &[u8; 6]) -> Node {
        Node {
            sequence: anchor.to_vec(),
            sequence_codes: anchor.iter().copied().map(aa_code).collect(),
            anchor: *anchor,
            codes: anchor.map(aa_code),
            combination_mask: 0b11_1111,
            frequency: 1,
        }
    }

    fn shares_key(a: &Node, b: &Node, seed: TerminalSeed) -> bool {
        let (ka, na) = composite_keys(a, seed);
        let (kb, nb) = composite_keys(b, seed);
        ka[..na].iter().any(|key| kb[..nb].contains(key))
    }

    #[test]
    fn shifted_example_shares_an_exact_composite_key() {
        let a = node(b"GALKLV");
        let b = node(b"ALVLVI");
        assert!(shares_key(&a, &b, TerminalSeed::Contiguous));
        assert!(shares_key(&a, &b, TerminalSeed::AllColumnPairs));
    }

    #[test]
    fn composite_key_counts_match_the_geometry() {
        let a = node(b"GALKLV");
        assert_eq!(composite_keys(&a, TerminalSeed::Contiguous).1, 4);
        assert_eq!(composite_keys(&a, TerminalSeed::AllColumnPairs).1, 9);
        // A homopolymeric terminus collapses to a single distinct key.
        let flat = node(b"AAAAAA");
        assert_eq!(composite_keys(&flat, TerminalSeed::Contiguous).1, 1);
        assert_eq!(composite_keys(&flat, TerminalSeed::AllColumnPairs).1, 1);
    }

    #[test]
    fn keys_are_sorted_and_deduplicated() {
        for seed in [TerminalSeed::Contiguous, TerminalSeed::AllColumnPairs] {
            let (keys, count) = composite_keys(&node(b"GALKLV"), seed);
            assert!(keys[..count].windows(2).all(|w| w[0] < w[1]));
        }
    }

    /// Both terminal columns of every contiguous dimer touch the middle
    /// residue, so a single middle-position substitution destroys both of them
    /// while the spaced (1,3) column pair survives it untouched. This is the
    /// mechanism behind the measured 36% candidate-recall loss.
    #[test]
    fn a_middle_terminal_mismatch_only_matches_through_the_spaced_column_pair() {
        let a = node(b"AWDKLV");
        let b = node(b"APDKLV");
        assert!(!shares_key(&a, &b, TerminalSeed::Contiguous));
        assert!(shares_key(&a, &b, TerminalSeed::AllColumnPairs));
    }
}
