use crate::kmer::{KmerSimilarityTable, N_DIMERS};
use crate::model::Node;

pub const N_COMPOSITE_KEYS: usize = N_DIMERS * N_DIMERS;

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

#[inline]
pub fn composite_keys(node: &Node) -> [u32; 4] {
    let front0 = node.codes[0] as u32 * 20 + node.codes[1] as u32;
    let front1 = node.codes[1] as u32 * 20 + node.codes[2] as u32;
    let back0 = node.codes[3] as u32 * 20 + node.codes[4] as u32;
    let back1 = node.codes[4] as u32 * 20 + node.codes[5] as u32;
    [
        front0 * N_DIMERS as u32 + back0,
        front0 * N_DIMERS as u32 + back1,
        front1 * N_DIMERS as u32 + back0,
        front1 * N_DIMERS as u32 + back1,
    ]
}

pub fn build_exact_index(nodes: &[Node]) -> Vec<Vec<u32>> {
    build_exact_index_subset(nodes, 0..nodes.len())
}

pub fn build_exact_index_subset<I>(nodes: &[Node], node_ids: I) -> Vec<Vec<u32>>
where
    I: IntoIterator<Item = usize>,
{
    let mut buckets = vec![Vec::<u32>::new(); N_COMPOSITE_KEYS];
    for node_id in node_ids {
        let node = &nodes[node_id];
        let mut keys = composite_keys(node);
        keys.sort_unstable();
        for i in 0..4 {
            if i == 0 || keys[i] != keys[i - 1] {
                buckets[keys[i] as usize].push(node_id as u32);
            }
        }
    }
    buckets
}

/// Retrieve distinct nodes sharing at least one similar front dimer and one
/// similar end dimer with `query`. The returned identifiers are canonical and
/// independent of posting-list traversal order.
pub fn retrieve_candidates(
    query_id: Option<u32>,
    query: &Node,
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
) -> Vec<u32> {
    let mut result = Vec::<u32>::new();
    let mut keys = composite_keys(query);
    keys.sort_unstable();
    for (position, &key) in keys.iter().enumerate() {
        if position > 0 && key == keys[position - 1] {
            continue;
        }
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

    #[test]
    fn shifted_example_shares_an_exact_composite_key() {
        let a = composite_keys(&node(b"GALKLV"));
        let b = composite_keys(&node(b"ALVLVI"));
        assert!(a.iter().any(|key| b.contains(key)));
    }
}
