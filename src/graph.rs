use crate::edge_store::read_edge;
use crate::fasta::DynError;
use crate::model::{Graph, Neighbor, Node};
use crate::scoring::Scorer;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::BufReader;
use std::path::Path;

pub fn load_graph(
    edge_file: &Path,
    node_count: usize,
    expected_edges: u64,
    max_memory_bytes: u64,
) -> Result<Graph, DynError> {
    let mut degrees = vec![0usize; node_count];
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, File::open(edge_file)?);
    let mut edge_count = 0u64;
    while let Some(edge) = read_edge(&mut reader)? {
        if edge.u as usize >= node_count || edge.v as usize >= node_count || edge.u >= edge.v {
            return Err("invalid edge record while constructing graph".into());
        }
        degrees[edge.u as usize] += 1;
        degrees[edge.v as usize] += 1;
        edge_count += 1;
    }
    if edge_count != expected_edges {
        return Err(format!(
            "edge-file count mismatch: expected {expected_edges}, observed {edge_count}"
        )
        .into());
    }
    let estimated = (node_count as u64 + 1) * std::mem::size_of::<usize>() as u64
        + edge_count * 2 * std::mem::size_of::<Neighbor>() as u64
        + node_count as u64 * std::mem::size_of::<usize>() as u64;
    if estimated > max_memory_bytes {
        return Err(format!(
            "graph needs approximately {:.2} GiB, above --max-memory-gb; increase the limit or use --index-only",
            estimated as f64 / 1024.0_f64.powi(3)
        )
        .into());
    }

    let mut offsets = vec![0usize; node_count + 1];
    for i in 0..node_count {
        offsets[i + 1] = offsets[i] + degrees[i];
    }
    let mut cursor = offsets[..node_count].to_vec();
    let mut neighbors = vec![Neighbor { node: 0, weight: 0 }; offsets[node_count]];
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, File::open(edge_file)?);
    while let Some(edge) = read_edge(&mut reader)? {
        let ui = edge.u as usize;
        let vi = edge.v as usize;
        neighbors[cursor[ui]] = Neighbor {
            node: edge.v,
            weight: edge.weight,
        };
        cursor[ui] += 1;
        neighbors[cursor[vi]] = Neighbor {
            node: edge.u,
            weight: edge.weight,
        };
        cursor[vi] += 1;
    }
    // The sorted edge stream normally creates sorted adjacency lists already.
    // Sorting explicitly protects binary lookup if the file format changes.
    for i in 0..node_count {
        neighbors[offsets[i]..offsets[i + 1]].sort_unstable_by_key(|n| n.node);
    }
    Ok(Graph {
        offsets,
        neighbors,
        edge_count,
    })
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct HeapEntry {
    coverage: u32,
    weight_sum: u64,
    frequency: u64,
    node: u32,
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.coverage
            .cmp(&other.coverage)
            .then(self.weight_sum.cmp(&other.weight_sum))
            .then(self.frequency.cmp(&other.frequency))
            // Smaller, alphabetically earlier node ids win the final tie.
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug)]
pub struct Clustering {
    pub cluster_of: Vec<u32>,
    pub representatives: Vec<u32>,
}

/// How the order of representative selection is decided.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RepresentativeOrder {
    /// Dynamic greedy set cover: repeatedly take the peptide covering the most
    /// still-unassigned peptides. Minimises the cluster count, but the key is a
    /// property of the whole dataset, so subsampling changes the selection.
    Coverage,
    /// Visit peptides in an order that depends only on the peptide itself. The
    /// order of any subset is the restriction of the full-dataset order, so
    /// selection no longer churns when the dataset composition changes; the
    /// residual difference comes only from representatives absent in the subset.
    Intrinsic,
}

impl RepresentativeOrder {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "coverage" => Ok(Self::Coverage),
            "intrinsic" => Ok(Self::Intrinsic),
            _ => Err(format!("invalid --representative-order: {value}")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Coverage => "coverage",
            Self::Intrinsic => "intrinsic",
        }
    }
}

/// Deterministic visiting order built only from per-peptide properties: longer
/// peptides first, then canonical sequence order. Input frequency is
/// deliberately excluded because it depends on which duplicates a sample
/// retained, which would reintroduce composition dependence.
pub fn intrinsic_order(nodes: &[Node]) -> Vec<u32> {
    let mut order: Vec<u32> = (0..nodes.len() as u32).collect();
    order.sort_unstable_by(|&a, &b| {
        let left = &nodes[a as usize];
        let right = &nodes[b as usize];
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    order
}

pub fn greedy_set_cover(nodes: &[Node], graph: &Graph, order: RepresentativeOrder) -> Clustering {
    greedy_set_cover_subset(nodes, graph, &vec![true; nodes.len()], order)
}

pub fn greedy_set_cover_subset(
    nodes: &[Node],
    graph: &Graph,
    eligible: &[bool],
    order: RepresentativeOrder,
) -> Clustering {
    if order == RepresentativeOrder::Intrinsic {
        return intrinsic_set_cover_subset(nodes, graph, eligible);
    }
    let n = nodes.len();
    assert_eq!(eligible.len(), n);
    let mut assigned: Vec<bool> = eligible.iter().map(|x| !*x).collect();
    let mut cluster_of = vec![u32::MAX; n];
    let mut coverage: Vec<u32> = (0..n)
        .map(|i| {
            if eligible[i] {
                1 + graph
                    .neighbors(i as u32)
                    .iter()
                    .filter(|x| eligible[x.node as usize])
                    .count() as u32
            } else {
                0
            }
        })
        .collect();
    let mut weight_sum: Vec<u64> = (0..n)
        .map(|i| {
            if eligible[i] {
                graph
                    .neighbors(i as u32)
                    .iter()
                    .filter(|x| eligible[x.node as usize])
                    .map(|x| x.weight as u64)
                    .sum()
            } else {
                0
            }
        })
        .collect();
    let mut heap = BinaryHeap::with_capacity(n);
    for node in 0..n {
        if eligible[node] {
            heap.push(HeapEntry {
                coverage: coverage[node],
                weight_sum: weight_sum[node],
                frequency: nodes[node].frequency,
                node: node as u32,
            });
        }
    }

    let mut representatives = Vec::new();
    while let Some(entry) = heap.pop() {
        let rep = entry.node as usize;
        if assigned[rep] || entry.coverage != coverage[rep] || entry.weight_sum != weight_sum[rep] {
            continue;
        }
        let cluster_id = representatives.len() as u32;
        representatives.push(rep as u32);
        let mut newly_assigned = Vec::with_capacity(coverage[rep] as usize);
        newly_assigned.push(rep as u32);
        for neighbor in graph.neighbors(rep as u32) {
            if !assigned[neighbor.node as usize] {
                newly_assigned.push(neighbor.node);
            }
        }
        newly_assigned.sort_unstable();
        newly_assigned.dedup();
        for &node in &newly_assigned {
            assigned[node as usize] = true;
            cluster_of[node as usize] = cluster_id;
        }
        for &removed in &newly_assigned {
            for neighbor in graph.neighbors(removed) {
                let candidate = neighbor.node as usize;
                if !assigned[candidate] {
                    coverage[candidate] -= 1;
                    weight_sum[candidate] -= neighbor.weight as u64;
                    heap.push(HeapEntry {
                        coverage: coverage[candidate],
                        weight_sum: weight_sum[candidate],
                        frequency: nodes[candidate].frequency,
                        node: candidate as u32,
                    });
                }
            }
        }
    }
    debug_assert!(cluster_of
        .iter()
        .enumerate()
        .all(|(i, cluster)| !eligible[i] || *cluster != u32::MAX));
    Clustering {
        cluster_of,
        representatives,
    }
}

/// Fixed-order greedy cover. Each peptide is visited once in `intrinsic_order`;
/// an unassigned peptide becomes a representative and absorbs every still
/// unassigned neighbour.
fn intrinsic_set_cover_subset(nodes: &[Node], graph: &Graph, eligible: &[bool]) -> Clustering {
    let n = nodes.len();
    assert_eq!(eligible.len(), n);
    let mut assigned: Vec<bool> = eligible.iter().map(|x| !*x).collect();
    let mut cluster_of = vec![u32::MAX; n];
    let mut representatives = Vec::new();
    for candidate in intrinsic_order(nodes) {
        if assigned[candidate as usize] {
            continue;
        }
        let cluster_id = representatives.len() as u32;
        representatives.push(candidate);
        assigned[candidate as usize] = true;
        cluster_of[candidate as usize] = cluster_id;
        for neighbor in graph.neighbors(candidate) {
            if !assigned[neighbor.node as usize] {
                assigned[neighbor.node as usize] = true;
                cluster_of[neighbor.node as usize] = cluster_id;
            }
        }
    }
    debug_assert!(cluster_of
        .iter()
        .enumerate()
        .all(|(i, cluster)| !eligible[i] || *cluster != u32::MAX));
    Clustering {
        cluster_of,
        representatives,
    }
}

/// Preserve non-singleton prefilter clusters and complete the remaining nodes
/// by set cover on the sensitive graph.
pub fn prefilter_masks(node_count: usize, prefilter: &Clustering) -> (Vec<bool>, Vec<bool>) {
    let members = cluster_members(&prefilter.cluster_of, prefilter.representatives.len());
    let mut unassigned = vec![true; node_count];
    let mut prefilter_representatives = vec![false; node_count];
    for (cluster, group) in members.iter().enumerate() {
        if group.len() > 1 {
            prefilter_representatives[prefilter.representatives[cluster] as usize] = true;
            for &node in group {
                unassigned[node as usize] = false;
            }
        }
    }
    (prefilter_representatives, unassigned)
}

fn cluster_members(cluster_of: &[u32], cluster_count: usize) -> Vec<Vec<u32>> {
    let mut members = vec![Vec::new(); cluster_count];
    for (node, &cluster) in cluster_of.iter().enumerate() {
        members[cluster as usize].push(node as u32);
    }
    members
}

fn update_representatives(
    nodes: &[Node],
    scorer: &Scorer,
    graph: &Graph,
    cluster_of: &[u32],
    current: &[u32],
) -> Vec<u32> {
    let members = cluster_members(cluster_of, current.len());
    members
        .par_iter()
        .enumerate()
        .map(|(cluster_id, cluster)| {
            let mut best = current[cluster_id];
            let mut best_key = (0u128, 0u128, 0u128, nodes[best as usize].frequency);
            for &candidate in cluster {
                let mut internal_neighbors = 0usize;
                let mut score = 0u128;
                let mut alignment = 0u128;
                let mut anchor = 0u128;
                for neighbor in graph.neighbors(candidate) {
                    if cluster_of[neighbor.node as usize] == cluster_id as u32 {
                        internal_neighbors += 1;
                        let frequency = nodes[neighbor.node as usize].frequency as u128;
                        score += neighbor.weight as u128 * frequency;
                        let components = scorer
                            .scores(&nodes[candidate as usize], &nodes[neighbor.node as usize]);
                        alignment += components.alignment as u128 * frequency;
                        anchor += components.anchor_combination as u128 * frequency;
                    }
                }
                if internal_neighbors + 1 != cluster.len() {
                    continue;
                }
                let frequency = nodes[candidate as usize].frequency;
                let key = (score, alignment, anchor, frequency);
                if key > best_key || (key == best_key && candidate < best) {
                    best = candidate;
                    best_key = key;
                }
            }
            best
        })
        .collect()
}

/// Synchronous reassignment.
///
/// `minimum_improvement` is hysteresis on the move decision: a peptide leaves
/// its current representative only when another beats it by more than this
/// margin, in the same thousandths as the ranking weight. Zero reproduces the
/// pre-0.5.0 behaviour, where any improvement at all — including an exact tie
/// broken by identifier — caused a move. Near-ties are the unstable case: which
/// representative wins a tie depends on which peptides happen to be in the
/// dataset, so with zero margin a small change in composition reshuffles a large
/// number of assignments.
pub fn refine(
    nodes: &[Node],
    scorer: &Scorer,
    graph: &Graph,
    mut clustering: Clustering,
    iterations: usize,
    minimum_improvement: u16,
) -> (Clustering, u64, usize) {
    let mut total_moves = 0u64;
    let mut completed = 0usize;
    for _ in 0..iterations {
        completed += 1;
        let mut representative_cluster = vec![u32::MAX; nodes.len()];
        for (cluster, &rep) in clustering.representatives.iter().enumerate() {
            representative_cluster[rep as usize] = cluster as u32;
        }
        let proposals: Vec<u32> = (0..nodes.len())
            .into_par_iter()
            .map(|node| {
                let current_cluster = clustering.cluster_of[node];
                if representative_cluster[node] != u32::MAX {
                    return current_cluster;
                }
                let current_rep = clustering.representatives[current_cluster as usize];
                let current_weight = graph.edge_weight(node as u32, current_rep).unwrap_or(0);
                let mut best_cluster = current_cluster;
                let mut best_rep = current_rep;
                let mut best_weight = current_weight;
                let current_scores = scorer
                    .eligible_pair_scores(&nodes[node], &nodes[current_rep as usize])
                    .expect("current graph representative must be eligible");
                let mut best_scores = current_scores;
                for neighbor in graph.neighbors(node as u32) {
                    let candidate_cluster = representative_cluster[neighbor.node as usize];
                    if candidate_cluster == u32::MAX {
                        continue;
                    }
                    let better = if neighbor.weight > best_weight {
                        true
                    } else if neighbor.weight < best_weight {
                        false
                    } else {
                        let scores = scorer
                            .eligible_pair_scores(&nodes[node], &nodes[neighbor.node as usize])
                            .expect("stored graph edge must be eligible");
                        let ordering = scorer.compare_scores(scores, best_scores);
                        if ordering.is_gt() || (ordering.is_eq() && neighbor.node < best_rep) {
                            best_scores = scores;
                            true
                        } else {
                            false
                        }
                    };
                    if better {
                        best_weight = neighbor.weight;
                        best_rep = neighbor.node;
                        best_cluster = candidate_cluster;
                        best_scores = scorer
                            .eligible_pair_scores(&nodes[node], &nodes[neighbor.node as usize])
                            .expect("stored graph edge must be eligible");
                    }
                }
                let component_improvement = minimum_improvement == 0
                    && best_weight == current_weight
                    && (scorer.compare_scores(best_scores, current_scores).is_gt()
                        || (scorer.compare_scores(best_scores, current_scores).is_eq()
                            && best_rep < current_rep));
                if best_cluster != current_cluster
                    && (best_weight as u32 > current_weight as u32 + minimum_improvement as u32
                        || component_improvement)
                {
                    best_cluster
                } else {
                    current_cluster
                }
            })
            .collect();
        let moves = proposals
            .iter()
            .zip(&clustering.cluster_of)
            .filter(|(a, b)| a != b)
            .count() as u64;
        clustering.cluster_of = proposals;
        let new_representatives = update_representatives(
            nodes,
            scorer,
            graph,
            &clustering.cluster_of,
            &clustering.representatives,
        );
        let rep_changes = new_representatives
            .iter()
            .zip(&clustering.representatives)
            .filter(|(a, b)| a != b)
            .count();
        clustering.representatives = new_representatives;
        total_moves += moves;
        if moves == 0 && rep_changes == 0 {
            break;
        }
    }
    (clustering, total_moves, completed)
}

pub fn strict_merge(
    graph: &Graph,
    mut clustering: Clustering,
    merge_cap: Option<usize>,
) -> (Clustering, usize) {
    let cluster_count = clustering.representatives.len();
    let mut members = cluster_members(&clustering.cluster_of, cluster_count);
    let mut active = vec![true; cluster_count];
    let mut order: Vec<usize> = (0..cluster_count).collect();
    order.sort_unstable_by(|&a, &b| {
        members[b]
            .len()
            .cmp(&members[a].len())
            .then(clustering.representatives[a].cmp(&clustering.representatives[b]))
    });
    let mut merge_count = 0usize;
    for &target in &order {
        if !active[target] {
            continue;
        }
        let rep = clustering.representatives[target];
        let mut counts: HashMap<u32, usize> = HashMap::new();
        for neighbor in graph.neighbors(rep) {
            let source = clustering.cluster_of[neighbor.node as usize];
            if source as usize != target && active[source as usize] {
                *counts.entry(source).or_insert(0) += 1;
            }
        }
        let mut candidates: Vec<usize> = counts
            .into_iter()
            .filter_map(|(source, count)| {
                let source = source as usize;
                (active[source]
                    && count == members[source].len()
                    && members[target].len() >= members[source].len())
                .then_some(source)
            })
            .collect();
        candidates.sort_unstable_by(|&a, &b| {
            members[b]
                .len()
                .cmp(&members[a].len())
                .then(clustering.representatives[a].cmp(&clustering.representatives[b]))
        });
        for source in candidates {
            if !active[source] || members[target].len() < members[source].len() {
                continue;
            }
            if let Some(cap) = merge_cap {
                if members[source]
                    .iter()
                    .take(cap)
                    .any(|node| graph.edge_weight(rep, *node).is_none())
                {
                    continue;
                }
            }
            // The cap can reject early, but successful merges are always
            // validated against every member.
            if members[source]
                .iter()
                .any(|node| graph.edge_weight(rep, *node).is_none())
            {
                continue;
            }
            let source_members = std::mem::take(&mut members[source]);
            for node in &source_members {
                clustering.cluster_of[*node as usize] = target as u32;
            }
            members[target].extend(source_members);
            active[source] = false;
            merge_count += 1;
        }
    }

    let mut surviving: Vec<usize> = (0..cluster_count).filter(|c| active[*c]).collect();
    surviving.sort_unstable_by_key(|c| clustering.representatives[*c]);
    let mut remap = vec![u32::MAX; cluster_count];
    let mut representatives = Vec::with_capacity(surviving.len());
    for (new, old) in surviving.into_iter().enumerate() {
        remap[old] = new as u32;
        representatives.push(clustering.representatives[old]);
    }
    for cluster in &mut clustering.cluster_of {
        *cluster = remap[*cluster as usize];
    }
    clustering.representatives = representatives;
    (clustering, merge_count)
}

#[derive(Clone, Copy, Debug, Default)]
pub struct IterationStats {
    pub iterations: usize,
    pub reassignment_moves: u64,
    pub representative_changes: u64,
    pub merges: usize,
    pub validation_failures: u64,
    pub converged: bool,
}

fn state_hash(clustering: &Clustering) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    clustering.cluster_of.hash(&mut hasher);
    clustering.representatives.hash(&mut hasher);
    hasher.finish()
}

#[allow(clippy::too_many_arguments)]
pub fn iterate_to_convergence(
    nodes: &[Node],
    scorer: &Scorer,
    graph: &Graph,
    mut clustering: Clustering,
    iteration_cap: Option<usize>,
    merge: bool,
    merge_cap: Option<usize>,
    minimum_improvement: u16,
) -> (Clustering, IterationStats) {
    let mut stats = IterationStats::default();
    let mut seen = HashSet::new();
    loop {
        if iteration_cap.is_some_and(|cap| stats.iterations >= cap) {
            break;
        }
        if !seen.insert(state_hash(&clustering)) {
            break;
        }
        stats.iterations += 1;
        let before_reps = clustering.representatives.clone();
        let (refined, moves, _) = refine(
            nodes,
            scorer,
            graph,
            clustering,
            1,
            minimum_improvement,
        );
        clustering = refined;
        let rep_changes = before_reps
            .iter()
            .zip(&clustering.representatives)
            .filter(|(a, b)| a != b)
            .count() as u64;
        let mut merges = 0usize;
        if merge {
            let result = strict_merge(graph, clustering, merge_cap);
            clustering = result.0;
            merges = result.1;
            if merges > 0 {
                let current = clustering.representatives.clone();
                clustering.representatives =
                    update_representatives(nodes, scorer, graph, &clustering.cluster_of, &current);
                stats.representative_changes += current
                    .iter()
                    .zip(&clustering.representatives)
                    .filter(|(a, b)| a != b)
                    .count() as u64;
            }
        }
        let (covered, total) = representative_coverage(graph, &clustering);
        let failures = total - covered;
        stats.validation_failures += failures;
        stats.reassignment_moves += moves;
        stats.representative_changes += rep_changes;
        stats.merges += merges;
        if moves == 0 && rep_changes == 0 && merges == 0 && failures == 0 {
            stats.converged = true;
            break;
        }
        // A failure should be impossible: every assignment, representative
        // update and accepted merge requires an eligible graph edge.
        if failures > 0 {
            break;
        }
    }
    (clustering, stats)
}

pub fn representative_coverage(graph: &Graph, clustering: &Clustering) -> (u64, u64) {
    let mut covered = 0u64;
    for (node, &cluster) in clustering.cluster_of.iter().enumerate() {
        let rep = clustering.representatives[cluster as usize];
        if graph.edge_weight(node as u32, rep).is_some() {
            covered += 1;
        }
    }
    (covered, clustering.cluster_of.len() as u64)
}

pub fn canonicalize(mut clustering: Clustering) -> Clustering {
    let mut order: Vec<usize> = (0..clustering.representatives.len()).collect();
    order.sort_unstable_by_key(|cluster| clustering.representatives[*cluster]);
    let mut remap = vec![0u32; order.len()];
    let mut representatives = Vec::with_capacity(order.len());
    for (new, old) in order.into_iter().enumerate() {
        remap[old] = new as u32;
        representatives.push(clustering.representatives[old]);
    }
    for cluster in &mut clustering.cluster_of {
        *cluster = remap[*cluster as usize];
    }
    clustering.representatives = representatives;
    clustering
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
    fn set_cover_makes_directly_covered_stars() {
        let nodes = vec![
            node(b"AAAAAA"),
            node(b"AAAAAR"),
            node(b"AAAAAN"),
            node(b"CCCCCC"),
        ];
        let graph = Graph {
            offsets: vec![0, 2, 3, 4, 4],
            neighbors: vec![
                Neighbor {
                    node: 1,
                    weight: 900,
                },
                Neighbor {
                    node: 2,
                    weight: 800,
                },
                Neighbor {
                    node: 0,
                    weight: 900,
                },
                Neighbor {
                    node: 0,
                    weight: 800,
                },
            ],
            edge_count: 2,
        };
        let clustering = greedy_set_cover(&nodes, &graph, RepresentativeOrder::Coverage);
        assert_eq!(clustering.representatives.len(), 2);
        assert_eq!(representative_coverage(&graph, &clustering), (4, 4));
        let fixed = greedy_set_cover(&nodes, &graph, RepresentativeOrder::Intrinsic);
        assert_eq!(fixed.representatives.len(), 2);
        assert_eq!(representative_coverage(&graph, &fixed), (4, 4));
    }

    /// The intrinsic order is a function of the peptide alone, so removing
    /// peptides never reorders the survivors. This is the property that makes
    /// subset selection stop churning.
    #[test]
    fn intrinsic_order_is_stable_under_subsetting() {
        let all = vec![
            node(b"AAAAAA"),
            node(b"CCCCCC"),
            node(b"AAAAAR"),
            node(b"AAAAAN"),
        ];
        let full: Vec<Vec<u8>> = intrinsic_order(&all)
            .into_iter()
            .map(|id| all[id as usize].sequence.clone())
            .collect();
        let kept: Vec<Node> = vec![all[1].clone(), all[3].clone()];
        let subset: Vec<Vec<u8>> = intrinsic_order(&kept)
            .into_iter()
            .map(|id| kept[id as usize].sequence.clone())
            .collect();
        let restricted: Vec<Vec<u8>> = full
            .into_iter()
            .filter(|sequence| subset.contains(sequence))
            .collect();
        assert_eq!(restricted, subset);
    }
}
