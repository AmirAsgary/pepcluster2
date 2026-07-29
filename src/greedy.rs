use crate::graph::{canonicalize, Clustering, IterationStats};
use crate::index::{build_exact_index_subset, retrieve_candidates};
use crate::kmer::KmerSimilarityTable;
use crate::model::Node;
use crate::scoring::{Scorer, SimilarityScores};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};
use std::hash::{Hash, Hasher};

#[derive(Clone, Copy, Debug, Default)]
pub struct GreedyRunStats {
    pub candidate_queries: u64,
    pub candidate_pairs_scored: u64,
    pub eligible_pairs: u64,
    pub representative_pair_scores: u64,
    pub merge_pair_scores: u64,
}

#[inline]
fn score_pair(nodes: &[Node], scorer: &Scorer, a: u32, b: u32) -> Option<u16> {
    crate::pair_trace::record(a, b);
    scorer.score_scalar(&nodes[a as usize], &nodes[b as usize])
}

#[inline]
fn pair_scores(nodes: &[Node], scorer: &Scorer, a: u32, b: u32) -> Option<SimilarityScores> {
    crate::pair_trace::record(a, b);
    scorer.eligible_pair_scores(&nodes[a as usize], &nodes[b as usize])
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct LazyEntry {
    coverage_upper_bound: u32,
    weight_sum_upper_bound: u64,
    frequency: u64,
    node: u32,
}

impl Ord for LazyEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.coverage_upper_bound
            .cmp(&other.coverage_upper_bound)
            .then(
                self.weight_sum_upper_bound
                    .cmp(&other.weight_sum_upper_bound),
            )
            .then(self.frequency.cmp(&other.frequency))
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for LazyEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn cluster_members(clustering: &Clustering) -> Vec<Vec<u32>> {
    let mut members = vec![Vec::new(); clustering.representatives.len()];
    for (node, &cluster) in clustering.cluster_of.iter().enumerate() {
        members[cluster as usize].push(node as u32);
    }
    members
}

fn state_hash(clustering: &Clustering) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    clustering.cluster_of.hash(&mut hasher);
    clustering.representatives.hash(&mut hasher);
    hasher.finish()
}

pub fn initial_clustering(
    nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
) -> (Clustering, GreedyRunStats) {
    let candidate_counts: Vec<usize> = nodes
        .par_iter()
        .enumerate()
        .map(|(node_id, node)| {
            retrieve_candidates(Some(node_id as u32), node, buckets, table).len()
        })
        .collect();
    let mut order: Vec<usize> = (0..nodes.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        candidate_counts[b]
            .cmp(&candidate_counts[a])
            .then(nodes[a].sequence.cmp(&nodes[b].sequence))
    });

    let mut assigned = vec![false; nodes.len()];
    let mut cluster_of = vec![u32::MAX; nodes.len()];
    let mut representatives = Vec::<u32>::new();
    let mut stats = GreedyRunStats {
        candidate_queries: nodes.len() as u64,
        ..GreedyRunStats::default()
    };

    for rep in order {
        if assigned[rep] {
            continue;
        }
        let cluster = representatives.len() as u32;
        representatives.push(rep as u32);
        assigned[rep] = true;
        cluster_of[rep] = cluster;
        let candidates = retrieve_candidates(Some(rep as u32), &nodes[rep], buckets, table);
        stats.candidate_queries += 1;
        let unassigned: Vec<u32> = candidates
            .into_iter()
            .filter(|candidate| !assigned[*candidate as usize])
            .collect();
        stats.candidate_pairs_scored += unassigned.len() as u64;
        let accepted: Vec<u32> = if unassigned.len() >= 1024 {
            unassigned
                .par_iter()
                .filter_map(|&candidate| {
                    score_pair(nodes, scorer, rep as u32, candidate).map(|_| candidate)
                })
                .collect()
        } else {
            unassigned
                .iter()
                .filter_map(|&candidate| {
                    score_pair(nodes, scorer, rep as u32, candidate).map(|_| candidate)
                })
                .collect()
        };
        stats.eligible_pairs += accepted.len() as u64;
        for candidate in accepted {
            assigned[candidate as usize] = true;
            cluster_of[candidate as usize] = cluster;
        }
    }
    (
        Clustering {
            cluster_of,
            representatives,
        },
        stats,
    )
}

/// Dynamic greedy set cover without materializing the eligible-edge graph.
/// K-mer candidate counts begin as safe coverage upper bounds. A candidate is
/// selected only after exact scoring proves that its current coverage key is
/// no smaller than every remaining upper bound.
pub fn initial_clustering_lazy_exact(
    nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
) -> (Clustering, GreedyRunStats) {
    let candidate_counts: Vec<usize> = nodes
        .par_iter()
        .enumerate()
        .map(|(node_id, node)| {
            retrieve_candidates(Some(node_id as u32), node, buckets, table).len()
        })
        .collect();
    let mut heap = BinaryHeap::with_capacity(nodes.len());
    for (node, &candidate_count) in candidate_counts.iter().enumerate() {
        heap.push(LazyEntry {
            coverage_upper_bound: 1 + candidate_count as u32,
            // Every mode-specific edge ranking weight is at most 1000.
            weight_sum_upper_bound: candidate_count as u64 * 1000,
            frequency: nodes[node].frequency,
            node: node as u32,
        });
    }

    let mut assigned = vec![false; nodes.len()];
    let mut cluster_of = vec![u32::MAX; nodes.len()];
    let mut representatives = Vec::<u32>::new();
    let mut stats = GreedyRunStats {
        candidate_queries: nodes.len() as u64,
        ..GreedyRunStats::default()
    };

    while let Some(entry) = heap.pop() {
        let rep = entry.node as usize;
        if assigned[rep] {
            continue;
        }
        let candidates = retrieve_candidates(Some(entry.node), &nodes[rep], buckets, table);
        stats.candidate_queries += 1;
        let unassigned: Vec<u32> = candidates
            .into_iter()
            .filter(|candidate| !assigned[*candidate as usize])
            .collect();
        stats.candidate_pairs_scored += unassigned.len() as u64;
        let accepted: Vec<(u32, u16)> = if unassigned.len() >= 1024 {
            unassigned
                .par_iter()
                .filter_map(|&candidate| {
                    score_pair(nodes, scorer, rep as u32, candidate)
                        .map(|weight| (candidate, weight))
                })
                .collect()
        } else {
            unassigned
                .iter()
                .filter_map(|&candidate| {
                    score_pair(nodes, scorer, rep as u32, candidate)
                        .map(|weight| (candidate, weight))
                })
                .collect()
        };
        stats.eligible_pairs += accepted.len() as u64;
        let exact = LazyEntry {
            coverage_upper_bound: 1 + accepted.len() as u32,
            weight_sum_upper_bound: accepted.iter().map(|item| item.1 as u64).sum(),
            frequency: nodes[rep].frequency,
            node: entry.node,
        };

        while heap
            .peek()
            .is_some_and(|candidate| assigned[candidate.node as usize])
        {
            heap.pop();
        }
        if heap.peek().is_some_and(|upper_bound| exact < *upper_bound) {
            heap.push(exact);
            continue;
        }

        let cluster = representatives.len() as u32;
        representatives.push(entry.node);
        assigned[rep] = true;
        cluster_of[rep] = cluster;
        for (candidate, _) in accepted {
            assigned[candidate as usize] = true;
            cluster_of[candidate as usize] = cluster;
        }
    }
    debug_assert!(assigned.iter().all(|value| *value));
    (
        Clustering {
            cluster_of,
            representatives,
        },
        stats,
    )
}

fn update_representatives(
    nodes: &[Node],
    scorer: &Scorer,
    clustering: &Clustering,
) -> (Vec<u32>, u64) {
    let members = cluster_members(clustering);
    let results: Vec<(u32, u64)> = members
        .par_iter()
        .enumerate()
        .map(|(cluster_id, cluster)| {
            let mut best = clustering.representatives[cluster_id];
            let mut best_key = (0u128, 0u128, 0u128, nodes[best as usize].frequency);
            let mut comparisons = 0u64;
            for &candidate in cluster {
                let mut aggregate_rank = 0u128;
                let mut aggregate_alignment = 0u128;
                let mut aggregate_anchor = 0u128;
                let mut valid = true;
                for &member in cluster {
                    if member == candidate {
                        continue;
                    }
                    comparisons += 1;
                    let Some(scores) = pair_scores(nodes, scorer, candidate, member) else {
                        valid = false;
                        break;
                    };
                    let frequency = nodes[member as usize].frequency as u128;
                    aggregate_rank += scores.ranking_weight as u128 * frequency;
                    aggregate_alignment += scores.alignment as u128 * frequency;
                    aggregate_anchor += scores.anchor_combination as u128 * frequency;
                }
                if !valid {
                    continue;
                }
                let key = (
                    aggregate_rank,
                    aggregate_alignment,
                    aggregate_anchor,
                    nodes[candidate as usize].frequency,
                );
                if key > best_key || (key == best_key && candidate < best) {
                    best = candidate;
                    best_key = key;
                }
            }
            (best, comparisons)
        })
        .collect();
    (
        results.iter().map(|result| result.0).collect(),
        results.iter().map(|result| result.1).sum(),
    )
}

fn reassign(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    clustering: &Clustering,
) -> (Vec<u32>, u64, u64) {
    let rep_buckets = build_exact_index_subset(
        nodes,
        clustering.representatives.iter().map(|rep| *rep as usize),
    );
    let mut rep_cluster = vec![u32::MAX; nodes.len()];
    for (cluster, &rep) in clustering.representatives.iter().enumerate() {
        rep_cluster[rep as usize] = cluster as u32;
    }
    let proposals: Vec<(u32, u64, u64)> = nodes
        .par_iter()
        .enumerate()
        .map(|(node_id, node)| {
            let current_cluster = clustering.cluster_of[node_id];
            let current_rep = clustering.representatives[current_cluster as usize];
            let mut candidates =
                retrieve_candidates(Some(node_id as u32), node, &rep_buckets, table);
            candidates.push(current_rep);
            if rep_cluster[node_id] != u32::MAX {
                candidates.push(node_id as u32);
            }
            candidates.sort_unstable();
            candidates.dedup();
            let mut best_rep = current_rep;
            let mut best_cluster = current_cluster;
            let mut best_scores = pair_scores(nodes, scorer, node_id as u32, current_rep)
                .expect("current representative must remain eligible");
            let mut scored = 0u64;
            let mut eligible = 0u64;
            for rep in candidates {
                let cluster = rep_cluster[rep as usize];
                if cluster == u32::MAX {
                    continue;
                }
                scored += 1;
                let Some(scores) = pair_scores(nodes, scorer, node_id as u32, rep) else {
                    continue;
                };
                eligible += 1;
                if scorer.compare_scores(scores, best_scores).is_gt()
                    || (scorer.compare_scores(scores, best_scores).is_eq() && rep < best_rep)
                {
                    best_rep = rep;
                    best_cluster = cluster;
                    best_scores = scores;
                }
            }
            (best_cluster, scored, eligible)
        })
        .collect();
    (
        proposals.iter().map(|proposal| proposal.0).collect(),
        proposals.iter().map(|proposal| proposal.1).sum(),
        proposals.iter().map(|proposal| proposal.2).sum(),
    )
}

fn compact_empty_clusters(mut clustering: Clustering) -> Clustering {
    let mut used = vec![false; clustering.representatives.len()];
    for &cluster in &clustering.cluster_of {
        used[cluster as usize] = true;
    }
    let mut remap = vec![u32::MAX; used.len()];
    let mut representatives = Vec::new();
    for (old, active) in used.into_iter().enumerate() {
        if active {
            remap[old] = representatives.len() as u32;
            representatives.push(clustering.representatives[old]);
        }
    }
    for cluster in &mut clustering.cluster_of {
        *cluster = remap[*cluster as usize];
    }
    clustering.representatives = representatives;
    clustering
}

fn strict_merge(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    mut clustering: Clustering,
    merge_cap: Option<usize>,
) -> (Clustering, usize, u64) {
    let mut members = cluster_members(&clustering);
    let mut active = vec![true; members.len()];
    let rep_buckets = build_exact_index_subset(
        nodes,
        clustering.representatives.iter().map(|rep| *rep as usize),
    );
    let mut rep_cluster = vec![u32::MAX; nodes.len()];
    for (cluster, &rep) in clustering.representatives.iter().enumerate() {
        rep_cluster[rep as usize] = cluster as u32;
    }
    let mut order: Vec<usize> = (0..members.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        members[b]
            .len()
            .cmp(&members[a].len())
            .then(clustering.representatives[a].cmp(&clustering.representatives[b]))
    });
    let mut merges = 0usize;
    let mut comparisons = 0u64;
    for &target in &order {
        if !active[target] {
            continue;
        }
        let target_rep = clustering.representatives[target];
        let candidates = retrieve_candidates(
            Some(target_rep),
            &nodes[target_rep as usize],
            &rep_buckets,
            table,
        );
        let mut sources: Vec<usize> = candidates
            .into_iter()
            .filter_map(|rep| {
                let cluster = rep_cluster[rep as usize];
                (cluster != u32::MAX).then_some(cluster as usize)
            })
            .filter(|&source| {
                source != target && active[source] && members[target].len() >= members[source].len()
            })
            .collect();
        sources.sort_unstable_by(|&a, &b| {
            members[b]
                .len()
                .cmp(&members[a].len())
                .then(clustering.representatives[a].cmp(&clustering.representatives[b]))
        });
        sources.dedup();
        for source in sources {
            if !active[source] || members[target].len() < members[source].len() {
                continue;
            }
            comparisons += 1;
            if score_pair(
                nodes,
                scorer,
                target_rep,
                clustering.representatives[source],
            )
            .is_none()
            {
                continue;
            }
            let source_rep = clustering.representatives[source];
            let early_failure = merge_cap.is_some_and(|cap| {
                members[source]
                    .iter()
                    .filter(|member| **member != source_rep)
                    .take(cap)
                    .any(|member| {
                        comparisons += 1;
                        score_pair(nodes, scorer, target_rep, *member).is_none()
                    })
            });
            if early_failure {
                continue;
            }
            let mut valid = true;
            for member in members[source]
                .iter()
                .filter(|member| **member != source_rep)
            {
                comparisons += 1;
                if score_pair(nodes, scorer, target_rep, *member).is_none() {
                    valid = false;
                    break;
                }
            }
            if !valid {
                continue;
            }
            let source_members = std::mem::take(&mut members[source]);
            for member in &source_members {
                clustering.cluster_of[*member as usize] = target as u32;
            }
            members[target].extend(source_members);
            active[source] = false;
            merges += 1;
        }
    }
    (compact_empty_clusters(clustering), merges, comparisons)
}

pub fn iterate_to_convergence(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    mut clustering: Clustering,
    iteration_cap: Option<usize>,
    merge: bool,
    merge_cap: Option<usize>,
    stats: &mut GreedyRunStats,
) -> (Clustering, IterationStats) {
    let mut iteration = IterationStats::default();
    let mut seen = HashSet::new();
    loop {
        if iteration_cap.is_some_and(|cap| iteration.iterations >= cap)
            || !seen.insert(state_hash(&clustering))
        {
            break;
        }
        iteration.iterations += 1;
        let old_assignments = clustering.cluster_of.clone();
        let old_representatives = clustering.representatives.clone();
        let (assignments, scored, eligible) = reassign(nodes, table, scorer, &clustering);
        stats.candidate_pairs_scored += scored;
        stats.eligible_pairs += eligible;
        iteration.reassignment_moves += assignments
            .iter()
            .zip(&clustering.cluster_of)
            .filter(|(a, b)| a != b)
            .count() as u64;
        clustering.cluster_of = assignments;
        clustering = compact_empty_clusters(clustering);
        let (representatives, comparisons) = update_representatives(nodes, scorer, &clustering);
        stats.representative_pair_scores += comparisons;
        iteration.representative_changes += representatives
            .iter()
            .zip(&clustering.representatives)
            .filter(|(a, b)| a != b)
            .count() as u64;
        clustering.representatives = representatives;

        let mut merges = 0usize;
        if merge {
            let result = strict_merge(nodes, table, scorer, clustering, merge_cap);
            clustering = result.0;
            merges = result.1;
            stats.merge_pair_scores += result.2;
            if merges > 0 {
                let (representatives, comparisons) =
                    update_representatives(nodes, scorer, &clustering);
                stats.representative_pair_scores += comparisons;
                iteration.representative_changes += representatives
                    .iter()
                    .zip(&clustering.representatives)
                    .filter(|(a, b)| a != b)
                    .count() as u64;
                clustering.representatives = representatives;
            }
        }
        iteration.merges += merges;

        let failures = clustering
            .cluster_of
            .iter()
            .enumerate()
            .filter(|(node, cluster)| {
                let rep = clustering.representatives[**cluster as usize];
                score_pair(nodes, scorer, *node as u32, rep).is_none()
            })
            .count() as u64;
        iteration.validation_failures += failures;
        if failures > 0 {
            break;
        }
        if clustering.cluster_of == old_assignments
            && clustering.representatives == old_representatives
            && merges == 0
        {
            iteration.converged = true;
            break;
        }
    }
    (canonicalize(clustering), iteration)
}

pub fn representative_coverage(
    nodes: &[Node],
    scorer: &Scorer,
    clustering: &Clustering,
) -> (u64, u64) {
    let covered = clustering
        .cluster_of
        .iter()
        .enumerate()
        .filter(|(node, cluster)| {
            let rep = clustering.representatives[**cluster as usize];
            score_pair(nodes, scorer, *node as u32, rep).is_some()
        })
        .count() as u64;
    (covered, nodes.len() as u64)
}
