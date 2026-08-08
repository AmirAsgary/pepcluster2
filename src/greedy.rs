use crate::graph::{canonicalize, intrinsic_order, Clustering, IterationStats, RepresentativeOrder};
use crate::index::{build_exact_index_subset, retrieve_candidates, TerminalSeed};
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
    pub index_candidate_occurrences: u64,
    pub anchor_bound_rejected: u64,
    pub candidate_pairs_scored: u64,
    pub eligible_pairs: u64,
    pub representative_pair_scores: u64,
    pub merge_pair_scores: u64,
    /// Reinserted nodes that reused a cached candidate list instead of
    /// re-retrieving and rescoring it.
    pub cache_hits: u64,
    /// Lists not retained because the cache budget was full. Costs time on a
    /// later pop, never correctness.
    pub cache_evictions: u64,
    /// High-water mark of cached pairs, for sizing --max-memory-gb.
    pub cache_peak_pairs: u64,
}

/// Index retrieval followed by the sound anchor upper bound. Both greedy
/// selections and the reassignment pass share this so their candidate lists are
/// pruned exactly like the graph path's.
fn bounded_candidates(
    query: u32,
    nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
) -> (Vec<u32>, u64) {
    let retrieved = retrieve_candidates(Some(query), &nodes[query as usize], buckets, table, seed);
    let total = retrieved.len() as u64;
    let kept: Vec<u32> = retrieved
        .into_iter()
        .filter(|candidate| {
            scorer.anchor_bound_passes(&nodes[query as usize], &nodes[*candidate as usize])
        })
        .collect();
    (kept, total)
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

#[allow(clippy::too_many_arguments)]
pub fn initial_clustering(
    nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
    representative_order: RepresentativeOrder,
) -> (Clustering, GreedyRunStats) {
    let mut stats = GreedyRunStats {
        candidate_queries: nodes.len() as u64,
        ..GreedyRunStats::default()
    };
    let order: Vec<usize> = match representative_order {
        RepresentativeOrder::Intrinsic => intrinsic_order(nodes)
            .into_iter()
            .map(|id| id as usize)
            .collect(),
        RepresentativeOrder::Coverage => {
            // Degree proxy: count only candidates that survive the sound anchor
            // bound, so the ordering tracks the eligible-neighbour count instead
            // of raw index traffic.
            let counts: Vec<(usize, u64)> = nodes
                .par_iter()
                .enumerate()
                .map(|(node_id, _)| {
                    let (kept, retrieved) =
                        bounded_candidates(node_id as u32, nodes, buckets, table, scorer, seed);
                    (kept.len(), retrieved)
                })
                .collect();
            stats.index_candidate_occurrences += counts.iter().map(|item| item.1).sum::<u64>();
            stats.anchor_bound_rejected += counts
                .iter()
                .map(|item| item.1 - item.0 as u64)
                .sum::<u64>();
            let candidate_counts: Vec<usize> = counts.iter().map(|item| item.0).collect();
            let mut order: Vec<usize> = (0..nodes.len()).collect();
            order.sort_unstable_by(|&a, &b| {
                candidate_counts[b]
                    .cmp(&candidate_counts[a])
                    .then(nodes[a].sequence.cmp(&nodes[b].sequence))
            });
            order
        }
    };

    let mut assigned = vec![false; nodes.len()];
    let mut cluster_of = vec![u32::MAX; nodes.len()];
    let mut representatives = Vec::<u32>::new();

    for rep in order {
        if assigned[rep] {
            continue;
        }
        let cluster = representatives.len() as u32;
        representatives.push(rep as u32);
        assigned[rep] = true;
        cluster_of[rep] = cluster;
        let (candidates, retrieved) =
            bounded_candidates(rep as u32, nodes, buckets, table, scorer, seed);
        stats.candidate_queries += 1;
        stats.index_candidate_occurrences += retrieved;
        stats.anchor_bound_rejected += retrieved - candidates.len() as u64;
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
    seed: TerminalSeed,
    cache_budget_bytes: u64,
) -> (Clustering, GreedyRunStats) {
    let mut stats = GreedyRunStats {
        candidate_queries: nodes.len() as u64,
        ..GreedyRunStats::default()
    };
    let bounds: Vec<(usize, u64)> = nodes
        .par_iter()
        .enumerate()
        .map(|(node_id, _)| {
            let (kept, retrieved) =
                bounded_candidates(node_id as u32, nodes, buckets, table, scorer, seed);
            (kept.len(), retrieved)
        })
        .collect();
    stats.index_candidate_occurrences += bounds.iter().map(|item| item.1).sum::<u64>();
    stats.anchor_bound_rejected += bounds
        .iter()
        .map(|item| item.1 - item.0 as u64)
        .sum::<u64>();
    let mut heap = BinaryHeap::with_capacity(nodes.len());
    for (node, &(candidate_count, _)) in bounds.iter().enumerate() {
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

    // Scored candidate lists for nodes waiting in the heap after reinsertion.
    //
    // A node's candidate set, and the score of every pair it takes part in, are
    // functions of immutable data; `assigned` only ever gains members. The
    // accepted list at a later pop is therefore exactly the list from the
    // earlier pop with newly assigned candidates removed - no retrieval and no
    // rescoring. Without this a reinserted node repeats its whole retrieval and
    // scoring pass, measured as 14x more pairs scored and 3.5x more index
    // retrievals than the static ordering needs.
    //
    // Only reinserted nodes are held. An entry is dropped as soon as its node
    // becomes a representative or is found already assigned, so the cache never
    // grows into the full edge graph that this path exists to avoid.
    let mut cache: Vec<Option<Vec<(u32, u16)>>> = vec![None; nodes.len()];
    // The cache is bounded. Retained lists grow with how many nodes are waiting
    // in the heap after reinsertion and with how many candidates each accepts,
    // and a low acceptance threshold inflates both: measured at one million
    // peptides, peak resident memory was 20.5 GB at threshold 0.40 and 79.5 GB
    // at 0.25. Past the budget a list is simply not retained and the next pop
    // recomputes it, so the budget trades time for memory and never changes the
    // result.
    let budget_pairs = cache_budget_bytes / std::mem::size_of::<(u32, u16)>() as u64;
    let mut cached_pairs: u64 = 0;

    while let Some(entry) = heap.pop() {
        let rep = entry.node as usize;
        if assigned[rep] {
            if let Some(stale) = cache[rep].take() {
                cached_pairs -= stale.len() as u64;
            }
            continue;
        }
        // Cached path: reuse the earlier scoring, keeping only candidates that
        // are still unassigned.
        if let Some(previous) = cache[rep].take() {
            cached_pairs -= previous.len() as u64;
            stats.cache_hits += 1;
            let accepted: Vec<(u32, u16)> = previous
                .into_iter()
                .filter(|(candidate, _)| !assigned[*candidate as usize])
                .collect();
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
                if cached_pairs + accepted.len() as u64 <= budget_pairs {
                    cached_pairs += accepted.len() as u64;
                    stats.cache_peak_pairs = stats.cache_peak_pairs.max(cached_pairs);
                    cache[rep] = Some(accepted);
                } else {
                    stats.cache_evictions += 1;
                }
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
            continue;
        }
        // Nothing cached for this node, so it needs a full retrieval and scoring
        // pass. Do it for a batch of heap candidates at once rather than one at
        // a time.
        //
        // This is exact. An exact value computed against the current `assigned`
        // snapshot stays a valid upper bound for every later state, because
        // assigning more peptides can only remove candidates from a node's
        // accepted list. Batch members therefore go back on the heap carrying
        // admissible bounds, and the commit rule - take a node only when its
        // exact value is at least every remaining upper bound - still selects
        // the node with the greatest true exact value. Bounds only ever tighten,
        // so a batched entry can move down the heap but never up.
        //
        // The main loop stays sequential; only the retrieval and scoring, which
        // is where the time goes, runs in parallel.
        // Batch only as far as the cache can still hold the results. Batching
        // pays for itself by amortising a node's evaluation across the later
        // pops that reuse it; once the budget is exhausted nothing is retained,
        // so a wide batch would re-evaluate dozens of nodes on every pop and get
        // slower the tighter memory is. Falling back to one node per pop makes an
        // exhausted cache degrade to the original serial behaviour instead.
        let batch_size = if budget_pairs.saturating_sub(cached_pairs) == 0 {
            1
        } else {
            rayon::current_num_threads().max(1) * 4
        };
        let mut batch: Vec<LazyEntry> = vec![entry];
        let mut deferred: Vec<LazyEntry> = Vec::new();
        while batch.len() < batch_size {
            match heap.pop() {
                None => break,
                Some(next) => {
                    if assigned[next.node as usize] {
                        if let Some(stale) = cache[next.node as usize].take() {
                            cached_pairs -= stale.len() as u64;
                        }
                        continue;
                    }
                    if cache[next.node as usize].is_some() {
                        // Already cheap to evaluate; leave it for the serial path.
                        deferred.push(next);
                        break;
                    }
                    batch.push(next);
                }
            }
        }
        let evaluated: Vec<(LazyEntry, Vec<(u32, u16)>, u64, u64, u64)> = batch
            .par_iter()
            .map(|member| {
                let node = member.node;
                // Cheapest test first. `bounded_candidates` applies the anchor
                // upper bound - 36 table lookups per pair - to everything the
                // index returns, including candidates that are already assigned
                // and will be dropped immediately afterwards. Discarding those
                // first skips the bound entirely for them. Both filters are
                // independent of each other, so the surviving set is unchanged.
                let retrieved_list =
                    retrieve_candidates(Some(node), &nodes[node as usize], buckets, table, seed);
                let retrieved = retrieved_list.len() as u64;
                let live: Vec<u32> = retrieved_list
                    .into_iter()
                    .filter(|candidate| !assigned[*candidate as usize])
                    .collect();
                let unassigned: Vec<u32> = live
                    .into_iter()
                    .filter(|candidate| {
                        scorer.anchor_bound_passes(
                            &nodes[node as usize],
                            &nodes[*candidate as usize],
                        )
                    })
                    .collect();
                let scored = unassigned.len() as u64;
                let rejected = retrieved - scored;
                let accepted: Vec<(u32, u16)> = unassigned
                    .into_iter()
                    .filter_map(|candidate| {
                        score_pair(nodes, scorer, node, candidate)
                            .map(|weight| (candidate, weight))
                    })
                    .collect();
                let exact = LazyEntry {
                    coverage_upper_bound: 1 + accepted.len() as u32,
                    weight_sum_upper_bound: accepted.iter().map(|item| item.1 as u64).sum(),
                    frequency: nodes[node as usize].frequency,
                    node,
                };
                (exact, accepted, retrieved, rejected, scored)
            })
            .collect();
        // The popped node is resolved here, not by pushing it back and waiting
        // for a cache hit to take the commit branch. Relying on the cache
        // deadlocks the moment the budget is exhausted: nothing is retained, so
        // the node is re-evaluated and pushed back forever and the loop never
        // makes progress. Every batch member other than the popped one goes back
        // on the heap carrying its exact value, which is an admissible bound.
        let mut own: Option<Vec<(u32, u16)>> = None;
        for (exact, accepted, retrieved, rejected, scored) in evaluated {
            stats.candidate_queries += 1;
            stats.index_candidate_occurrences += retrieved;
            stats.anchor_bound_rejected += rejected;
            stats.candidate_pairs_scored += scored;
            stats.eligible_pairs += accepted.len() as u64;
            if exact.node as usize == rep {
                own = Some(accepted);
                heap.push(exact);
                continue;
            }
            if cached_pairs + accepted.len() as u64 <= budget_pairs {
                cached_pairs += accepted.len() as u64;
                stats.cache_peak_pairs = stats.cache_peak_pairs.max(cached_pairs);
                cache[exact.node as usize] = Some(accepted);
            } else {
                stats.cache_evictions += 1;
            }
            heap.push(exact);
        }
        for entry in deferred {
            heap.push(entry);
        }
        let Some(accepted) = own else { continue };
        // Retake the entry just pushed for this node so the commit test below
        // compares it against the rest of the heap, exactly as the serial code
        // did before batching.
        let exact = match heap.pop() {
            Some(top) if top.node as usize == rep => top,
            Some(other) => {
                heap.push(other);
                if cached_pairs + accepted.len() as u64 <= budget_pairs {
                    cached_pairs += accepted.len() as u64;
                    stats.cache_peak_pairs = stats.cache_peak_pairs.max(cached_pairs);
                    cache[rep] = Some(accepted);
                }
                continue;
            }
            None => continue,
        };
        while heap
            .peek()
            .is_some_and(|candidate| assigned[candidate.node as usize])
        {
            heap.pop();
        }
        if heap.peek().is_some_and(|upper_bound| exact < *upper_bound) {
            if cached_pairs + accepted.len() as u64 <= budget_pairs {
                cached_pairs += accepted.len() as u64;
                stats.cache_peak_pairs = stats.cache_peak_pairs.max(cached_pairs);
                cache[rep] = Some(accepted);
            } else {
                stats.cache_evictions += 1;
            }
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

/// Synchronous reassignment against the representative index. The
/// `minimum_improvement` hysteresis matches `graph::refine`: a peptide leaves
/// its current representative only when another beats it by more than that
/// margin.
fn reassign(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    clustering: &Clustering,
    seed: TerminalSeed,
    minimum_improvement: u16,
) -> (Vec<u32>, u64, u64) {
    let rep_buckets = build_exact_index_subset(
        nodes,
        clustering.representatives.iter().map(|rep| *rep as usize),
        seed,
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
                retrieve_candidates(Some(node_id as u32), node, &rep_buckets, table, seed);
            candidates.retain(|rep| {
                scorer.anchor_bound_passes(node, &nodes[*rep as usize])
            });
            candidates.push(current_rep);
            if rep_cluster[node_id] != u32::MAX {
                candidates.push(node_id as u32);
            }
            candidates.sort_unstable();
            candidates.dedup();
            let mut best_rep = current_rep;
            let mut best_cluster = current_cluster;
            let current_scores = pair_scores(nodes, scorer, node_id as u32, current_rep)
                .expect("current representative must remain eligible");
            let mut best_scores = current_scores;
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
            // Hysteresis: keep the current representative unless the winner
            // beats it by more than the margin. Applied only to the decision to
            // leave, so ties among the alternatives still resolve
            // deterministically above.
            if best_rep != current_rep
                && best_scores.ranking_weight as u32
                    <= current_scores.ranking_weight as u32 + minimum_improvement as u32
            {
                best_cluster = current_cluster;
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
    seed: TerminalSeed,
) -> (Clustering, usize, u64) {
    let mut members = cluster_members(&clustering);
    let mut active = vec![true; members.len()];
    let rep_buckets = build_exact_index_subset(
        nodes,
        clustering.representatives.iter().map(|rep| *rep as usize),
        seed,
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
            seed,
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

#[allow(clippy::too_many_arguments)]
pub fn iterate_to_convergence(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    mut clustering: Clustering,
    iteration_cap: Option<usize>,
    merge: bool,
    merge_cap: Option<usize>,
    seed: TerminalSeed,
    minimum_improvement: u16,
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
        let (assignments, scored, eligible) =
            reassign(nodes, table, scorer, &clustering, seed, minimum_improvement);
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
            let result = strict_merge(nodes, table, scorer, clustering, merge_cap, seed);
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
