//! Optional motif-level layer above the similarity clustering.
//!
//! # Why this stage exists
//!
//! A PepCluster2 cluster is a similarity ball: every member passes the selected
//! scoring rule against its representative. A binding motif is a different
//! object. It is a product of per-position residue preferences, narrow at the
//! anchor positions and close to flat everywhere else, so the region of sequence
//! space it occupies is strongly anisotropic. A ball in an additive similarity
//! cannot cover such a region, and lowering the threshold widens the ball along
//! every axis at once rather than along the tolerant ones only. One motif
//! therefore fragments into many similarity clusters, and no single threshold
//! repairs that.
//!
//! This module reads the finished similarity clusters as amino-acid profiles and
//! merges the pairs for which a Dirichlet-multinomial model prefers one shared
//! profile to two separate ones, optionally followed by expectation-maximization
//! refinement of a mixture of position weight matrices.
//!
//! # What it does not do
//!
//! The motif partition does **not** satisfy the representative-to-member
//! invariant of Section 13.4 of the algorithm specification: two peptides in one
//! motif need not pass the scoring rule against any common representative, which
//! is exactly the point. It is therefore reported as a separate output layer and
//! never overwrites the similarity clustering.
//!
//! Merging can only coarsen a partition. It cannot move a peptide out of a
//! cluster it should not have been in, so contamination present in the input
//! fragments propagates. The EM stage can move peptides and is the only part
//! that repairs such errors.
//!
//! # Determinism
//!
//! Every quantity is `f64`, accumulated in a fixed order. The agglomeration
//! argmax breaks ties by cluster index, EM is seeded from the merge result
//! rather than at random, and the parallel accumulation uses fixed chunk
//! boundaries combined in index order. Results do not depend on thread count or
//! on FASTA record order.

use crate::graph::Clustering;
use crate::model::Node;
use rayon::prelude::*;

/// Columns in the motif frame.
pub const MOTIF_COLUMNS: usize = 9;
/// Canonical amino acids.
pub const ALPHABET: usize = 20;
/// Frame column with no residue, used for the central column of an 8-mer.
pub const MISSING: u8 = 255;

/// Natural logarithm of the gamma function (Lanczos, g = 7, n = 9).
///
/// Accurate to roughly 15 significant digits over the positive reals, which is
/// what the marginal likelihood needs: the Bayes factor is a difference of three
/// such values and cancellation would otherwise dominate for large clusters.
pub fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection formula. Not reached by this module's own callers, which
        // only ever pass positive concentrations and non-negative counts, but
        // keeping it makes the function correct in isolation.
        std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().abs().ln()
            - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let t = x + 7.5;
        let mut a = C[0];
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Project a peptide onto the nine motif columns.
///
/// ```text
/// L >= 9   columns 1-4 <- positions 1..4        columns 5-9 <- positions L-4..L
/// L == 8   columns 1-4 <- positions 1..4        column  5   <- gap
///                                              columns 6-9 <- positions 5..8
/// ```
///
/// For a 9-mer this is the identity, so the dominant MHC-I length loses nothing.
/// For longer peptides the central residues are dropped: they bulge out of the
/// binding groove, make little contact with the MHC, and carry correspondingly
/// little allele-specific signal. An 8-mer leaves the central column unobserved
/// rather than shifting its C-terminal residues, which would put the dominant
/// C-terminal anchor in the wrong column.
///
/// Returns `MISSING` in a column with no residue, and for any non-canonical
/// residue code.
pub fn frame(sequence_codes: &[u8]) -> [u8; MOTIF_COLUMNS] {
    let mut out = [MISSING; MOTIF_COLUMNS];
    let n = sequence_codes.len();
    let take = |code: u8| if (code as usize) < ALPHABET { code } else { MISSING };
    if n >= MOTIF_COLUMNS {
        for c in 0..4 {
            out[c] = take(sequence_codes[c]);
        }
        for (offset, c) in (4..MOTIF_COLUMNS).enumerate() {
            out[c] = take(sequence_codes[n - 5 + offset]);
        }
    } else if n == 8 {
        for c in 0..4 {
            out[c] = take(sequence_codes[c]);
        }
        for (offset, c) in (5..MOTIF_COLUMNS).enumerate() {
            out[c] = take(sequence_codes[4 + offset]);
        }
    }
    out
}

/// Amino-acid counts for one cluster, laid out column-major as
/// `column * ALPHABET + residue`.
type Profile = [f64; MOTIF_COLUMNS * ALPHABET];

fn empty_profile() -> Profile {
    [0.0; MOTIF_COLUMNS * ALPHABET]
}

/// Dirichlet prior shared by every column, and the constants its marginal
/// likelihood needs.
struct Prior {
    alpha: [f64; ALPHABET],
    /// `sum_a alpha_a`.
    total: f64,
    /// `ln Gamma(total)`.
    ln_gamma_total: f64,
    /// `sum_a ln Gamma(alpha_a)`.
    ln_gamma_alpha_sum: f64,
}

impl Prior {
    /// `concentration` sets the total pseudocount per column: the prior is
    /// `alpha_a = concentration * ALPHABET * background_a`, so `concentration =
    /// 1` puts one pseudocount per residue on average and larger values smooth
    /// harder. Spreading over the observed background rather than uniformly
    /// stops a rare residue from being treated as equally expected everywhere.
    fn new(concentration: f64, background: &[f64; ALPHABET]) -> Self {
        let mut alpha = [0.0; ALPHABET];
        for a in 0..ALPHABET {
            // The floor keeps ln Gamma finite when a residue is absent from the
            // dataset entirely.
            alpha[a] = (concentration * ALPHABET as f64 * background[a]).max(1e-9);
        }
        let total = alpha.iter().sum::<f64>();
        Self {
            ln_gamma_total: ln_gamma(total),
            ln_gamma_alpha_sum: alpha.iter().map(|&a| ln_gamma(a)).sum(),
            alpha,
            total,
        }
    }

    /// Log marginal likelihood of the labelled residues behind `counts` under a
    /// single profile, with the per-column amino-acid distribution integrated
    /// out rather than fitted:
    ///
    /// ```text
    /// log L = sum_j [ ln G(A0) - ln G(A0 + N_j)
    ///                 + sum_a ( ln G(n_ja + alpha_a) - ln G(alpha_a) ) ]
    /// ```
    ///
    /// There is no multinomial coefficient because the data are the labelled
    /// observations, not the unordered counts. Including one would not cancel in
    /// the Bayes factor and would silently change the criterion.
    fn log_marginal(&self, counts: &Profile) -> f64 {
        let mut total = 0.0;
        for column in 0..MOTIF_COLUMNS {
            let slice = &counts[column * ALPHABET..(column + 1) * ALPHABET];
            let observed: f64 = slice.iter().sum();
            let mut term = self.ln_gamma_total - ln_gamma(self.total + observed);
            for a in 0..ALPHABET {
                term += ln_gamma(slice[a] + self.alpha[a]);
            }
            total += term - self.ln_gamma_alpha_sum;
        }
        total
    }
}

/// Tunable parameters of the motif stage.
#[derive(Clone, Copy, Debug)]
pub struct MotifParams {
    /// Dirichlet prior concentration used when scoring merges.
    pub prior_concentration: f64,
    /// Stop merging when exactly this many groups remain, overriding
    /// `merge_threshold`. In practice the alleles in a sample are usually known
    /// from typing, so supplying the count is ordinary use rather than an oracle.
    /// EM may still empty a component afterwards, so the final motif count can be
    /// lower.
    pub target_count: Option<usize>,
    /// Merge while the best available log Bayes factor exceeds this. Zero means
    /// "merge whenever the evidence favours one profile over two". A positive
    /// value is a flat per-cluster penalty: requiring `log BF > t` is the same
    /// as a prior over partitions proportional to `exp(-t * k)`. Because the
    /// same constant applies to every candidate pair it shifts where the
    /// agglomeration stops without changing the order in which pairs merge.
    pub merge_threshold: f64,
    /// Run EM refinement after merging.
    pub em: bool,
    /// Dirichlet prior concentration used to smooth the EM profiles.
    pub em_prior_concentration: f64,
    pub em_max_iterations: usize,
    /// Relative change in log likelihood below which EM stops.
    pub em_tolerance: f64,
}

/// Outcome of the motif stage.
pub struct MotifResult {
    /// Motif index per node, parallel to `Clustering::cluster_of`.
    pub motif_of: Vec<u32>,
    /// Motifs holding at least one peptide.
    pub motif_count: usize,
    /// Motifs after merging, before EM.
    pub merged_count: usize,
    pub merges: usize,
    pub em_iterations: usize,
    pub em_converged: bool,
    /// Final probability profiles, one per merged motif, laid out
    /// `column * ALPHABET + residue`. Present whether or not EM ran.
    pub profiles: Vec<Vec<f64>>,
    /// Peptide count per merged motif, weighted by input frequency.
    pub weights: Vec<f64>,
}

/// Upper-triangle index for `i < j` over `k` items.
#[inline]
fn tri(i: usize, j: usize, k: usize) -> usize {
    i * k - i * (i + 1) / 2 + (j - i - 1)
}

/// Background amino-acid frequencies over the whole dataset, frequency-weighted.
fn background(nodes: &[Node], frames: &[[u8; MOTIF_COLUMNS]]) -> [f64; ALPHABET] {
    let mut counts = [0.0f64; ALPHABET];
    for (node, frame) in nodes.iter().zip(frames) {
        let weight = node.frequency as f64;
        for &residue in frame.iter() {
            if (residue as usize) < ALPHABET {
                counts[residue as usize] += weight;
            }
        }
    }
    let total: f64 = counts.iter().sum();
    if total <= 0.0 {
        return [1.0 / ALPHABET as f64; ALPHABET];
    }
    for c in counts.iter_mut() {
        // A uniform floor keeps the prior proper for residues absent from the
        // dataset; without it their pseudocount would be zero and any later
        // observation would produce an infinite log density.
        *c = (*c / total).max(1e-6);
    }
    let renormalise: f64 = counts.iter().sum();
    for c in counts.iter_mut() {
        *c /= renormalise;
    }
    counts
}

#[inline]
fn merged_profile(a: &Profile, b: &Profile) -> Profile {
    let mut out = *a;
    for (target, value) in out.iter_mut().zip(b.iter()) {
        *target += value;
    }
    out
}

/// Greedy agglomeration by log Bayes factor.
///
/// Returns the group index of each input cluster and the resulting group count.
///
/// Cost is `O(K^2)` marginal-likelihood evaluations to fill the score matrix,
/// then `O(K)` per merge to refresh the row of the surviving cluster. A cached
/// best partner per row keeps the argmax at `O(K)` instead of `O(K^2)`; a row is
/// recomputed only when its cached partner has just been consumed by a merge.
/// Memory is the upper triangle of the score matrix, `K(K-1)/2` doubles.
///
/// Every cost here scales with the number of clusters, not the number of
/// peptides. The peptides are read once, when the profiles are built.
fn agglomerate(profiles: &mut [Profile], prior: &Prior, threshold: f64,
               target: Option<usize>) -> (Vec<u32>, usize) {
    let k = profiles.len();
    let mut group: Vec<u32> = (0..k as u32).collect();
    if k < 2 {
        return (group, k);
    }
    let mut remaining = k;
    let mut single: Vec<f64> = profiles.iter().map(|p| prior.log_marginal(p)).collect();
    let mut alive = vec![true; k];

    // Row-wise so each task owns a contiguous span and the layout is fixed.
    let rows: Vec<Vec<f64>> = (0..k)
        .into_par_iter()
        .map(|i| {
            let source = &profiles[i];
            ((i + 1)..k)
                .map(|j| {
                    prior.log_marginal(&merged_profile(source, &profiles[j]))
                        - single[i]
                        - single[j]
                })
                .collect()
        })
        .collect();
    let mut score = Vec::with_capacity(k * (k - 1) / 2);
    for row in rows {
        score.extend_from_slice(&row);
    }

    // best[i] = (score, partner) over live j > i. Ties resolve to the smallest
    // partner index, so the merge sequence is reproducible.
    let refresh = |score: &[f64], alive: &[bool], i: usize| -> (f64, usize) {
        let mut best = f64::NEG_INFINITY;
        let mut partner = usize::MAX;
        for j in (i + 1)..k {
            if alive[j] && score[tri(i, j, k)] > best {
                best = score[tri(i, j, k)];
                partner = j;
            }
        }
        (best, partner)
    };
    let mut best: Vec<(f64, usize)> = (0..k).map(|i| refresh(&score, &alive, i)).collect();

    loop {
        let mut top = f64::NEG_INFINITY;
        let mut pair = None;
        for i in 0..k {
            if alive[i] && best[i].1 != usize::MAX && best[i].0 > top {
                top = best[i].0;
                pair = Some((i, best[i].1));
            }
        }
        let Some((i, j)) = pair else { break };
        // A requested count overrides the evidence threshold: the merge order is
        // still the Bayes factor, only the stopping rule changes.
        match target {
            Some(want) if remaining <= want.max(1) => break,
            Some(_) => {}
            None if !(top > threshold) => break,
            None => {}
        }

        profiles[i] = merged_profile(&profiles[i], &profiles[j]);
        single[i] = prior.log_marginal(&profiles[i]);
        alive[j] = false;
        for g in group.iter_mut() {
            if *g == j as u32 {
                *g = i as u32;
            }
        }
        remaining -= 1;

        // Refresh every score touching the surviving cluster.
        let base = profiles[i];
        let single_i = single[i];
        let updates: Vec<(usize, f64)> = (0..k)
            .into_par_iter()
            .filter(|&other| other != i && alive[other])
            .map(|other| {
                (
                    other,
                    prior.log_marginal(&merged_profile(&base, &profiles[other])) - single_i
                        - single[other],
                )
            })
            .collect();
        for &(other, value) in &updates {
            let (a, b) = if other < i { (other, i) } else { (i, other) };
            score[tri(a, b, k)] = value;
        }

        // A cached partner is stale only if it was just consumed, or if the
        // refreshed score against i now beats it.
        best[i] = refresh(&score, &alive, i);
        for other in 0..k {
            if other == i || !alive[other] {
                continue;
            }
            if best[other].1 == i || best[other].1 == j {
                best[other] = refresh(&score, &alive, other);
            } else if other < i {
                let value = score[tri(other, i, k)];
                if value > best[other].0 {
                    best[other] = (value, i);
                }
            }
        }
    }

    // Compact group labels to a dense 0..merged range, in first-appearance order
    // of the surviving cluster index so the labelling is deterministic.
    let mut remap = vec![u32::MAX; k];
    let mut next = 0u32;
    for g in group.iter_mut() {
        let root = *g as usize;
        if remap[root] == u32::MAX {
            remap[root] = next;
            next += 1;
        }
        *g = remap[root];
    }
    (group, next as usize)
}

/// Sum frequency-weighted residue counts for each group.
fn group_profiles(
    nodes: &[Node],
    frames: &[[u8; MOTIF_COLUMNS]],
    assignment: &[u32],
    groups: usize,
) -> Vec<Profile> {
    let mut out = vec![empty_profile(); groups];
    for (index, node) in nodes.iter().enumerate() {
        let target = &mut out[assignment[index] as usize];
        let weight = node.frequency as f64;
        for (column, &residue) in frames[index].iter().enumerate() {
            if (residue as usize) < ALPHABET {
                target[column * ALPHABET + residue as usize] += weight;
            }
        }
    }
    out
}

/// Convert counts to a smoothed probability profile.
fn normalise(counts: &Profile, prior: &Prior) -> Vec<f64> {
    let mut out = vec![0.0; MOTIF_COLUMNS * ALPHABET];
    for column in 0..MOTIF_COLUMNS {
        let base = column * ALPHABET;
        let mut total = 0.0;
        for a in 0..ALPHABET {
            let value = counts[base + a] + prior.alpha[a];
            out[base + a] = value;
            total += value;
        }
        for a in 0..ALPHABET {
            out[base + a] /= total;
        }
    }
    out
}

/// EM for a mixture of position weight matrices, seeded from `assignment`.
///
/// Positions are independent given the component, which is the same assumption
/// the merge criterion makes and the reason this is a motif model rather than a
/// homology model. Peptides missing a column simply contribute nothing to it,
/// so an 8-mer informs the other eight columns normally.
///
/// Returns the hard assignment, the profiles, the component weights, the number
/// of iterations and whether the log likelihood converged.
#[allow(clippy::type_complexity)]
fn expectation_maximization(
    nodes: &[Node],
    frames: &[[u8; MOTIF_COLUMNS]],
    assignment: &[u32],
    components: usize,
    prior: &Prior,
    max_iterations: usize,
    tolerance: f64,
) -> (Vec<u32>, Vec<Vec<f64>>, Vec<f64>, usize, bool) {
    const CHUNK: usize = 4096;
    let total_weight: f64 = nodes.iter().map(|n| n.frequency as f64).sum();
    let mut counts = group_profiles(nodes, frames, assignment, components);
    let mut mixing = vec![0.0f64; components];
    for (index, node) in nodes.iter().enumerate() {
        mixing[assignment[index] as usize] += node.frequency as f64;
    }
    for w in mixing.iter_mut() {
        *w = (*w / total_weight).max(1e-12);
    }

    let mut hard = assignment.to_vec();
    let mut previous = f64::NEG_INFINITY;
    let mut iterations = 0usize;
    let mut converged = false;

    for _ in 0..max_iterations {
        iterations += 1;
        let log_profiles: Vec<Vec<f64>> = counts
            .iter()
            .map(|c| normalise(c, prior).iter().map(|p| p.ln()).collect())
            .collect();
        let log_mixing: Vec<f64> = mixing.iter().map(|w| w.max(1e-300).ln()).collect();

        // Chunked so the accumulation order is fixed regardless of thread count.
        let chunks: Vec<(f64, Vec<Profile>, Vec<f64>, Vec<u32>)> = frames
            .par_chunks(CHUNK)
            .enumerate()
            .map(|(chunk_index, chunk)| {
                let start = chunk_index * CHUNK;
                let mut local_counts = vec![empty_profile(); components];
                let mut local_mixing = vec![0.0f64; components];
                let mut local_hard = Vec::with_capacity(chunk.len());
                let mut objective = 0.0;
                let mut responsibility = vec![0.0f64; components];
                for (offset, frame) in chunk.iter().enumerate() {
                    let weight = nodes[start + offset].frequency as f64;
                    let mut best = f64::NEG_INFINITY;
                    let mut best_component = 0u32;
                    for component in 0..components {
                        let profile = &log_profiles[component];
                        let mut value = log_mixing[component];
                        for (column, &residue) in frame.iter().enumerate() {
                            if (residue as usize) < ALPHABET {
                                value += profile[column * ALPHABET + residue as usize];
                            }
                        }
                        responsibility[component] = value;
                        if value > best {
                            best = value;
                            best_component = component as u32;
                        }
                    }
                    let mut sum = 0.0;
                    for value in responsibility.iter_mut() {
                        *value = (*value - best).exp();
                        sum += *value;
                    }
                    objective += weight * (best + sum.ln());
                    local_hard.push(best_component);
                    for component in 0..components {
                        let share = weight * responsibility[component] / sum;
                        if share == 0.0 {
                            continue;
                        }
                        local_mixing[component] += share;
                        let target = &mut local_counts[component];
                        for (column, &residue) in frame.iter().enumerate() {
                            if (residue as usize) < ALPHABET {
                                target[column * ALPHABET + residue as usize] += share;
                            }
                        }
                    }
                }
                (objective, local_counts, local_mixing, local_hard)
            })
            .collect();

        let mut objective = 0.0;
        let mut next_counts = vec![empty_profile(); components];
        let mut next_mixing = vec![0.0f64; components];
        hard.clear();
        for (chunk_objective, chunk_counts, chunk_mixing, chunk_hard) in chunks {
            objective += chunk_objective;
            for component in 0..components {
                for (target, value) in next_counts[component]
                    .iter_mut()
                    .zip(chunk_counts[component].iter())
                {
                    *target += value;
                }
                next_mixing[component] += chunk_mixing[component];
            }
            hard.extend_from_slice(&chunk_hard);
        }
        let mixing_total: f64 = next_mixing.iter().sum();
        for w in next_mixing.iter_mut() {
            *w = (*w / mixing_total).max(1e-12);
        }
        counts = next_counts;
        mixing = next_mixing;

        if previous.is_finite() && (objective - previous).abs() <= tolerance * previous.abs() {
            converged = true;
            break;
        }
        previous = objective;
    }

    let profiles = counts.iter().map(|c| normalise(c, prior)).collect();
    (hard, profiles, mixing, iterations, converged)
}

/// Run the motif stage over a finished similarity clustering.
pub fn build_motifs(
    nodes: &[Node],
    clustering: &Clustering,
    params: &MotifParams,
) -> MotifResult {
    let frames: Vec<[u8; MOTIF_COLUMNS]> = nodes
        .iter()
        .map(|node| frame(&node.sequence_codes))
        .collect();
    let background = background(nodes, &frames);
    let merge_prior = Prior::new(params.prior_concentration, &background);

    let cluster_count = clustering.representatives.len();
    let mut profiles = group_profiles(nodes, &frames, &clustering.cluster_of, cluster_count);
    let (group_of_cluster, merged_count) = agglomerate(
        &mut profiles, &merge_prior, params.merge_threshold, params.target_count);
    let merges = cluster_count.saturating_sub(merged_count);

    let mut motif_of: Vec<u32> = clustering
        .cluster_of
        .iter()
        .map(|&cluster| group_of_cluster[cluster as usize])
        .collect();

    let em_prior = Prior::new(params.em_prior_concentration, &background);
    let (final_profiles, weights, em_iterations, em_converged) =
        if params.em && merged_count > 1 {
            let (hard, profiles, mixing, iterations, converged) = expectation_maximization(
                nodes,
                &frames,
                &motif_of,
                merged_count,
                &em_prior,
                params.em_max_iterations,
                params.em_tolerance,
            );
            motif_of = hard;
            (profiles, mixing, iterations, converged)
        } else {
            let counts = group_profiles(nodes, &frames, &motif_of, merged_count.max(1));
            let mut weights = vec![0.0f64; merged_count.max(1)];
            for (index, node) in nodes.iter().enumerate() {
                weights[motif_of[index] as usize] += node.frequency as f64;
            }
            let total: f64 = weights.iter().sum();
            if total > 0.0 {
                for w in weights.iter_mut() {
                    *w /= total;
                }
            }
            (
                counts.iter().map(|c| normalise(c, &em_prior)).collect(),
                weights,
                0,
                true,
            )
        };

    let mut occupied = vec![false; merged_count.max(1)];
    for &motif in &motif_of {
        occupied[motif as usize] = true;
    }
    MotifResult {
        motif_count: occupied.iter().filter(|x| **x).count(),
        motif_of,
        merged_count,
        merges,
        em_iterations,
        em_converged,
        profiles: final_profiles,
        weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::aa_code;

    fn codes(peptide: &str) -> Vec<u8> {
        peptide.bytes().map(aa_code).collect()
    }

    #[test]
    fn ln_gamma_matches_factorials() {
        for (x, expected) in [
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 2.0f64.ln()),
            (5.0, 24.0f64.ln()),
            (10.0, 362_880.0f64.ln()),
        ] {
            assert!(
                (ln_gamma(x) - expected).abs() < 1e-10,
                "ln_gamma({x}) = {} expected {expected}",
                ln_gamma(x)
            );
        }
        // Half-integer check against ln(sqrt(pi)).
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-10);
    }

    #[test]
    fn nine_mer_frame_is_the_identity() {
        let peptide = codes("APRTAAFLL");
        assert_eq!(frame(&peptide).to_vec(), peptide);
    }

    #[test]
    fn eight_mer_gaps_the_central_column() {
        let framed = frame(&codes("APRTAAFL"));
        assert_eq!(framed[4], MISSING);
        // First four and last four residues keep their order around the gap.
        assert_eq!(framed[0], aa_code(b'A'));
        assert_eq!(framed[3], aa_code(b'T'));
        assert_eq!(framed[5], aa_code(b'A'));
        assert_eq!(framed[8], aa_code(b'L'));
    }

    #[test]
    fn long_peptide_drops_the_centre_and_keeps_both_termini() {
        // A(1)P(2)R(3)T(4) A(5)A(6) F(7)L(8)L(9)Q(10)K(11)
        let framed = frame(&codes("APRTAAFLLQK"));
        assert_eq!(framed[0], aa_code(b'A'));
        assert_eq!(framed[3], aa_code(b'T'));
        assert_eq!(framed[4], aa_code(b'F'));
        assert_eq!(framed[8], aa_code(b'K'));
        assert!(framed.iter().all(|&c| c != MISSING));
    }

    #[test]
    fn bayes_factor_prefers_one_profile_for_like_columns() {
        // One column, flat prior. Two clusters that differ only by sampling
        // noise should merge; one dominated by a different residue should not.
        let uniform = [1.0 / ALPHABET as f64; ALPHABET];
        let prior = Prior::new(1.0, &uniform);
        let leu = aa_code(b'L') as usize;
        let val = aa_code(b'V') as usize;
        let lys = aa_code(b'K') as usize;

        let mut a = empty_profile();
        a[leu] = 8.0;
        a[val] = 2.0;
        let mut b = empty_profile();
        b[leu] = 7.0;
        b[val] = 3.0;
        let mut c = empty_profile();
        c[lys] = 9.0;
        c[val] = 1.0;

        let mut ab = a;
        for (m, v) in ab.iter_mut().zip(b.iter()) {
            *m += v;
        }
        let mut ac = a;
        for (m, v) in ac.iter_mut().zip(c.iter()) {
            *m += v;
        }
        let like = prior.log_marginal(&ab) - prior.log_marginal(&a) - prior.log_marginal(&b);
        let unlike = prior.log_marginal(&ac) - prior.log_marginal(&a) - prior.log_marginal(&c);
        assert!(like > 0.0, "like-for-like columns should merge, got {like}");
        assert!(unlike < 0.0, "different columns should not merge, got {unlike}");
        assert!(like > unlike);
    }

    #[test]
    fn marginal_likelihood_is_symmetric_in_cluster_order() {
        let uniform = [1.0 / ALPHABET as f64; ALPHABET];
        let prior = Prior::new(2.5, &uniform);
        let mut a = empty_profile();
        let mut b = empty_profile();
        for column in 0..MOTIF_COLUMNS {
            a[column * ALPHABET + column % ALPHABET] = 5.0;
            b[column * ALPHABET + (column + 3) % ALPHABET] = 4.0;
        }
        let mut ab = a;
        for (m, v) in ab.iter_mut().zip(b.iter()) {
            *m += v;
        }
        let mut ba = b;
        for (m, v) in ba.iter_mut().zip(a.iter()) {
            *m += v;
        }
        assert!((prior.log_marginal(&ab) - prior.log_marginal(&ba)).abs() < 1e-9);
    }

    #[test]
    fn tri_index_is_a_bijection_on_the_upper_triangle() {
        let k = 7;
        let mut seen = vec![false; k * (k - 1) / 2];
        for i in 0..k {
            for j in (i + 1)..k {
                let index = tri(i, j, k);
                assert!(!seen[index], "duplicate index for ({i},{j})");
                seen[index] = true;
            }
        }
        assert!(seen.into_iter().all(|x| x));
    }
}
