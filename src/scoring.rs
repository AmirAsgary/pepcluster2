use crate::model::{Edge, Node};
use std::cell::RefCell;
use std::cmp::Ordering;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SimdMode {
    Auto,
    On,
    Off,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScoringMode {
    CombinedKmerAnchor,
    CombinedFullAnchor,
    SeparateAlnAnchor,
    /// Terminal k-mer similarity and anchor-combination similarity, each
    /// against its own threshold. Reads only the six terminal residues, so no
    /// alignment is computed and the core never contributes.
    SeparateKmerAnchor,
}

impl ScoringMode {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "combined_kmer_anchor" => Ok(Self::CombinedKmerAnchor),
            "combined_full_anchor" => Ok(Self::CombinedFullAnchor),
            "separate_aln_anchor" => Ok(Self::SeparateAlnAnchor),
            "separate_kmer_anchor" => Ok(Self::SeparateKmerAnchor),
            _ => Err(format!("invalid --mode: {value}")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::CombinedKmerAnchor => "combined_kmer_anchor",
            Self::CombinedFullAnchor => "combined_full_anchor",
            Self::SeparateAlnAnchor => "separate_aln_anchor",
            Self::SeparateKmerAnchor => "separate_kmer_anchor",
        }
    }

    /// Whether the mode needs the constrained full-peptide alignment. The two
    /// k-mer modes read only the terminal residues.
    pub fn uses_alignment(self) -> bool {
        !matches!(self, Self::CombinedKmerAnchor | Self::SeparateKmerAnchor)
    }

    /// Whether the mode thresholds its two components independently.
    pub fn is_separate(self) -> bool {
        matches!(self, Self::SeparateAlnAnchor | Self::SeparateKmerAnchor)
    }
}

#[rustfmt::skip]
const BLOSUM62: [i8; 400] = [
     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,
];

const COMBINATIONS: [(usize, usize); 6] = [(0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 5)];
const RAW_SCALE: i32 = 1000;
const NEGATIVE_INFINITY: i32 = i32::MIN / 8;

thread_local! {
    static ALIGNMENT_WORKSPACE: RefCell<Vec<i32>> = const { RefCell::new(Vec::new()) };
}

/// Number of constrained-alignment dynamic programs run in this process. The
/// alignment dominates scoring cost, so it is reported separately from candidate
/// volume; a relaxed counter is far cheaper than the program it counts.
static ALIGNMENT_EVALUATIONS: AtomicU64 = AtomicU64::new(0);

pub fn alignment_evaluations() -> u64 {
    ALIGNMENT_EVALUATIONS.load(AtomicOrdering::Relaxed)
}

pub fn normalized_residue_scores() -> [i32; 400] {
    let mut result = [0i32; 400];
    for a in 0..20 {
        for b in 0..20 {
            let denominator =
                ((BLOSUM62[a * 20 + a] as f64) * (BLOSUM62[b * 20 + b] as f64)).sqrt();
            result[a * 20 + b] =
                ((BLOSUM62[a * 20 + b] as f64 / denominator) * 1000.0).round() as i32;
        }
    }
    result
}

pub struct Scorer {
    residue_similarity: [i16; 400],
    pair_similarity: Vec<i16>,
    mode: ScoringMode,
    threshold_q: i32,
    alignment_threshold_q: i32,
    anchor_threshold_q: i32,
    kmer_threshold_q: i32,
    prefilter_threshold_q: i32,
    prefilter_alignment_threshold_q: i32,
    prefilter_anchor_threshold_q: i32,
    prefilter_kmer_threshold_q: i32,
    gap_open_q: i32,
    gap_extension_q: i32,
    terminal_gap_open_q: i32,
    terminal_gap_extension_q: i32,
    minimum_terminal_match_length: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct SimilarityScores {
    pub terminal_kmer: u16,
    pub anchor_combination: u16,
    pub alignment: u16,
    pub combined: u16,
    pub ranking_weight: u16,
}

impl Scorer {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mode: ScoringMode,
        threshold: f64,
        alignment_threshold: f64,
        anchor_threshold: f64,
        kmer_threshold: f64,
        gap_open: f64,
        gap_extension: f64,
        terminal_gap_open: f64,
        terminal_gap_extension: f64,
        minimum_terminal_match_length: usize,
        _simd_mode: SimdMode,
    ) -> Result<Self, String> {
        for value in [threshold, alignment_threshold, anchor_threshold, kmer_threshold] {
            if !(0.0..=1.0).contains(&value) {
                return Err("similarity thresholds must be between 0 and 1".into());
            }
        }
        if minimum_terminal_match_length == 0 || minimum_terminal_match_length > 3 {
            return Err("minimum terminal match length must be between 1 and 3".into());
        }
        let residue = normalized_residue_scores();
        let residue_similarity = residue.map(|value| value as i16);
        let mut pair_similarity = vec![0i16; 400 * 400];
        for first in 0..400 {
            let a = first / 20;
            let b = first % 20;
            for second in 0..400 {
                let c = second / 20;
                let d = second % 20;
                pair_similarity[first * 400 + second] =
                    (residue[a * 20 + c] + residue[b * 20 + d]) as i16;
            }
        }
        let q = |value: f64| (value * 1000.0).round() as i32;
        let gap_q = |value: f64| (value * 2.0 * RAW_SCALE as f64).round() as i32;
        Ok(Self {
            residue_similarity,
            pair_similarity,
            mode,
            threshold_q: q(threshold),
            alignment_threshold_q: q(alignment_threshold),
            anchor_threshold_q: q(anchor_threshold),
            kmer_threshold_q: q(kmer_threshold),
            prefilter_threshold_q: q(threshold.max(0.75)),
            prefilter_alignment_threshold_q: q(alignment_threshold.max(0.75)),
            prefilter_anchor_threshold_q: q(anchor_threshold.max(0.75)),
            prefilter_kmer_threshold_q: q(kmer_threshold.max(0.75)),
            gap_open_q: gap_q(gap_open),
            gap_extension_q: gap_q(gap_extension),
            terminal_gap_open_q: gap_q(terminal_gap_open),
            terminal_gap_extension_q: gap_q(terminal_gap_extension),
            minimum_terminal_match_length,
        })
    }

    pub fn mode(&self) -> ScoringMode {
        self.mode
    }

    pub fn simd_name(&self) -> &'static str {
        "integer_affine_dp_and_bitmask_assignment"
    }

    fn active_pairs(node: &Node) -> ([u16; 6], usize) {
        let mut result = [0u16; 6];
        let mut count = 0usize;
        for (bit, (left, right)) in COMBINATIONS.iter().copied().enumerate() {
            if node.combination_mask & (1 << bit) != 0 {
                result[count] = node.codes[left] as u16 * 20 + node.codes[right] as u16;
                count += 1;
            }
        }
        (result, count)
    }

    fn best_assignment_mean(
        &self,
        mut left: [u16; 6],
        mut n_left: usize,
        mut right: [u16; 6],
        mut n_right: usize,
    ) -> u16 {
        if n_left == 0 || n_right == 0 {
            return 0;
        }
        if n_left > n_right {
            std::mem::swap(&mut left, &mut right);
            std::mem::swap(&mut n_left, &mut n_right);
        }
        let states = 1usize << n_right;
        let mut dp = [NEGATIVE_INFINITY; 64];
        dp[0] = 0;
        for row in 0..n_left {
            let mut next = [NEGATIVE_INFINITY; 64];
            for mask in 0..states {
                if dp[mask] == NEGATIVE_INFINITY || mask.count_ones() as usize != row {
                    continue;
                }
                for column in 0..n_right {
                    if mask & (1 << column) == 0 {
                        let next_mask = mask | (1 << column);
                        let pair_score = self.pair_similarity
                            [left[row] as usize * 400 + right[column] as usize]
                            as i32;
                        next[next_mask] = next[next_mask].max(dp[mask] + pair_score);
                    }
                }
            }
            dp = next;
        }
        let best = (0..states)
            .filter(|mask| mask.count_ones() as usize == n_left)
            .map(|mask| dp[mask])
            .max()
            .unwrap_or(NEGATIVE_INFINITY);
        (best.max(0) / (2 * n_left as i32)).clamp(0, 1000) as u16
    }

    /// Upper bound on `anchor_combination_similarity` obtained by dropping the
    /// one-to-one constraint: every hypothesis of the smaller set takes its best
    /// partner independently. Any assignment the exact bit-mask program can
    /// choose is also feasible here, so this value is never smaller than the
    /// exact score. It costs at most 36 table lookups and no dynamic program.
    fn anchor_upper_bound_q(&self, a: &Node, b: &Node) -> u16 {
        let (mut left, mut n_left) = Self::active_pairs(a);
        let (mut right, mut n_right) = Self::active_pairs(b);
        if n_left == 0 || n_right == 0 {
            return 0;
        }
        if n_left > n_right {
            std::mem::swap(&mut left, &mut right);
            std::mem::swap(&mut n_left, &mut n_right);
        }
        let mut total = 0i32;
        for row in 0..n_left {
            let mut best = i32::MIN;
            for column in 0..n_right {
                let score =
                    self.pair_similarity[left[row] as usize * 400 + right[column] as usize] as i32;
                best = best.max(score);
            }
            total += best;
        }
        (total.max(0) / (2 * n_left as i32)).clamp(0, 1000) as u16
    }

    /// Sound rejection test for candidate generation. A pair failing this test
    /// cannot satisfy the mode's acceptance rule, so discarding it loses no
    /// eligible relationship.
    ///
    /// In `separate_aln_anchor` the anchor threshold must be met on its own. In
    /// the combined modes the other component is at most 1000, so the combined
    /// score is at most `(anchor_upper_bound + 1000) / 2`.
    pub fn anchor_bound_passes(&self, a: &Node, b: &Node) -> bool {
        let bound = self.anchor_upper_bound_q(a, b) as i32;
        match self.mode {
            ScoringMode::SeparateAlnAnchor | ScoringMode::SeparateKmerAnchor => {
                bound >= self.anchor_threshold_q
            }
            ScoringMode::CombinedKmerAnchor | ScoringMode::CombinedFullAnchor => {
                (bound + 1000 + 1) / 2 >= self.threshold_q
            }
        }
    }

    fn aligned_three_mer_mean(&self, a: &Node, b: &Node, offset: usize) -> u16 {
        let mut sum = 0i32;
        for position in offset..offset + 3 {
            sum += self.residue_similarity
                [a.codes[position] as usize * 20 + b.codes[position] as usize]
                as i32;
        }
        sum.max(0).div_euclid(3).clamp(0, 1000) as u16
    }

    fn anchor_combination_q(&self, a: &Node, b: &Node) -> u16 {
        let (left, n_left) = Self::active_pairs(a);
        let (right, n_right) = Self::active_pairs(b);
        self.best_assignment_mean(left, n_left, right, n_right)
    }

    #[inline]
    fn residue_weight(index: usize, length: usize) -> i32 {
        if index < 3 || index >= length - 3 {
            4
        } else {
            1
        }
    }

    #[inline]
    fn dp_index(
        i: usize,
        j: usize,
        n: usize,
        terminal_states: usize,
        n_match: usize,
        c_match: usize,
        state: usize,
    ) -> usize {
        (((i * (n + 1) + j) * terminal_states + n_match) * terminal_states + c_match) * 3 + state
    }

    fn constrained_alignment_q(&self, a: &Node, b: &Node) -> u16 {
        ALIGNMENT_EVALUATIONS.fetch_add(1, AtomicOrdering::Relaxed);
        let aa = &a.sequence;
        let bb = &b.sequence;
        let m = aa.len();
        let n = bb.len();
        let required = self.minimum_terminal_match_length;
        let terminal_states = required + 1;
        let state_count = (m + 1) * (n + 1) * terminal_states * terminal_states * 3;

        let raw = ALIGNMENT_WORKSPACE.with(|workspace| {
            let mut dp = workspace.borrow_mut();
            dp.resize(state_count, NEGATIVE_INFINITY);
            dp.fill(NEGATIVE_INFINITY);
            let origin = Self::dp_index(0, 0, n, terminal_states, 0, 0, 0);
            dp[origin] = 0;

            for i in 0..=m {
                for j in 0..=n {
                    for n_match in 0..=required {
                        for c_match in 0..=required {
                            for state in 0..3 {
                                let idx = Self::dp_index(
                                    i,
                                    j,
                                    n,
                                    terminal_states,
                                    n_match,
                                    c_match,
                                    state,
                                );
                                let current = dp[idx];
                                if current == NEGATIVE_INFINITY {
                                    continue;
                                }

                                if i < m && j < n {
                                    let next_n = if i < 3 && j < 3 {
                                        (n_match + 1).min(required)
                                    } else {
                                        n_match
                                    };
                                    let next_c = if i >= m - 3 && j >= n - 3 {
                                        (c_match + 1).min(required)
                                    } else {
                                        c_match
                                    };
                                    let wa = Self::residue_weight(i, m);
                                    let wb = Self::residue_weight(j, n);
                                    let substitution = (wa + wb)
                                        * BLOSUM62[a.sequence_codes[i] as usize * 20
                                            + b.sequence_codes[j] as usize]
                                            as i32
                                        * RAW_SCALE;
                                    let next = Self::dp_index(
                                        i + 1,
                                        j + 1,
                                        n,
                                        terminal_states,
                                        next_n,
                                        next_c,
                                        0,
                                    );
                                    dp[next] = dp[next].max(current.saturating_add(substitution));
                                }

                                if i < m && state != 2 {
                                    let terminal = j == 0 || j == n;
                                    let penalty = if state == 1 {
                                        if terminal {
                                            self.terminal_gap_extension_q
                                        } else {
                                            self.gap_extension_q
                                        }
                                    } else if terminal {
                                        self.terminal_gap_open_q
                                    } else {
                                        self.gap_open_q
                                    };
                                    let next = Self::dp_index(
                                        i + 1,
                                        j,
                                        n,
                                        terminal_states,
                                        n_match,
                                        c_match,
                                        1,
                                    );
                                    dp[next] = dp[next].max(current.saturating_add(penalty));
                                }
                                if j < n && state != 1 {
                                    let terminal = i == 0 || i == m;
                                    let penalty = if state == 2 {
                                        if terminal {
                                            self.terminal_gap_extension_q
                                        } else {
                                            self.gap_extension_q
                                        }
                                    } else if terminal {
                                        self.terminal_gap_open_q
                                    } else {
                                        self.gap_open_q
                                    };
                                    let next = Self::dp_index(
                                        i,
                                        j + 1,
                                        n,
                                        terminal_states,
                                        n_match,
                                        c_match,
                                        2,
                                    );
                                    dp[next] = dp[next].max(current.saturating_add(penalty));
                                }
                            }
                        }
                    }
                }
            }
            (0..3)
                .map(|state| {
                    dp[Self::dp_index(m, n, n, terminal_states, required, required, state)]
                })
                .max()
                .unwrap_or(NEGATIVE_INFINITY)
        });

        if raw <= 0 {
            return 0;
        }
        let self_score = |node: &Node| -> i64 {
            node.sequence
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    let code = node.sequence_codes[i] as usize;
                    2i64 * Self::residue_weight(i, node.sequence.len()) as i64
                        * BLOSUM62[code * 20 + code] as i64
                        * RAW_SCALE as i64
                })
                .sum()
        };
        let denominator = ((self_score(a) as f64) * (self_score(b) as f64)).sqrt();
        ((raw as f64 / denominator) * 1000.0)
            .round()
            .clamp(0.0, 1000.0) as u16
    }

    pub fn scores(&self, a: &Node, b: &Node) -> SimilarityScores {
        let front = self.aligned_three_mer_mean(a, b, 0);
        let end = self.aligned_three_mer_mean(a, b, 3);
        let terminal_kmer = ((front as u32 + end as u32 + 1) / 2) as u16;
        let anchor_combination = self.anchor_combination_q(a, b);
        let alignment = if self.mode.uses_alignment() {
            self.constrained_alignment_q(a, b)
        } else {
            0
        };
        let component = if self.mode.uses_alignment() {
            alignment
        } else {
            terminal_kmer
        };
        let combined = ((component as u32 + anchor_combination as u32 + 1) / 2) as u16;
        let ranking_weight = match self.mode {
            ScoringMode::CombinedKmerAnchor | ScoringMode::CombinedFullAnchor => combined,
            ScoringMode::SeparateAlnAnchor => (alignment as i32 - self.alignment_threshold_q)
                .min(anchor_combination as i32 - self.anchor_threshold_q)
                .max(0) as u16,
            ScoringMode::SeparateKmerAnchor => (terminal_kmer as i32 - self.kmer_threshold_q)
                .min(anchor_combination as i32 - self.anchor_threshold_q)
                .max(0) as u16,
        };
        SimilarityScores {
            terminal_kmer,
            anchor_combination,
            alignment,
            combined,
            ranking_weight,
        }
    }

    fn eligible_scores(&self, scores: SimilarityScores, prefilter: bool) -> bool {
        match self.mode {
            ScoringMode::CombinedKmerAnchor | ScoringMode::CombinedFullAnchor => {
                scores.combined as i32
                    >= if prefilter {
                        self.prefilter_threshold_q
                    } else {
                        self.threshold_q
                    }
            }
            ScoringMode::SeparateKmerAnchor => {
                scores.terminal_kmer as i32
                    >= if prefilter {
                        self.prefilter_kmer_threshold_q
                    } else {
                        self.kmer_threshold_q
                    }
                    && scores.anchor_combination as i32
                        >= if prefilter {
                            self.prefilter_anchor_threshold_q
                        } else {
                            self.anchor_threshold_q
                        }
            }
            ScoringMode::SeparateAlnAnchor => {
                scores.alignment as i32
                    >= if prefilter {
                        self.prefilter_alignment_threshold_q
                    } else {
                        self.alignment_threshold_q
                    }
                    && scores.anchor_combination as i32
                        >= if prefilter {
                            self.prefilter_anchor_threshold_q
                        } else {
                            self.anchor_threshold_q
                        }
            }
        }
    }

    pub fn is_eligible(&self, scores: SimilarityScores) -> bool {
        self.eligible_scores(scores, false)
    }

    /// Return complete scores for an eligible pair. Separate-threshold mode
    /// rejects on the inexpensive anchor score before running alignment.
    pub fn eligible_pair_scores(&self, a: &Node, b: &Node) -> Option<SimilarityScores> {
        if self.mode == ScoringMode::SeparateAlnAnchor {
            let anchor_combination = self.anchor_combination_q(a, b);
            if (anchor_combination as i32) < self.anchor_threshold_q {
                return None;
            }
            let alignment = self.constrained_alignment_q(a, b);
            if (alignment as i32) < self.alignment_threshold_q {
                return None;
            }
            let ranking_weight = (alignment as i32 - self.alignment_threshold_q)
                .min(anchor_combination as i32 - self.anchor_threshold_q)
                .max(0) as u16;
            return Some(SimilarityScores {
                terminal_kmer: 0,
                anchor_combination,
                alignment,
                combined: ((alignment as u32 + anchor_combination as u32 + 1) / 2) as u16,
                ranking_weight,
            });
        }
        if self.mode == ScoringMode::SeparateKmerAnchor {
            // Cheap component first. Terminal k-mer similarity is six table
            // lookups; anchor-combination similarity runs a bit-mask assignment
            // dynamic program over up to six hypotheses per peptide. Computing
            // both and testing afterwards, as the generic path below does, pays
            // for the dynamic program on every pair that the k-mer test alone
            // would have rejected. The separate-alignment branch above already
            // short-circuits this way; this mode did not.
            let front = self.aligned_three_mer_mean(a, b, 0);
            let end = self.aligned_three_mer_mean(a, b, 3);
            let terminal_kmer = ((front as u32 + end as u32 + 1) / 2) as u16;
            if (terminal_kmer as i32) < self.kmer_threshold_q {
                return None;
            }
            let anchor_combination = self.anchor_combination_q(a, b);
            if (anchor_combination as i32) < self.anchor_threshold_q {
                return None;
            }
            let combined =
                ((terminal_kmer as u32 + anchor_combination as u32 + 1) / 2) as u16;
            let ranking_weight = (terminal_kmer as i32 - self.kmer_threshold_q)
                .min(anchor_combination as i32 - self.anchor_threshold_q)
                .max(0) as u16;
            return Some(SimilarityScores {
                terminal_kmer,
                anchor_combination,
                alignment: 0,
                combined,
                ranking_weight,
            });
        }
        let scores = self.scores(a, b);
        self.eligible_scores(scores, false).then_some(scores)
    }

    pub fn score_scalar(&self, a: &Node, b: &Node) -> Option<u16> {
        self.eligible_pair_scores(a, b)
            .map(|scores| scores.ranking_weight)
    }

    pub fn compare_scores(&self, left: SimilarityScores, right: SimilarityScores) -> Ordering {
        match self.mode {
            ScoringMode::CombinedKmerAnchor | ScoringMode::CombinedFullAnchor => {
                left.combined.cmp(&right.combined)
            }
            ScoringMode::SeparateAlnAnchor => left
                .ranking_weight
                .cmp(&right.ranking_weight)
                .then(left.alignment.cmp(&right.alignment))
                .then(left.anchor_combination.cmp(&right.anchor_combination)),
            ScoringMode::SeparateKmerAnchor => left
                .ranking_weight
                .cmp(&right.ranking_weight)
                .then(left.terminal_kmer.cmp(&right.terminal_kmer))
                .then(left.anchor_combination.cmp(&right.anchor_combination)),
        }
    }

    pub fn distinct_shared_anchor_types(a: &Node, b: &Node) -> usize {
        let (mut left, n_left) = Self::active_pairs(a);
        let (mut right, n_right) = Self::active_pairs(b);
        left[..n_left].sort_unstable();
        right[..n_right].sort_unstable();
        let mut i = 0usize;
        let mut j = 0usize;
        let mut shared = 0usize;
        let mut previous = None;
        while i < n_left && j < n_right {
            if left[i] < right[j] {
                i += 1;
            } else if left[i] > right[j] {
                j += 1;
            } else {
                if previous != Some(left[i]) {
                    shared += 1;
                    previous = Some(left[i]);
                }
                let value = left[i];
                while i < n_left && left[i] == value {
                    i += 1;
                }
                while j < n_right && right[j] == value {
                    j += 1;
                }
            }
        }
        shared
    }

    pub fn distinct_anchor_types(node: &Node) -> Vec<u16> {
        let (mut pairs, count) = Self::active_pairs(node);
        pairs[..count].sort_unstable();
        let mut result = pairs[..count].to_vec();
        result.dedup();
        result
    }

    pub fn score_prefilter(&self, a: &Node, b: &Node) -> Option<u16> {
        if Self::distinct_shared_anchor_types(a, b) < 2 {
            return None;
        }
        if self.mode == ScoringMode::SeparateAlnAnchor {
            let anchor = self.anchor_combination_q(a, b);
            if (anchor as i32) < self.prefilter_anchor_threshold_q {
                return None;
            }
            let alignment = self.constrained_alignment_q(a, b);
            if (alignment as i32) < self.prefilter_alignment_threshold_q {
                return None;
            }
            return Some(
                (alignment as i32 - self.alignment_threshold_q)
                    .min(anchor as i32 - self.anchor_threshold_q)
                    .max(0) as u16,
            );
        }
        let scores = self.scores(a, b);
        self.eligible_scores(scores, true)
            .then_some(scores.ranking_weight)
    }

    pub fn filter_pairs(
        &self,
        nodes: &[Node],
        pairs: &[(u32, u32)],
        prefilter: bool,
        out: &mut Vec<Edge>,
    ) {
        crate::pair_trace::record_many(pairs);
        for &(u, v) in pairs {
            let a = &nodes[u as usize];
            let b = &nodes[v as usize];
            let weight = if prefilter {
                self.score_prefilter(a, b)
            } else {
                self.score_scalar(a, b)
            };
            if let Some(weight) = weight {
                out.push(Edge { u, v, weight });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{aa_code, combination_mask, extract_anchor};

    fn node(sequence: &[u8]) -> Node {
        let anchor = extract_anchor(sequence).unwrap();
        Node {
            sequence: sequence.to_vec(),
            sequence_codes: sequence.iter().copied().map(aa_code).collect(),
            anchor,
            codes: anchor.map(aa_code),
            combination_mask: combination_mask(sequence.len()),
            frequency: 1,
        }
    }

    fn scorer(mode: ScoringMode) -> Scorer {
        Scorer::new(
            mode,
            0.6,
            0.6,
            0.6,
            0.6,
            -4.0,
            -1.0,
            -2.0,
            -1.0,
            2,
            SimdMode::Off,
        )
        .unwrap()
    }

    #[test]
    fn identical_alignment_is_one_and_symmetric() {
        let s = scorer(ScoringMode::SeparateAlnAnchor);
        let a = node(b"ALVQRVKRF");
        let b = node(b"GALVQRVKRF");
        assert_eq!(s.scores(&a, &a).alignment, 1000);
        assert_eq!(s.scores(&a, &b).alignment, s.scores(&b, &a).alignment);
    }

    #[test]
    fn terminal_weight_prioritizes_preserved_termini() {
        let s = scorer(ScoringMode::SeparateAlnAnchor);
        let reference = node(b"ALVQRVKRF");
        let shifted = node(b"GALVQRVKRFT");
        let core_change = node(b"ALVSSSKRF");
        let reordered = node(b"VALQRVFKR");
        let shift_score = s.scores(&reference, &shifted).alignment;
        assert!(s.scores(&reference, &core_change).alignment > shift_score);
        assert!(shift_score > s.scores(&reference, &reordered).alignment);
    }

    #[test]
    fn agreed_alignment_examples_have_expected_scores() {
        let s = scorer(ScoringMode::SeparateAlnAnchor);
        let reference = node(b"ALVQRVKRF");
        let examples: [(&[u8], u16); 8] = [
            (b"ALVQRVKRF", 1000),
            (b"GALVQRVKRF", 895),
            (b"ALVQRVKRFT", 907),
            (b"GALVQRVKRFT", 808),
            (b"ALVSSSKRF", 872),
            (b"VALQRVKRF", 825),
            (b"GLVQRVKRF", 847),
            (b"VALQRVFKR", 587),
        ];
        for (sequence, expected) in examples {
            assert_eq!(s.scores(&reference, &node(sequence)).alignment, expected);
        }
    }

    #[test]
    fn separate_mode_requires_both_components() {
        let s = scorer(ScoringMode::SeparateAlnAnchor);
        let a = node(b"ALVQRVKRF");
        let b = node(b"GALVQRVKRFT");
        let scores = s.scores(&a, &b);
        assert_eq!(
            s.score_scalar(&a, &b).is_some(),
            scores.alignment >= 600 && scores.anchor_combination >= 600
        );
    }

    /// The k-mer similarity is a position-wise comparison of the two terminal
    /// 3-mers only: position i against position i, no alignment, core excluded.
    #[test]
    fn kmer_similarity_is_positionwise_over_the_two_termini() {
        let s = scorer(ScoringMode::SeparateKmerAnchor);
        // Identical termini, completely different cores -> perfect k-mer score.
        let a = node(b"AAACCCCCWWW");
        let b = node(b"AAAKKKKKWWW");
        assert_eq!(s.scores(&a, &b).terminal_kmer, 1000);
        // The core is never read, so lengthening it changes nothing.
        let c = node(b"AAAKKKKKKKKKKWWW");
        assert_eq!(s.scores(&a, &c).terminal_kmer, 1000);
        // No alignment is computed in this mode.
        assert_eq!(s.scores(&a, &b).alignment, 0);
        // Positions are not realigned: shifting the N-terminus by one breaks it.
        let shifted = node(b"GAAACCCCCWWW");
        assert!(s.scores(&a, &shifted).terminal_kmer < 1000);
    }

    /// Each terminus is the mean of three residue similarities, and the two
    /// termini are averaged with equal weight.
    #[test]
    fn kmer_similarity_matches_the_documented_formula() {
        let s = scorer(ScoringMode::SeparateKmerAnchor);
        let residue = normalized_residue_scores();
        let a = node(b"AWDKKKKLMN");
        let b = node(b"APDKKKKLMN");
        let code = |c: u8| aa_code(c) as usize;
        let front: i32 = [(b'A', b'A'), (b'W', b'P'), (b'D', b'D')]
            .iter()
            .map(|(x, y)| residue[code(*x) * 20 + code(*y)])
            .sum();
        let end: i32 = [(b'L', b'L'), (b'M', b'M'), (b'N', b'N')]
            .iter()
            .map(|(x, y)| residue[code(*x) * 20 + code(*y)])
            .sum();
        let expected_front = front.max(0) / 3;
        let expected_end = end.max(0) / 3;
        let expected = (expected_front + expected_end + 1) / 2;
        assert_eq!(s.scores(&a, &b).terminal_kmer as i32, expected);
    }

    /// Both components must pass independently, and neither can compensate for
    /// the other; a zero threshold disables its component.
    #[test]
    fn separate_kmer_mode_thresholds_both_components() {
        let nodes = sample_nodes(24);
        let strict = Scorer::new(ScoringMode::SeparateKmerAnchor, 0.60, 0.50, 0.60, 0.60,
                                 -4.0, -1.0, -2.0, -1.0, 2, SimdMode::Auto).unwrap();
        let kmer_only = Scorer::new(ScoringMode::SeparateKmerAnchor, 0.60, 0.50, 0.0, 0.60,
                                    -4.0, -1.0, -2.0, -1.0, 2, SimdMode::Auto).unwrap();
        let mut both = 0usize;
        let mut relaxed = 0usize;
        for (index, a) in nodes.iter().enumerate() {
            for b in &nodes[index + 1..] {
                let scores = strict.scores(a, b);
                let expected = scores.terminal_kmer >= 600 && scores.anchor_combination >= 600;
                assert_eq!(strict.eligible_pair_scores(a, b).is_some(), expected);
                both += usize::from(expected);
                // Anchor threshold 0 leaves only the k-mer condition.
                let only = kmer_only.eligible_pair_scores(a, b).is_some();
                assert_eq!(only, scores.terminal_kmer >= 600);
                relaxed += usize::from(only);
            }
        }
        assert!(both > 0, "no eligible pair to test");
        assert!(relaxed >= both, "dropping the anchor threshold must not remove edges");
    }

    /// With neither component threshold supplied, --threshold governs both.
    #[test]
    fn threshold_precedence_for_the_kmer_mode() {
        let a = node(b"AWDKKKKLMN");
        let b = node(b"APDKKKKLMN");
        let shared = Scorer::new(ScoringMode::SeparateKmerAnchor, 0.30, 0.30, 0.30, 0.30,
                                 -4.0, -1.0, -2.0, -1.0, 2, SimdMode::Auto).unwrap();
        let split = Scorer::new(ScoringMode::SeparateKmerAnchor, 0.99, 0.50, 0.30, 0.30,
                                -4.0, -1.0, -2.0, -1.0, 2, SimdMode::Auto).unwrap();
        assert_eq!(shared.eligible_pair_scores(&a, &b).is_some(),
                   split.eligible_pair_scores(&a, &b).is_some());
    }

    #[test]
    fn repeated_anchor_values_count_once_for_prefilter() {
        let a = node(b"AAAAAALLL");
        assert_eq!(Scorer::distinct_shared_anchor_types(&a, &a), 1);
    }

    /// Deterministic peptide families: an unrelated base sequence plus close
    /// relatives produced by point substitutions and one-residue terminal
    /// shifts. Lengths span 8..=14 so both anchor-hypothesis counts (3 for
    /// length 8, 6 otherwise) are covered, and related members guarantee that
    /// eligible pairs actually occur at the default thresholds.
    fn sample_nodes(families: usize) -> Vec<Node> {
        const ALPHABET: &[u8] = b"ARNDCQEGHILKMFPSTWYV";
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) as usize
        };
        let mut nodes = Vec::new();
        for _ in 0..families {
            let length = 9 + next() % 6;
            let base: Vec<u8> = (0..length)
                .map(|_| ALPHABET[next() % ALPHABET.len()])
                .collect();
            nodes.push(node(&base));
            // Point substitutions, including inside a terminal 3-mer.
            for _ in 0..3 {
                let mut variant = base.clone();
                let position = next() % variant.len();
                variant[position] = ALPHABET[next() % ALPHABET.len()];
                nodes.push(node(&variant));
            }
            // One-residue terminal extension and trim: the shifted-terminal case.
            let mut extended = vec![ALPHABET[next() % ALPHABET.len()]];
            extended.extend_from_slice(&base);
            nodes.push(node(&extended));
            if base.len() > 8 {
                nodes.push(node(&base[1..]));
            }
        }
        nodes
    }

    #[test]
    fn anchor_upper_bound_never_underestimates_the_exact_assignment() {
        let s = scorer(ScoringMode::SeparateAlnAnchor);
        let nodes = sample_nodes(40);
        for a in &nodes {
            for b in &nodes {
                let bound = s.anchor_upper_bound_q(a, b);
                let exact = s.anchor_combination_q(a, b);
                assert!(
                    bound >= exact,
                    "relaxed bound {bound} below exact anchor score {exact} for {} vs {}",
                    String::from_utf8_lossy(&a.sequence),
                    String::from_utf8_lossy(&b.sequence)
                );
            }
        }
    }

    /// The candidate-generation gate must never discard a pair the mode would
    /// accept; this is what makes the pruning lossless rather than heuristic.
    #[test]
    fn every_eligible_pair_passes_the_candidate_bound() {
        let nodes = sample_nodes(24);
        for mode in [
            ScoringMode::SeparateAlnAnchor,
            ScoringMode::CombinedFullAnchor,
            ScoringMode::CombinedKmerAnchor,
        ] {
            let s = scorer(mode);
            let mut eligible = 0usize;
            for (index, a) in nodes.iter().enumerate() {
                for b in &nodes[index + 1..] {
                    if s.eligible_pair_scores(a, b).is_some() {
                        eligible += 1;
                        assert!(
                            s.anchor_bound_passes(a, b),
                            "mode {} discarded eligible pair {} vs {}",
                            mode.name(),
                            String::from_utf8_lossy(&a.sequence),
                            String::from_utf8_lossy(&b.sequence)
                        );
                    }
                }
            }
            assert!(eligible > 0, "mode {} produced no eligible pair", mode.name());
        }
    }
}
