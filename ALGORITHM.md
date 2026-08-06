# PepCluster2 algorithm specification

## 1. Purpose and scope

PepCluster2 clusters peptides using information from both peptide termini and,
in the full-alignment modes, the central peptide region. Small terminal shifts
are allowed because biologically similar peptides may differ by one or more
terminal residues.

The command-line option `--mode` selects one of three scoring definitions:

| Mode | Edge-acceptance rule |
|---|---|
| `combined_kmer_anchor` | Mean of terminal k-mer similarity and anchor-combination similarity passes one threshold |
| `combined_full_anchor` | Mean of constrained full-alignment similarity and anchor-combination similarity passes one threshold |
| `separate_aln_anchor` | Constrained full-alignment similarity and anchor-combination similarity pass separate thresholds using `AND` |

`separate_aln_anchor` is the default biological model. The other modes
are retained for reproducibility, comparison, and parameter evaluation. This
document is the intended algorithm specification; a documented mode must not be
reported as available by the program until its implementation and tests exist.

Peptides shorter than eight amino acids are invalid and are excluded. All
ordering and tie-breaking rules must be independent of FASTA input order.

## 2. Terminology

PepCluster2 uses the following distinct quantities. The unqualified word
“similarity” should be avoided in reports and logs.

1. **K-mer seed score:** a thresholded dimer score used only to retrieve
   candidate peptide pairs.
2. **Terminal k-mer similarity:** an ungapped comparison of the first three and
   last three residues. It is used only by `combined_kmer_anchor`.
3. **Anchor-combination similarity:** similarity between the possible terminal
   anchor hypotheses of two peptides.
4. **Constrained full-alignment similarity:** a weighted, normalized BLOSUM62
   score for a terminally anchored full-peptide alignment.
5. **Combined score:** the arithmetic mean used only by the two `combined_*`
   compatibility modes. It is not used by `separate_aln_anchor`.
6. **Representative-ranking margin:** a ranking quantity used only after an
   edge has passed both separate thresholds. It is not an edge-acceptance
   similarity.

## 3. Terminal representation and anchor hypotheses

For a peptide of length at least eight, extract its first three and last three
residues. For terminal residues `abc...def`:

```text
front dimers = (ab, bc)
end dimers   = (de, ef)
```

The six possible ordered N-terminal/C-terminal anchor-position hypotheses are:

```text
(N1,C1), (N1,C2), (N1,C3), (N2,C2), (N2,C3), (N3,C3)
```

An anchor-position hypothesis is retained only when its two positions in the
full peptide are at least six positions apart. Each retained hypothesis is
represented by its ordered pair of residue values. Repeated values at different
valid position pairs remain separate hypotheses for similarity scoring.

Peptides with the same sequence representation and valid-position geometry may
be collapsed into one weighted node. Their input frequency is retained.

## 4. Normalized BLOSUM62 scores

Let `B(a,b)` be the BLOSUM62 score for amino acids `a` and `b`. Define normalized
residue similarity as

```text
r(a,b) = B(a,b) / sqrt(B(a,a) * B(b,b)).
```

For ordered dimers `ab` and `cd`, define

```text
d(ab,cd) = [r(a,c) + r(b,d)] / 2.
```

The 400 x 400 dimer scores are precomputed, quantized to increments of 0.001,
stored in a compact binary file, and memory-mapped read-only. Negative scores
are retained while assignments are optimized. A final reported terminal k-mer
similarity or anchor-combination similarity is clamped to `[0,1]`.

## 5. K-mer candidate retrieval

### 5.1 Seed geometry

The seed must enumerate the residue columns that the accepted alignment is
required to contain. Section 8.1 requires at least
`--minimum-terminal-match-length` residue-to-residue columns drawn from the
first three residues of *both* peptides, and the same at the C terminus. For the
default of two, those two columns may sit at any ordered position pair
`(i1 < i2)` of one peptide against any ordered pair `(j1 < j2)` of the other.
There are three such pairs inside a 3-mer, so each terminus of each peptide
contributes three ordered 2-mer keys:

```text
front keys = (N1N2, N1N3, N2N3)
end keys   = (C1C2, C1C3, C2C3)
```

`--terminal-seed all-column-pairs` (default) indexes all three.
`--terminal-seed contiguous` indexes only `(N1N2, N2N3)` and `(C1C2, C2C3)`.

The contiguous geometry is not merely less sensitive; it is misaligned with the
acceptance rule. Both contiguous dimers of a terminus contain the middle
residue, so a single substitution at that position destroys both of them while
leaving the spaced `(1,3)` pair untouched. Measured on 20 datasets of 10,000
peptides, the contiguous geometry retrieved 64% of the pairs that pass both
exact thresholds, and the failures split evenly and almost disjointly between
the two termini, which is the signature of a conjunction of two independently
insensitive tests rather than one broken terminus.

### 5.2 Seed rule

Two 2-mers are seed neighbours when

```text
d(kmer1,kmer2) >= --kmer-seed-threshold
```

The default k-mer seed threshold is 0.40. A peptide pair is retrieved as a
sensitive candidate only if it has at least one neighbouring front key pair and
at least one neighbouring end key pair. Keys of both termini are combined into
composite `front * 400 + end` index keys, so a node occupies at most nine
buckets and retrieval expands the neighbours of both components. An inverted
index and externally sorted pair records are used so the complete candidate set
need not be held in memory.

The k-mer seed score controls candidate recall and computational cost only. A
k-mer seed hit never creates a cluster edge by itself. This remains true for
all three modes.

### 5.3 Lossless anchor rejection during candidate generation

Anchor-combination similarity (Section 7) is the *mean* of the selected
one-to-one hypothesis pairs. Dropping the one-to-one constraint, so that every
hypothesis of the smaller set independently takes its best partner, can only
increase the optimum. That relaxed value is therefore an upper bound on the
exact anchor-combination similarity, and it costs at most 36 table lookups with
no dynamic program.

Candidate generation discards any pair whose relaxed bound cannot reach the
mode's acceptance rule:

```text
separate_aln_anchor:
    anchor_upper_bound >= resolved_anchor_threshold

combined_kmer_anchor, combined_full_anchor:
    (anchor_upper_bound + 1000) / 2 >= --threshold
```

The combined form uses the fact that the other component is at most 1. Because
the bound can never fall below the exact score, this rejection cannot discard an
eligible pair: it removes work, not relationships. It is applied before pairs
are written to temporary storage, so it reduces temporary disk as well as
scoring time, and it is what allows the more sensitive geometry of Section 5.1
to remain inexpensive.

Validation must report candidate cost as a decomposition — index hits, pairs
rejected by the bound, distinct pairs exactly scored, and constrained-alignment
evaluations — because these are different costs and a single figure hides which
one dominates.

## 6. Terminal k-mer similarity

This quantity is used only by `combined_kmer_anchor`. Compare the first three
residues of peptide A with the first three residues of peptide B at identical
relative positions, without gaps. Repeat for the last three residues:

```text
front_similarity(A,B)
    = mean(r(A_N1,B_N1), r(A_N2,B_N2), r(A_N3,B_N3))

end_similarity(A,B)
    = mean(r(A_C1,B_C1), r(A_C2,B_C2), r(A_C3,B_C3))

terminal_kmer_similarity(A,B)
    = [front_similarity(A,B) + end_similarity(A,B)] / 2.
```

This is the aligned-terminal 3-mer formulation. It is not the older
highest-scoring one-to-one assignment of two dimers, and it is not the average
of all four dimer cross-comparisons. Those experimental definitions are not
algorithm modes.

## 7. Anchor-combination similarity

For every valid anchor hypothesis of peptide A and every valid anchor hypothesis
of peptide B, compute normalized dimer score `d`. Find the highest-scoring
one-to-one assignment using a small bit-mask dynamic program. No hypothesis may
be reused. If the peptides contain different numbers of hypotheses, every
hypothesis in the smaller set is assigned to a distinct hypothesis in the
larger set.

```text
anchor_combination_similarity(A,B)
    = best_one_to_one_score_sum / number_of_selected_pairs.
```

The final value is clamped to `[0,1]`.

## 8. Constrained full-alignment similarity

This quantity is used by `combined_full_anchor` and `separate_aln_anchor`. It is
a single constrained alignment, not three independently concatenated
alignments.

### 8.1 Alignment construction

For each peptide, distinguish the first three residues, the core, and the last
three residues. A single affine-gap dynamic program globally optimizes the full
alignment while tracking how many residue-to-residue columns jointly originate
from the two N-terminal 3-mers and from the two C-terminal 3-mers. A valid path
must contain at least `--minimum-terminal-match-length` such columns at each
terminus; proposed default: 2.

This formulation directly optimizes the fixed-terminal idea discussed during
development: the terminal regions constrain the allowed full alignment, while
the central region is optimized conditional on satisfying both terminal
constraints. It does not independently optimize and concatenate three
alignments. Residues before the first aligned column or after the last aligned
column are retained as terminal overhangs and receive the lower terminal gap
penalties. They are never free.

Only the optimum score is required for clustering. All dynamic-programming
states, traversal order, arithmetic, normalization, and score quantization are
fixed, making the reported score deterministic and symmetric with respect to
peptide order.

### 8.2 Column weights and gaps

A residue has weight 4 when it belongs to the original first or last three
positions of its peptide and weight 1 when it belongs to the core. For a
residue-to-residue column containing residues `a` from A and `b` from B, use

```text
w_column = [w_A(a) + w_B(b)] / 2
column_score = w_column * B(a,b).
```

Gap penalties are added separately and are not multiplied by residue weights.
For a gap of length `L >= 1`:

```text
gap_cost(L) = gap_open + (L - 1) * gap_extension.
```

Proposed adjustable defaults are:

```text
--gap-open                    -4
--gap-extension               -1
--terminal-overhang-gap-open  -2
--terminal-overhang-gap-extension -1
--minimum-terminal-match-length 2
```

Thus, terminal shifts are penalized but less severely than internal insertions
or deletions. All gap parameters must remain user-adjustable for calibration on
labelled peptide–MHC data.

### 8.3 Symmetric normalization

Let `S_AB` be the weighted substitution sum plus gap penalties for the selected
constrained alignment. Let `S_AA` and `S_BB` be the weighted, gap-free self
scores of the complete peptides:

```text
S_AA = sum_i w_A(i) * B(A_i,A_i)
S_BB = sum_j w_B(j) * B(B_j,B_j).
```

Define

```text
constrained_full_alignment_similarity(A,B)
    = clamp(S_AB / sqrt(S_AA * S_BB), 0, 1).
```

This normalization is symmetric, accounts for peptide length and composition,
and prevents a terminal overhang from retaining a perfect score.

## 9. Scoring modes and threshold semantics

### 9.1 `combined_kmer_anchor`

```text
combined_score
    = [terminal_kmer_similarity + anchor_combination_similarity] / 2

eligible edge iff combined_score >= --threshold.
```

### 9.2 `combined_full_anchor`

```text
combined_score
    = [constrained_full_alignment_similarity
       + anchor_combination_similarity] / 2

eligible edge iff combined_score >= --threshold.
```

This mode measures how the former averaging rule behaves after terminal k-mer
similarity is replaced by the constrained full alignment.

### 9.3 `separate_aln_anchor`

```text
eligible edge iff
    constrained_full_alignment_similarity
        >= --alignment-similarity-threshold
AND anchor_combination_similarity
        >= --anchor-combination-similarity-threshold.
```

The two components must be reported separately. One cannot compensate for
failure of the other. The resolved defaults are alignment similarity 0.50 and
anchor-combination similarity 0.60.

For this mode, CLI precedence is:

1. With no threshold flags, alignment resolves to 0.50 and anchor-combination
   similarity resolves to 0.60.
2. An explicitly supplied `--threshold X` sets both component thresholds to
   `X`.
3. An explicitly supplied `--alignment-similarity-threshold A` overrides only
   the alignment threshold.
4. An explicitly supplied `--anchor-combination-similarity-threshold H`
   overrides only the anchor-combination threshold.

All thresholds range from 0 to 1. The resolved values, their source (default,
`--threshold`, or component override), and the selected mode must be written to
the run configuration and summary. Equal component thresholds are a convenient
starting point, not an assumption that their optimal calibrated values are
equal.

## 10. Clustering methods

Scoring and clustering are independent choices. `--mode` selects how a
candidate peptide pair is accepted. `--clustering-method` selects how accepted
relationships are used to construct representative-centred clusters:

```text
--clustering-method graph
--clustering-method greedy
```

Both methods must finish with the same mode-specific representative-to-member
invariant. They are not expected to produce identical partitions because their
initial representative-selection objectives differ.

### 10.1 Graph method

The graph method materializes every accepted candidate pair as an undirected
edge. The stored graph is reused for representative update, reassignment,
merging, and validation. Initial representatives are selected in one of two
orders, chosen by `--representative-order`.

`coverage` (default) is dynamic greedy set cover: at each selection, the method
prefers the node covering the greatest number of still unassigned nodes, then the
greatest accepted-edge weight, then frequency, then canonical sequence order.
This minimizes the cluster count.

`intrinsic` visits peptides once, in an order computed only from the peptide
itself — longer peptides first, then canonical sequence order — and lets each
unassigned peptide become a representative and absorb its still unassigned
neighbours. Input frequency is deliberately excluded from the key because it
depends on which duplicates a sample retained.

The distinction matters for dataset-size dependence. A coverage key is a
property of the whole dataset, so subsampling reorders selection and the
resulting partitions move; because the key is a degree, a small change in the
edge set is amplified into a large change in the partition. The intrinsic order
of a subset is exactly the restriction of the full-dataset order, so selection
does not churn and the only remaining difference is representatives absent from
the subset. The trade-off is compactness: the intrinsic order produces more
clusters. Neither order is uniformly better and both must be reported.

The graph method may use the optional high-confidence prefilter in Section 11.
It is appropriate when complete accepted-edge topology or dynamic set-cover
selection is desired and the graph fits the configured disk and memory limits.

### 10.2 Greedy method

The greedy method never materializes the full accepted-edge graph. It retains
the compact terminal k-mer inverted index and retrieves candidates on demand.
The initial representative rule is selected independently:

```text
--greedy-selection kmer-degree
--greedy-selection lazy-exact
```

The default `kmer-degree` pass is:

1. For every peptide, retrieve candidates satisfying at least one front and one
   end k-mer seed hit.
2. Count **distinct candidate peptides**, not raw k-mer-hit occurrences.
3. Sort possible representatives by distinct candidate count descending, then
   full peptide sequence ascending.
4. Visit possible representatives in that fixed order. Skip a peptide already
   assigned to an earlier representative.
5. For each new representative, retrieve its unassigned candidates, calculate
   the selected mode's exact score or scores, and synchronously attach every
   candidate passing the mode-specific acceptance rule.

The `lazy-exact` pass uses k-mer candidate counts only as safe
upper bounds. When a peptide reaches the top of the representative queue, its
currently unassigned candidates are scored exactly. It becomes a
representative only when its exact eligible-neighbour count and edge-weight sum
are no smaller than every remaining upper bound; otherwise it is reinserted
with its tighter exact bound. Thus it implements the graph method's dynamic
set-cover selection without storing the accepted-edge graph. Its memory use
remains bounded by the index, heap, and one retrieved candidate list, but a
candidate list may be rescored and runtime can exceed both the static greedy
and materialized-graph paths.

Candidate retrieval and scoring within one representative may run in parallel,
but representatives are committed in deterministic sorted order. Temporary
candidate identifiers are sorted and deduplicated per query. Memory therefore
scales with the peptide table, inverted index, thread-local alignment workspace,
and largest retrieved candidate list rather than the total accepted-edge count.

After the initial pass, iterative refinement performs:

1. Freeze the current representatives and build a k-mer index containing only
   those representatives.
2. For every peptide, retrieve all candidate representatives, calculate exact
   mode-specific scores, and select the best eligible representative. In
   `separate_aln_anchor`, “best” means greatest representative-ranking margin.
3. Apply all assignments simultaneously.
4. Recalculate representatives. A candidate representative must pass the
   mode-specific rule against every member; valid candidates are ranked by the
   aggregate score or aggregate representative-ranking margin.
5. Generate representative-to-representative merge candidates through the
   representative index, validate proposed unions exactly, and recalculate the
   representative after each accepted merge.
6. Validate every member against its representative and repeat until no
   assignment, representative, merge, or validation result changes, or until
   `--iteration-cap` is reached.

The greedy method has no high-confidence prefilter. `--force-prefilter` and
`--full-sensitive-after-prefilter` are invalid with
`--clustering-method greedy`; `--no-prefilter` is accepted as an explicit
statement of the greedy method's fixed behavior.

The representative-only refinement lookup can use a more permissive k-mer seed
threshold in future, analogous to PepCluster's upper-bound-guided neighbouring
block search. Such lookup remains a retrieval heuristic unless a mathematical
upper bound proves that omitted representatives cannot pass the exact scoring
rule. Validation must distinguish k-mer retrieval recall from clustering-method
agreement.

## 11. Optional high-confidence prefilter for the graph method

For `--clustering-method graph`, the prefilter is enabled by
`--force-prefilter`, disabled by `--no-prefilter`,
or selected automatically when the estimated temporary disk requirement of the
non-prefilter path exceeds the configured safe fraction of available space.

A candidate pair first must:

1. have at least one front and one end k-mer seed hit; and
2. share at least two **distinct ordered anchor-residue-pair values** among its
   valid anchor hypotheses. Repeated values such as `AL` count once.

It then must pass the high-confidence rule for the selected mode:

```text
combined_kmer_anchor:
    combined_score >= max(0.75, --threshold)

combined_full_anchor:
    combined_score >= max(0.75, --threshold)

separate_aln_anchor:
    constrained_full_alignment_similarity
        >= max(0.75, resolved_alignment_threshold)
    AND anchor_combination_similarity
        >= max(0.75, resolved_anchor_threshold).
```

Accepted pairs are **valid prefilter edges** and define provisional
high-confidence clusters. They are not final clusters. The prefilter is an
acceleration strategy and must not silently change the requested biological
edge definition.

Without prefiltering, all k-mer-seed candidates proceed directly to exact
mode-specific scoring. The k-mer inverted index is still used; “without
prefilter” does not mean exhaustive all-pairs comparison.

## 12. Graph sensitive completion and iterative clustering

After prefiltering, sensitive candidate generation includes:

- every provisional representative versus all peptides;
- every unassigned peptide versus all peptides, including assigned
  non-representatives; and
- provisional representative versus provisional representative.

A sensitive candidate becomes an **eligible sensitive edge** only when it
passes the edge rule of the selected mode at the ordinary, non-prefilter
threshold or thresholds.

The provisional partition defines candidate scope only. After sensitive edges
are built, the provisional partition is discarded and clusters are rebuilt by
deterministic greedy set cover. The following process then repeats:

```text
┌──────────────────────────────────────────────────────────────┐
│ 1. Build sequence representations and anchor hypotheses     │
└──────────────────────────────┬───────────────────────────────┘
                               ▼
┌──────────────────────────────────────────────────────────────┐
│ 2. Estimate temporary disk requirement and choose path      │
└──────────────────────┬───────────────────────┬───────────────┘
                       │ prefilter             │ no prefilter
                       ▼                       │
┌──────────────────────────────────────────┐   │
│ 3. Generate valid prefilter edges        │   │
│ 4. Build provisional clusters            │   │
└──────────────────────┬───────────────────┘   │
                       └──────────────┬────────┘
                                      ▼
┌──────────────────────────────────────────────────────────────┐
│ 5. Generate sensitive candidates and eligible edges         │
│ 6. Build initial clusters by deterministic greedy set cover │
└──────────────────────────────┬───────────────────────────────┘
                               ▼
        ┌──────────────── iterative section ─────────────────┐
        │ 7. Freeze current representatives                  │
        │ 8. Reassign every peptide synchronously            │
        │ 9. Recalculate constrained representatives         │
        │10. Propose representative-to-representative merges │
        │11. Validate proposed merged unions                 │
        │12. Recalculate representatives after merging       │
        │13. Verify every member against its representative  │
        │14. Reassign failures; cluster remaining unassigned │
        └───────────────────────┬─────────────────────────────┘
                                │
                changed ────────┴────── no change
                   ▲                         ▼
                   └──────── repeat     final clusters
```

`--iteration-cap INTEGER` optionally limits repeats of the iterative section.
It is unset by default, so iteration otherwise continues to convergence.

## 13. Representative selection, reassignment, and merging

### 13.1 Eligibility

Every representative-to-member relationship must pass the selected mode's
edge-acceptance rule. In `separate_aln_anchor`, both component thresholds must
pass independently during initial clustering, reassignment, merge validation,
and final validation.

### 13.2 Ranking eligible representatives

In a combined mode, eligible representatives are ranked first by combined
score. In `separate_aln_anchor`, define the weakest threshold margin

```text
representative_ranking_margin(A,B)
    = min(
        constrained_full_alignment_similarity(A,B)
            - resolved_alignment_threshold,
        anchor_combination_similarity(A,B)
            - resolved_anchor_threshold
      ).
```

This margin ranks pairs only after both thresholds have passed. It cannot make
an ineligible pair eligible and it is not reported as a biological similarity.
Prefer the representative with the greatest margin; break ties by greater
alignment similarity, then greater anchor-combination similarity, then greater
representative frequency, then canonical peptide sequence. This rule is also
used for synchronous reassignment.

### 13.3 Representative calculation

A representative candidate is valid only if it passes the selected mode's
edge rule against every cluster member. In combined modes, rank valid candidates
by aggregate member-to-candidate combined score. In `separate_aln_anchor`, rank
them by aggregate representative-ranking margin, followed by aggregate
alignment similarity, aggregate anchor-combination similarity, input frequency,
and canonical peptide sequence.

### 13.4 Merging and validation

A representative-to-representative eligible edge proposes a merge. A merge is
accepted only if at least one candidate representative passes the mode-specific
edge rule against every member of the proposed union. If `--merge-cap` is set,
its deterministic sample may reject a proposal early, but a successful merge
still requires exact validation of every union member.

Representatives are recalculated after merging. Every final member is then
scored exactly against its new representative. Reassignment freezes current
representatives, examines every peptide, and applies all moves simultaneously.
The process repeats until stable or until `--iteration-cap` is reached.

Reassignment is hysteretic. A peptide leaves its current representative only when
another exceeds it by more than `--reassignment-margin`, expressed in the same
units as the representative-ranking margin; the proposed default is 0.01. Without
hysteresis, any improvement moves a peptide, including an exact tie resolved by
identifier. Near-ties are precisely the composition-sensitive case: which
representative wins depends on which peptides the dataset happens to contain, so a
small change in composition reshuffles many assignments. Measured stage by stage
on a complete edge set, reassignment rather than merging is where subset stability
is lost, so the margin applies where it matters. It is a stability parameter and
must not be described as a similarity threshold: it cannot make an ineligible pair
eligible, and the representative-to-member invariant of Section 13.4 holds
regardless of its value.

The final invariant is:

> Every peptide passes the selected mode's acceptance rule against the reported
> representative of its cluster.

This representative-centred invariant does not imply that every pair of cluster
members passes the threshold or thresholds.

## 14. Determinism and output

Canonical sequence order, integer-quantized lookup scores, globally sorted
candidate pairs, synchronous reassignment, and complete tie-breaking make the
result independent of FASTA order.

Each normal run writes:

- `clusters.tsv`: one row per input peptide and its final cluster;
- `cluster_representatives.tsv`: one row per cluster representative;
- `cluster_summary.tsv`: cluster sizes and representative-level summaries;
- `run_summary.txt`: resolved mode, settings, candidate and edge counts,
  iterations, clusters, cluster-size distribution, timings, and diagnostics;
- `run_config.txt`, `command.txt`, and `run_stats.json`: reproducibility files.

Outputs for the combined modes report both component similarities and their
combined score. Outputs for `separate_aln_anchor` report constrained
full-alignment similarity, anchor-combination similarity, both resolved
thresholds, and the representative-ranking margin. Per-cluster FASTA files are
optional through `--write-cluster-fastas`.

## 15. Optional motif layer

Sections 1 to 14 define a **similarity** partition: every peptide passes the
selected mode's rule against the representative of its cluster. This section
defines a second, optional partition over the same peptides, enabled by
`--merge-motifs`. It is off by default and does not alter any quantity defined
above.

### 15.1 Why a second partition is needed

A cluster of Section 13 is a ball of a fixed radius around a representative. A
binding motif is a product of per-position residue preferences: narrow at the
anchor positions, close to flat elsewhere. The region of sequence space it
occupies is therefore strongly anisotropic, and a ball in an additive similarity
cannot cover it. Lowering the acceptance threshold widens the ball along every
axis simultaneously, crossing anchor boundaries before it spans the tolerant
positions, so no single threshold recovers a motif. One motif fragments into many
clusters as a matter of geometry, not of calibration.

### 15.2 Frame

A motif is a fixed object of nine columns, independent of peptide length. Length
is absorbed by a projection from peptide positions onto columns:

```text
L >= 9   columns 1-4 <- positions 1..4     columns 5-9 <- positions L-4..L
L == 8   columns 1-4 <- positions 1..4     column  5   <- unobserved
                                           columns 6-9 <- positions 5..8
```

For a 9-mer the projection is the identity. For longer peptides the central
residues are discarded: they bulge out of the binding groove, contact the MHC
weakly, and carry correspondingly little allele-specific information. An 8-mer
leaves the central column unobserved rather than shifting its C-terminal residues
inward, which would place the dominant C-terminal anchor in the wrong column. A
column with no residue contributes to no likelihood, so a peptide informs the
columns it does occupy and no others.

### 15.3 Profiles and marginal likelihood

Each cluster is summarised by a nine-by-twenty matrix of frequency-weighted
residue counts. Place an independent Dirichlet prior `Dir(alpha)` on each
column's amino-acid distribution, with

```text
alpha_a = --motif-prior-concentration * 20 * background_a
```

where `background` is the frequency-weighted residue composition of the dataset.
Spreading the prior over the background rather than uniformly stops a residue
that is rare overall from being treated as equally expected at every column.

With the column distribution integrated out rather than fitted, the log marginal
likelihood of the labelled residues behind a count matrix `n` is

```text
log L(n) = sum_j [ log G(A0) - log G(A0 + N_j)
                   + sum_a ( log G(n_ja + alpha_a) - log G(alpha_a) ) ]
```

with `A0 = sum_a alpha_a` and `N_j` the residues observed at column `j`. There is
no multinomial coefficient: the data are the labelled observations, not the
unordered counts. A coefficient would not cancel in Section 15.4 and would
silently change the criterion.

Positions are assumed independent given the motif. This is the standard position
weight matrix assumption, and it is what makes the stage recover motifs rather
than homologues; its cost is that two clusters with matching per-position
marginals merge even when their joint residue distributions are disjoint.

### 15.4 Merging

For clusters `A` and `B`, compare one shared profile against two separate ones:

```text
log BF(A,B) = log L(n_A + n_B) - log L(n_A) - log L(n_B).
```

Merge greedily, always taking the pair of greatest `log BF`, while that value
exceeds `--motif-merge-threshold`. Counts add exactly, so a merged profile is the
elementwise sum and no re-reading of peptides is required.

The threshold is a prior over partitions, not a similarity. Requiring
`log BF > t` is equivalent to weighting a partition of `k` clusters by
`exp(-t * k)`, a flat per-cluster penalty. Because the same constant applies to
every candidate pair, it moves where agglomeration stops without changing the
order in which pairs merge.

A size-scaled partition prior of Chinese-restaurant form was evaluated and
rejected. Its per-merge term grows like `n log 2` for two clusters of size `n`,
while the evidence against merging grows like `n` times the per-peptide
divergence between the two motifs; measured over nine columns that divergence is
frequently below `log 2` for alleles of one supertype, so the prior outvoted the
data and the agglomeration collapsed. The flat penalty above does not have this
failure because it does not scale with cluster size.

Ties in the argmax resolve to the smallest cluster index, so the merge sequence
is reproducible.

### 15.5 Refinement

Unless `--no-motif-em` is given, a mixture of position weight matrices is fitted
by expectation-maximization, seeded from the merged partition and the
corresponding mixing weights. Peptides are assigned by maximum responsibility.

Merging can only coarsen a partition: it cannot move a peptide out of a cluster
it should not have joined, so contamination present in the input clusters
propagates through Section 15.4 unchanged. Refinement is the only stage that
repairs such errors, and it is therefore the only stage that can exceed the
accuracy obtainable by merging alone.

Seeding from Section 15.4 rather than at random makes the result deterministic
and removes the restart-and-select procedure a randomly initialised mixture
requires. It is not, on the measurements available, the source of the method's
accuracy: a randomly initialised fit of the same model at the same component
count reaches materially similar agreement with allele labels. The seed buys
reproducibility and the automatic choice of component count, not accuracy.

### 15.6 Reporting and status

The motif partition **does not** satisfy the invariant of Section 13.4. Two
peptides sharing a motif need not pass the mode's rule against any common
representative. It is therefore written to `motif_clusters.tsv` and
`motif_profiles.tsv` and never merged into the similarity outputs.

All parameters of this section are provisional. The values compiled as defaults
were selected by a sweep scored on the same held-out split used to report the
result, which is not a valid selection protocol; a nested cross-validated
selection is required before they can be described as calibrated. There is no
background or outlier component, so contaminant peptides are forced into a real
motif.

## 16. Interpretation and limitations

The existing scoped sensitive stage can omit eligible edges between two
assigned non-representatives. It therefore does not guarantee that prefilter and
non-prefilter runs recover identical clusters. `--full-sensitive-after-prefilter`
is a correctness diagnostic, but reconstructing the complete sensitive graph
removes much of the disk-saving advantage. Any discrepancy must be measured and
reported as prefilter loss.

The terminal k-mer lookup remains a heuristic on the threshold, not on the
geometry. Section 5.1 enumerates every column pair the alignment could use, so a
pair is missed only when the residues in those columns score below
`--kmer-seed-threshold` at one terminus while the rest of the alignment carries
the score. The anchor rejection of Section 5.3 is lossless and cannot contribute.
Candidate-retrieval recall must therefore still be validated separately from the
scoring thresholds, and the validation must attribute every missed pair to the
conjunct that rejected it.

The greedy clustering methods never score all retrievable pairs: they only score
representative-to-unassigned pairs. Their pair-level search recall is therefore
not comparable with the graph method's, and a validation that compares the two
numbers directly is measuring different things.

Comparing a run against a reference that stops at greedy set cover confounds two
different errors, because the run's own iterative section moves the partition
away from set cover even when its candidate search is perfect. A validation must
therefore also build a reference that applies the same iterative section to the
complete edge set. Agreement with that reference is what isolates
candidate-search loss; the difference between the two references is the size of
the refinement effect.

The algorithm guarantees representative-to-member consistency, not complete
pairwise cluster consistency. All seed, alignment, anchor, gap, and prefilter
parameters require calibration against labelled peptide–MHC data. The three
modes must be compared using the same candidate pairs when their biological
scoring behaviour is evaluated.
