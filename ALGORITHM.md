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

Two dimers are seed neighbours when

```text
d(dimer1,dimer2) >= --kmer-seed-threshold
```

The proposed default k-mer seed threshold is 0.50. A peptide pair is retrieved
as a sensitive candidate only if it has at least one neighbouring front-dimer
pair and at least one neighbouring end-dimer pair. An inverted index and
externally sorted pair records are used so the complete candidate set need not
be held in memory.

The k-mer seed score controls candidate recall and computational cost only. A
k-mer seed hit never creates a cluster edge by itself. This remains true for
all three modes.

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
edge. Initial representatives are selected by dynamic greedy set cover: at each
selection, the method prefers the node covering the greatest number of still
unassigned nodes, then the greatest accepted-edge weight, then frequency, then
canonical sequence order. The stored graph is reused for representative update,
reassignment, merging, and validation.

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

## 15. Interpretation and limitations

The existing scoped sensitive stage can omit eligible edges between two
assigned non-representatives. It therefore does not guarantee that prefilter and
non-prefilter runs recover identical clusters. `--full-sensitive-after-prefilter`
is a correctness diagnostic, but reconstructing the complete sensitive graph
removes much of the disk-saving advantage. Any discrepancy must be measured and
reported as prefilter loss.

The terminal k-mer lookup can also miss a biologically valid pair before exact
scoring if it has no qualifying front or end dimer seed. Candidate-retrieval
recall must therefore be validated separately from the scoring thresholds.

The algorithm guarantees representative-to-member consistency, not complete
pairwise cluster consistency. All seed, alignment, anchor, gap, and prefilter
parameters require calibration against labelled peptide–MHC data. The three
modes must be compared using the same candidate pairs when their biological
scoring behaviour is evaluated.
