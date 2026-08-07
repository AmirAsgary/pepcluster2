# PepCluster2

PepCluster2 is a deterministic Rust program for clustering
MHC-I peptides. It combines terminal k-mer candidate retrieval with exact
anchor-combination and constrained full-peptide alignment scores. It is local
validation software distributed as the `pepcluster2` command-line program.

The academic algorithm specification is in [ALGORITHM.md](ALGORITHM.md).

## Current recommended model

The default scoring mode is `separate_aln_anchor`. A peptide pair is eligible
only when both conditions pass:

```text
constrained full-alignment similarity >= 0.50
AND
anchor-combination similarity >= 0.60
```

The first and last three residues have alignment weight 4; core residues have
weight 1. The exact alignment is global and affine-gap, requires at least two
matched residue columns within each pair of terminal 3-mers, and gives terminal
overhangs lower penalties than internal gaps.

The k-mer seed threshold defaults to 0.40. It retrieves candidate pairs only;
it never accepts a clustering relationship.

## Candidate retrieval

The seed indexes all three ordered residue-column pairs of each terminal 3-mer,
`(1,2)`, `(1,3)` and `(2,3)`, because the accepted alignment is required to
contain at least two matched columns drawn from the first three residues of both
peptides, and those columns may be shifted when the peptides differ in length.
Retrieval requires at least one neighbouring front column pair and at least one
neighbouring end column pair.

`--terminal-seed contiguous` indexes only the contiguous dimers. Both contiguous
dimers of a terminus contain the middle residue, so one substitution there
destroys both: on 20 datasets of 10,000 peptides that geometry retrieves 64% of
the pairs passing both exact thresholds, against 97% for the default.

Candidate generation additionally discards pairs whose relaxed
anchor-combination bound cannot reach the threshold. The relaxation drops the
one-to-one constraint of the anchor assignment, so it is an upper bound on the
exact score and the rejection cannot lose an eligible pair. It is applied before
pairs reach temporary storage, so the more sensitive geometry costs little: the
default scores 0.97% of all possible pairs.

`--threshold X` explicitly sets both component thresholds to `X` in separate
mode. Component flags take precedence:

```text
--alignment-similarity-threshold FLOAT
--anchor-combination-similarity-threshold FLOAT
```

Two compatibility scoring modes remain available:

- `combined_kmer_anchor`: mean terminal-3-mer and anchor score;
- `combined_full_anchor`: mean constrained-alignment and anchor score;
- `separate_aln_anchor`: independent alignment and anchor thresholds using
  `AND`.

## Clustering paths

### Graph

`--clustering-method graph` materializes accepted edges, selects initial
representatives, and reuses the graph during reassignment, representative
updates, validated merging, and final validation.

`--reassignment-margin` (default 0.01) adds hysteresis to reassignment: a peptide
leaves its representative only when another beats it by more than that margin.
Zero moves a peptide on any improvement at all, including an exact tie broken by
identifier, which makes assignments flip whenever the dataset composition
changes. Reassignment, not merging, is the stage that costs subset stability.

`--representative-order` selects how representatives are chosen.
`coverage` (default) is dynamic greedy set cover and gives the fewest clusters.
`intrinsic` visits peptides in an order derived only from the peptide itself, so
the order of a subset is the restriction of the full-dataset order; it produces
more clusters but is markedly less sensitive to dataset composition and to small
changes in the edge set.

Graph prefiltering is optional. Automatic selection uses the estimated
temporary-disk requirement; it can be controlled with `--force-prefilter` or
`--no-prefilter`. Scoped prefilter completion is approximate and can differ
from the non-prefilter graph.

### Greedy

`--clustering-method greedy` retrieves and scores candidates on demand without
storing the accepted-edge graph. It supports two representative-selection
rules:

- `--greedy-selection kmer-degree`: fast static ordering by the number of
  distinct k-mer candidates;
- `--greedy-selection lazy-exact`: dynamic exact eligible-neighbour coverage,
  closely approximating graph set cover while retaining bounded memory.

Both greedy variants perform iterative synchronous reassignment, constrained
representative recalculation, strict cluster merging, and exact final
representative-to-member validation. Prefilter flags do not apply to greedy.

Every successful path guarantees:

> Every final member passes the selected scoring rule against its reported
> cluster representative.

This does not imply that every pair of members within a cluster passes.

## Motif layer (optional, off by default)

`--merge-motifs` adds a second partition above the similarity clusters.

A cluster produced by the sections above is a similarity ball: every member
passes the scoring rule against its representative. A binding motif is a
different object — a product of per-position residue preferences, narrow at the
anchors and near-flat elsewhere. A ball cannot cover such a region, and lowering
the threshold widens it along every axis rather than only the tolerant ones, so
one motif fragments into many clusters. On the peptide-MHC benchmark the
similarity clustering reaches BCubed recall 0.062 at ~175 clusters per pool: the
clusters are enriched for their allele but far too small.

The stage summarises each cluster as amino-acid counts on a nine-column frame,
optionally merges pairs that a Dirichlet-multinomial marginal likelihood says
came from one profile,

```text
log BF = log L(counts_A + counts_B) - log L(counts_A) - log L(counts_B)
```

then fits a mixture of position weight matrices by EM, seeded deterministically
rather than at random.

### Three variants

```bash
# 1. merge and refine (default): the motif count is chosen by the data
pepcluster2 -i peptides.fasta -o out --merge-motifs

# 2. refine only: skip the merge, finer partition, count still chosen by the data
pepcluster2 -i peptides.fasta -o out --merge-motifs --no-motif-merge

# 3. given count: return exactly K motifs
pepcluster2 -i peptides.fasta -o out --merge-motifs --motif-count 12
```

On 48 independent test pools:

| Variant | AMI | Purity | Precision | Recall | F1 | Motifs |
|---|---:|---:|---:|---:|---:|---:|
| Merge and refine | 0.596 | 0.462 | 0.512 | 0.716 | 0.568 | 7.4 |
| Refine only | 0.606 | 0.518 | 0.562 | 0.652 | 0.579 | 10.0 |
| Given count | 0.620 | 0.511 | 0.556 | 0.690 | 0.600 | exactly K |
| *similarity clustering alone* | *0.333* | *0.497* | *0.547* | *0.062* | *0.110* | *174.9* |
| MixMHCp, as documented | 0.392 | 0.222 | 0.299 | 0.828 | 0.392 | 4.1 |
| MixMHCp, given the true count | 0.492 | 0.418 | 0.471 | 0.484 | 0.473 | 12.4 |
| GibbsCluster, as documented | 0.180 | 0.100 | 0.195 | 0.605 | 0.266 | 2.9 |
| GibbsCluster, given the true count | 0.252 | 0.186 | 0.271 | 0.258 | 0.263 | 12.4 |

Variants 1 and 2 are statistically indistinguishable on AMI and F1 (p = 0.21 and
0.22 paired); they differ in granularity, variant 2 trading recall for precision.
Neither is declared correct.

Variant 3 is not an oracle in normal use: a sample's alleles are usually known
from typing. It helps mainly above ~12 alleles, where the automatic count
saturates, and costs nothing below.

### Cost

Serial runs on an exclusive node, nine pools from 995 to 11,656 peptides. Median
CPU seconds per pool:

| | CPU s | Relative |
|---|---:|---:|
| similarity clustering alone | 0.88 | 0.57× |
| merge and refine | 1.55 | 1.00× |
| given count | 0.85 | 0.55× |
| MixMHCp | 7.87 | 5.06× |
| GibbsCluster | 459.52 | 296× |

The frame is nine columns. For `L >= 9` it takes peptide positions 1–4 and
`L-4..L`, so a 9-mer maps identically and the centre of a longer peptide is
dropped — it bulges out of the binding groove and carries little allele-specific
signal. An 8-mer fills columns 1–4 and 6–9 and leaves column 5 unobserved.

Writes `motif_clusters.tsv` and `motif_profiles.tsv` alongside the usual output.

Three things to know before relying on it:

- **The motif partition does not satisfy the representative-to-member
  invariant.** Two peptides in one motif need not pass the scoring rule against
  any common representative. That is the point of the stage, and it is why the
  motif layer is written to its own files and never replaces `clusters.tsv`.
- **Refinement carries the method, not merging.** Merging alone reaches AMI
  0.430; adding EM reaches 0.596. `--motif-em-prior-concentration` moves AMI by
  0.34 across its range and is the parameter to tune first on new data.
- **A single Dirichlet prior, not a mixture.** The criterion treats the twenty
  amino acids as unordered categories, so it is blind to chemical similarity
  between residues.

Cost scales with the number of clusters, not the number of peptides: the peptides
are read once to build the profiles, after which the merge is `O(K^2)`
marginal-likelihood evaluations.

## Install

From PyPI:

```bash
pip install pepcluster2
pepcluster2 --version
```

From source:

```bash
cd /home/amir/amir/ParseFold/Pepcluster2
cargo build --release
target/release/pepcluster2 --version
```

## Examples

Non-prefilter graph with the recommended defaults:

```bash
target/release/pepcluster2 \
    --input peptides.fasta \
    --output-dir results/graph \
    --mode separate_aln_anchor \
    --clustering-method graph \
    --no-prefilter \
    --threads 0
```

Memory-saving lazy greedy:

```bash
target/release/pepcluster2 \
    --input peptides.fasta \
    --output-dir results/greedy_lazy \
    --mode separate_aln_anchor \
    --clustering-method greedy \
    --greedy-selection lazy-exact \
    --threads 0
```

Explicit non-default thresholds:

```bash
target/release/pepcluster2 \
    --input peptides.fasta \
    --output-dir results/custom \
    --alignment-similarity-threshold 0.55 \
    --anchor-combination-similarity-threshold 0.65
```

## Important options

```text
--mode combined_kmer_anchor|combined_full_anchor|separate_aln_anchor
--clustering-method graph|greedy
--greedy-selection kmer-degree|lazy-exact
--representative-order coverage|intrinsic
--reassignment-margin FLOAT
--terminal-seed all-column-pairs|contiguous
--threshold FLOAT
--alignment-similarity-threshold FLOAT
--anchor-combination-similarity-threshold FLOAT
--kmer-seed-threshold FLOAT
--gap-open FLOAT
--gap-extension FLOAT
--terminal-overhang-gap-open FLOAT
--terminal-overhang-gap-extension FLOAT
--minimum-terminal-match-length 1|2|3
--iteration-cap INT
--merge-cap INT
--no-merge
--merge-motifs
--no-motif-merge
--motif-count INT
--motif-prior-concentration FLOAT
--motif-merge-threshold FLOAT
--motif-em-prior-concentration FLOAT
--motif-em-max-iterations INT
--motif-em-tolerance FLOAT
--threads INT
--tmp-dir PATH
--compact-output
--write-cluster-fastas
```

Graph-only options:

```text
--force-prefilter
--no-prefilter
--full-sensitive-after-prefilter
--candidate-buffer-mb INT
--max-memory-gb FLOAT
--keep-tmp
--write-edges
```

`--threads 0` uses all available CPUs. Score arithmetic, candidate ordering,
representative ordering, and tie-breaking are deterministic and independent of
FASTA record order.

## Input

Input is a plain-text FASTA file. Peptides must contain at least eight canonical
amino acids. Invalid records are excluded and counted by default; `--strict`
stops at the first invalid record. Identical peptide sequences are represented
once internally with their input frequency retained.

## Outputs

Normal runs write:

- `clusters.tsv`: every input record, cluster, representative, and exact
  representative-to-member component scores;
- `cluster_representatives.tsv`: one representative per cluster;
- `cluster_summary.tsv`: cluster sizes and representative summaries;
- `anchor_clusters.tsv`: one row per unique peptide sequence;
- `run_summary.txt`: readable counts, iterations, merges, timings, and
  diagnostics;
- `run_config.txt`, `command.txt`, and `run_stats.json`: resolved settings and
  machine-readable reproducibility data;
- optional `cluster_fastas/` and graph `edges.tsv`;
- with `--merge-motifs`, `motif_clusters.tsv` (motif and similarity cluster per
  sequence, plus the framed peptide) and `motif_profiles.tsv` (fitted profile and
  mixing weight per motif).

`--compact-output` writes `node_clusters.tsv` instead of rescanning the input
FASTA and is intended for high-replicate validation.

Temporary files use `--tmp-dir`; the default is `<output-dir>/tmp`. They are
removed after successful graph runs unless `--keep-tmp` is supplied.

## Performance guidance

The graph path is usually fastest when its candidate edges fit the configured
disk and memory limits. Run `--index-only` to estimate graph expansion before a
large clustering. The scoped graph prefilter reduces storage at the cost of a
measurable possibility of missing eligible relationships.

Static greedy uses less edge storage but can create more clusters because its
representative order is based on approximate k-mer degree. Lazy-exact greedy
substantially improves graph agreement without materializing the graph, but can
rescore candidate lists and therefore be slower.

## Validation status

The validation configuration uses 20 independently sampled datasets of 10,000
peptides, alignment threshold 0.50, anchor threshold 0.60, terminal/core weights
4:1, and compares graph, forced-prefilter graph, static greedy, and lazy-exact
greedy against two exhaustive references. Validation code, exact settings, logs,
CSV tables, figures, and reports are stored together under `validation/` for
reproducibility. See [VALIDATION.md](VALIDATION.md) for the current results.

Reference choice matters when reading agreement numbers. The reference scores
every pair exactly and then applies the identical clustering procedure, so a run
differs from it only through candidate search. A reference that stopped earlier
in the procedure could not be reached by a run that completes it.

All thresholds and gap parameters remain subject to biological calibration on
labelled peptide–MHC data. High computational agreement does not by itself
establish biological cluster purity.
