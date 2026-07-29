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

The k-mer seed threshold defaults to 0.50. It retrieves candidate pairs only;
it never accepts a clustering relationship.

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
representatives by dynamic greedy set cover, and reuses the graph during
reassignment, representative updates, validated merging, and final validation.

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
- optional `cluster_fastas/` and graph `edges.tsv`.

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

The 0.4.3 validation configuration uses 20 independently sampled datasets of
10,000 peptides, alignment threshold 0.50, anchor threshold 0.60, terminal/core
weights 4:1, and compares graph, forced-prefilter graph, static greedy, and
lazy-exact greedy. Validation code, exact settings, logs, CSV tables, figures,
and reports are stored together under `validation/` for reproducibility.

All thresholds and gap parameters remain subject to biological calibration on
labelled peptide–MHC data. High computational agreement does not by itself
establish biological cluster purity.
