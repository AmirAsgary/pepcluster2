# Changelog

## 0.5.0

Candidate search redesigned. Recall against an exhaustive all-versus-all
reference rises from 0.64 to about 0.97 while *fewer* candidate pairs are
written to disk than in 0.4.3.

### Fixed

- **Terminal seeding missed roughly a third of eligible pairs.** The index built
  keys from contiguous dimers only (`N1N2`, `N2N3` and `C1C2`, `C2C3`), but a
  pair is accepted on a constrained alignment whose required terminal columns may
  sit at any two of the three positions in each terminal corner — including the
  `(1,3)` combination the seed never indexed. Each terminus was therefore about
  81% sensitive, and a pair had to pass both, giving 0.81² ≈ 0.64 overall. The
  seed now enumerates all three column-pairs per terminus. Peptides of differing
  length were hit hardest, which is 87% of true pairs.

- **Reads of the spilled pair file could return short.** `Read::read` is not
  obliged to fill its buffer; the pair reader assumed it did and aborted with
  "truncated pair record". This affected every run above roughly 500k peptides.
  Records are now read in a fill loop, with a regression test.

### Added

- `separate_kmer_anchor` scoring mode, with `--kmer-similarity-threshold`.
  Similarity is positionwise BLOSUM62 over the first and last three residues,
  normalised per position as `B62(x,y) / sqrt(B62(x,x)·B62(y,y))`; the middle of
  the peptide does not contribute. No alignment is performed, so it is markedly
  faster than `separate_aln_anchor`.

- An anchor upper bound applied before a candidate is spilled to disk. Relaxing
  the one-to-one anchor assignment to a per-row maximum drops injectivity and so
  can only over-estimate the exact score, making rejection provably safe. Costs
  about 36 table lookups and no dynamic programming. Search precision rises from
  1.4% to 3.4%. A test asserts every eligible pair passes the bound.

- `--terminal-seed contiguous|all-column-pairs` (default `all-column-pairs`).
  `contiguous` reproduces 0.4.3 seeding.

- `--representative-order coverage|intrinsic`. `coverage` keeps the dynamic
  greedy set cover. `intrinsic` orders on peptide properties alone, so the
  representative order of a subset is the restriction of the full-data order,
  which makes subset partitions nested-stable at some cost in agreement.

- `--reassignment-margin`, a hysteresis threshold on refinement moves.

### Changed

- `--kmer-seed-threshold` default 0.50 → **0.40**, which is where measured recall
  reaches about 0.97.

Both default changes alter clustering output relative to 0.4.3. To reproduce
0.4.3 exactly, pass `--terminal-seed contiguous --kmer-seed-threshold 0.50`.

### Validation

`validation/2026-07-30_0.5.0-dev_search_redesign/` holds the supporting runs,
including a peptide-MHC benchmark: nested cross-validation over 56,160
configurations per scoring mode, with allele-holdout outer folds. See
`START_HERE.md` there.

## 0.4.3

Previous release.
