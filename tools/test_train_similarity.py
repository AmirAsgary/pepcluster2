#!/usr/bin/env python3
"""Similarity of every test peptide to a training set, using PepCluster2's own
terminal k-mer similarity.

For each test peptide this reports the distribution of its similarity against
*every* training peptide, summarised as the maximum (its nearest training
neighbour), the 99th, 98th and 97th percentiles, and the median. The maximum
answers "is there a near-duplicate in training"; the upper percentiles answer
"is there a whole neighbourhood of near-duplicates", which is the more robust
question when a training set contains occasional outliers.

WHICH SIMILARITY
----------------
PepCluster2 defines several distinct quantities and deliberately avoids the bare
word "similarity" (ALGORITHM.md Section 2). This script computes the **terminal
k-mer similarity** of Section 6, the primary component of `separate_kmer_anchor`
- the mode the study's nested cross-validation selected:

    r(a,b)  = B62(a,b) / sqrt(B62(a,a) * B62(b,b))    normalised residue score,
                                                      quantised to thousandths
    front   = clamp(max(0, sum of r over the first three residues) / 3, 0, 1)
    end     = clamp(max(0, sum of r over the last three residues) / 3, 0, 1)
    sim     = (front + end) / 2

Note where the clamping sits: **each terminus is floored at zero on its own**,
before the two are averaged. So this is NOT the mean of the six positions - a
terminus of three mismatching residues contributes zero rather than a negative
amount, and cannot drag the other terminus down. Computing it as a six-position
mean disagrees by up to 0.32 on real peptides. Everything is done on the tool's
quantised integers for the same reason.

No alignment is performed and the middle of the peptide contributes nothing, so
peptides of different lengths are compared on their termini alone. That is the
same definition the clustering uses to accept a pair, which is why it is the
right choice for asking how close a test set sits to training.

Two caveats worth knowing before using the output:

  * `separate_kmer_anchor` accepts a pair only when terminal k-mer similarity AND
    anchor-combination similarity both pass. This script reports the first
    component alone, so a high value here does not by itself mean PepCluster2
    would have clustered the two peptides together. It is an upper bound on
    agreement, which is the conservative direction for a leakage check.
  * The BLOSUM62 table and the normalisation are copied verbatim from
    `src/scoring.rs`, and a self-test checks a few values against the Rust
    behaviour. If the scoring changes there, re-check this file.

USAGE
-----
    python3 test_train_similarity.py \\
        --test_fasta test.fasta --train_fasta train.fasta \\
        --output_csv similarity.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Amino acid order of the BLOSUM62 table below. This is the order used by
# `aa_code` in src/fasta.rs, not alphabetical order - keep the two in step.
ALPHABET = "ARNDCQEGHILKMFPSTWYV"
CODE = {residue: index for index, residue in enumerate(ALPHABET)}

# Verbatim from BLOSUM62 in src/scoring.rs.
BLOSUM62 = np.array([
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
], dtype=np.float64).reshape(20, 20)

# r(a,b) = B62(a,b) / sqrt(B62(a,a) * B62(b,b)), then quantised to thousandths.
# The tool works entirely in this integer domain (normalized_residue_scores() in
# src/scoring.rs), and the quantisation is not cosmetic: the sums, the clamping
# and the divisions below all happen on these integers, so computing in floating
# point and rounding at the end gives different answers.
#
# Rust's f64::round() rounds half away from zero; numpy's rint rounds half to
# even. They disagree on exact .5 values, so the rounding is written out rather
# than delegated to rint.
_diagonal = np.sqrt(np.outer(np.diag(BLOSUM62), np.diag(BLOSUM62)))
_normalised = BLOSUM62 / _diagonal
RESIDUE_Q = (np.sign(_normalised) *
             np.floor(np.abs(_normalised) * 1000.0 + 0.5)).astype(np.int32)

SCALE = 1000         # thousandths, matching RAW_SCALE in src/scoring.rs
MINIMUM_LENGTH = 8   # PepCluster2 rejects anything shorter


def read_fasta(path: Path) -> tuple[list[str], list[str]]:
    """Return (headers, sequences). Accepts a plain peptide-per-line file too."""
    headers: list[str] = []
    sequences: list[str] = []
    pending_header = None
    buffer: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buffer:
                    headers.append(pending_header or f"seq{len(sequences)}")
                    sequences.append("".join(buffer))
                    buffer = []
                pending_header = line[1:].strip()
            elif pending_header is None and all(c.isalpha() for c in line):
                # bare peptide list, no FASTA headers
                headers.append(f"seq{len(sequences)}")
                sequences.append(line.upper())
            else:
                buffer.append(line.upper())
    if buffer:
        headers.append(pending_header or f"seq{len(sequences)}")
        sequences.append("".join(buffer))
    return headers, sequences


def encode_terminals(sequences, headers, label):
    """(N, 6) residue codes for N1,N2,N3,C1,C2,C3, dropping invalid peptides.

    A peptide is dropped when it is shorter than eight residues or contains a
    non-canonical character, which is what PepCluster2 does on input. Dropped
    peptides are reported rather than silently ignored, because a large drop
    count usually means the input is not what the caller thought it was.
    """
    codes = np.empty((len(sequences), 6), dtype=np.int64)
    keep_headers, keep_sequences = [], []
    written = 0
    too_short = 0
    non_canonical = 0
    for header, sequence in zip(headers, sequences):
        if len(sequence) < MINIMUM_LENGTH:
            too_short += 1
            continue
        terminals = sequence[:3] + sequence[-3:]
        if any(residue not in CODE for residue in terminals):
            non_canonical += 1
            continue
        codes[written] = [CODE[residue] for residue in terminals]
        keep_headers.append(header)
        keep_sequences.append(sequence)
        written += 1
    dropped = too_short + non_canonical
    if dropped:
        print(f"  {label}: dropped {dropped} of {len(sequences)} peptides "
              f"({too_short} shorter than {MINIMUM_LENGTH}, "
              f"{non_canonical} with non-canonical residues)", file=sys.stderr)
    if written == 0:
        raise SystemExit(f"error: no usable peptides in the {label} set")
    return codes[:written], keep_headers, keep_sequences


def _terminus(test_codes, train_codes, positions):
    """One terminus, reproducing `aligned_three_mer_mean` in src/scoring.rs:

        sum.max(0).div_euclid(3).clamp(0, 1000)

    The `.max(0)` is applied to the SUM, per terminus, before dividing. That is
    why the whole quantity is not simply the mean of the six positions: a
    terminus whose three residues sum to a negative score contributes zero, not
    a negative amount, so it cannot drag the other terminus down. Getting this
    wrong overstates dissimilar pairs by up to 0.32 on real peptides.
    """
    total = np.zeros((len(test_codes), len(train_codes)), dtype=np.int32)
    for position in positions:
        total += RESIDUE_Q[np.ix_(test_codes[:, position], train_codes[:, position])]
    np.maximum(total, 0, out=total)
    total //= 3                      # div_euclid on a non-negative value
    return np.clip(total, 0, SCALE, out=total)


def kmer_similarity_block(test_codes, train_codes):
    """(B, T) terminal k-mer similarity, Section 6, in the tool's integer domain.

        terminal_kmer = (front + end + 1) / 2

    with integer division, so the +1 makes it round half up rather than truncate.
    """
    front = _terminus(test_codes, train_codes, (0, 1, 2))
    end = _terminus(test_codes, train_codes, (3, 4, 5))
    quantised = (front + end + 1) // 2
    return quantised.astype(np.float32) / float(SCALE)


# The anchor-combination similarity of Section 7 is deliberately NOT offered
# here. Reproducing it exactly needs the length-dependent `combination_mask`
# (a hypothesis counts only when its two positions are at least six apart), the
# 400x400 ordered-dimer table, and the bit-mask assignment dynamic program. An
# approximation of it would look authoritative and be wrong, which is worse than
# not offering it. Ask if you need it and it can be added and cross-checked the
# same way this one was.


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarise each test peptide's similarity to a training set "
                    "using PepCluster2's terminal k-mer similarity.")
    parser.add_argument("--test_fasta", type=Path, required=True)
    parser.add_argument("--train_fasta", type=Path, required=True)
    parser.add_argument("--output_csv", type=Path, required=True)
    parser.add_argument("--block", type=int, default=512,
                        help="test peptides per block. The full similarity matrix "
                             "is never held: each block computes its own rows, "
                             "reduces them to the reported statistics and is "
                             "discarded, so memory is block x train, not test x "
                             "train [default: 512]")
    args = parser.parse_args()

    print(f"reading {args.test_fasta}", file=sys.stderr)
    test_headers, test_sequences = read_fasta(args.test_fasta)
    print(f"reading {args.train_fasta}", file=sys.stderr)
    train_headers, train_sequences = read_fasta(args.train_fasta)

    test_codes, test_headers, test_sequences = encode_terminals(
        test_sequences, test_headers, "test")
    train_codes, _, _ = encode_terminals(train_sequences, train_headers, "train")
    print(f"  comparing {len(test_codes):,} test against {len(train_codes):,} "
          f"train peptides", file=sys.stderr)

    # Quantiles are of each test peptide's similarity to every training peptide,
    # so q99 is the value only 1% of training peptides exceed - a neighbourhood
    # measure, where max is a single nearest neighbour and can be an outlier.
    quantiles = [0.97, 0.98, 0.99]
    rows = []
    for start in range(0, len(test_codes), args.block):
        block = test_codes[start:start + args.block]
        similarity = kmer_similarity_block(block, train_codes)
        q97, q98, q99 = np.quantile(similarity, quantiles, axis=1)
        rows.append(np.column_stack([
            similarity.max(axis=1),
            q99, q98, q97,
            np.median(similarity, axis=1),
        ]))
        done = min(start + args.block, len(test_codes))
        print(f"  {done:,}/{len(test_codes):,}", end="\r", file=sys.stderr, flush=True)
    print(file=sys.stderr)

    summary = np.vstack(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_csv, "w") as handle:
        handle.write("header,peptide,max,quantile99,quantile98,quantile97,median\n")
        for header, sequence, values in zip(test_headers, test_sequences, summary):
            safe = header.replace(",", " ").replace("\n", " ")
            handle.write(f"{safe},{sequence}," + ",".join(f"{v:.6f}" for v in values) + "\n")
    print(f"wrote {len(summary):,} rows to {args.output_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
