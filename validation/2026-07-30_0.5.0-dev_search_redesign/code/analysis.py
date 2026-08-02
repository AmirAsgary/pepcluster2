#!/usr/bin/env python3
"""Shared analysis primitives for the PepCluster2 search-redesign validation.

Three things live here:

* partition comparison (ARI, NMI, pairwise Jaccard, co-association recall and
  precision), identical in definition to the 0.4.3 validation so numbers stay
  comparable;
* search-rule metrics against the exhaustive true-pair set, reported together
  with the cost decomposition the tool now emits;
* the missed-pair audit, which attributes every eligible pair the search failed
  to retrieve to the conjunct that rejected it.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import struct
from collections import Counter
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# scoring tables, mirrored from src/scoring.rs so the audit uses identical
# quantisation to the Rust implementation
# --------------------------------------------------------------------------

AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"
CODE = {residue: index for index, residue in enumerate(AMINO_ACIDS)}
_ANCHOR_COMBINATIONS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


def _blosum62(source: Path) -> np.ndarray:
    text = source.read_text()
    match = re.search(r"const BLOSUM62: \[i8; 400\] = \[(.*?)\];", text, re.S)
    if match is None:
        raise ValueError(f"BLOSUM62 table not found in {source}")
    values = [int(x) for x in re.findall(r"-?\d+", match.group(1))]
    if len(values) != 400:
        raise ValueError(f"expected 400 BLOSUM62 entries, found {len(values)}")
    return np.asarray(values, dtype=np.float64).reshape(20, 20)


class ScoreTables:
    """Quantised residue and ordered-2-mer similarity, as the tool computes them."""

    def __init__(self, scoring_source: Path):
        blosum = _blosum62(scoring_source)
        diagonal = np.sqrt(np.outer(np.diag(blosum), np.diag(blosum)))
        self.residue = np.rint(blosum / diagonal * 1000.0).astype(np.int32)
        summed = (
            self.residue[:, None, :, None] + self.residue[None, :, None, :]
        ).reshape(400, 400)
        # kmer.rs rounds half away from zero.
        self.dimer = np.where(
            summed >= 0, (summed + 1) // 2, -((-summed + 1) // 2)
        ).astype(np.int32)
        self.pair_sum = summed.astype(np.int32)


# --------------------------------------------------------------------------
# input/output helpers
# --------------------------------------------------------------------------


def read_fasta_sequences(path: Path) -> list[str]:
    sequences: list[str] = []
    current: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
            else:
                current.append(line)
    if current:
        sequences.append("".join(current))
    return sequences


def node_sequences(path: Path) -> list[str]:
    """Node identifiers are indices into the sorted unique sequence list, which
    is how `fasta::load_nodes` orders nodes."""
    return sorted(set(read_fasta_sequences(path)))


def read_pair_file(path: Path, magic: bytes) -> np.ndarray:
    payload = path.read_bytes()
    if payload[:8] != magic:
        raise ValueError(f"invalid pair file {path}: magic {payload[:8]!r}")
    count = struct.unpack_from("<Q", payload, 8)[0]
    return np.frombuffer(payload, dtype="<u8", count=count, offset=16)


def unpack_pairs(packed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (
        (packed >> np.uint64(32)).astype(np.int64),
        (packed & np.uint64(0xFFFFFFFF)).astype(np.int64),
    )


def read_partition(path: Path, sequence_column: str = "sequence") -> dict[str, str]:
    with path.open(newline="") as handle:
        return {
            row[sequence_column]: row["cluster_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_stats(path: Path) -> dict:
    return json.loads(path.read_text())


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------
# partition comparison
# --------------------------------------------------------------------------


def _choose2(value: int) -> int:
    return value * (value - 1) // 2


def partition_metrics(
    reference: dict[str, str], query: dict[str, str]
) -> dict[str, float]:
    keys = sorted(set(reference) & set(query))
    reference_counts = Counter(reference[key] for key in keys)
    query_counts = Counter(query[key] for key in keys)
    cells = Counter((reference[key], query[key]) for key in keys)
    reference_pairs = sum(_choose2(value) for value in reference_counts.values())
    query_pairs = sum(_choose2(value) for value in query_counts.values())
    shared = sum(_choose2(value) for value in cells.values())
    total_pairs = _choose2(len(keys))
    expected = reference_pairs * query_pairs / total_pairs if total_pairs else 0.0
    denominator = 0.5 * (reference_pairs + query_pairs) - expected
    ari = (shared - expected) / denominator if denominator else 1.0
    n = len(keys)
    mutual_information = sum(
        count / n * math.log(n * count / (reference_counts[ref] * query_counts[qry]))
        for (ref, qry), count in cells.items()
    )
    reference_entropy = -sum(
        count / n * math.log(count / n) for count in reference_counts.values()
    )
    query_entropy = -sum(
        count / n * math.log(count / n) for count in query_counts.values()
    )
    entropy_sum = reference_entropy + query_entropy
    nmi = 2 * mutual_information / entropy_sum if entropy_sum else 1.0
    union = reference_pairs + query_pairs - shared
    return {
        "ari": ari,
        "nmi": nmi,
        "pairwise_jaccard": shared / union if union else 1.0,
        "coassociation_recall": shared / reference_pairs if reference_pairs else 1.0,
        "coassociation_precision": shared / query_pairs if query_pairs else 1.0,
        "reference_cocluster_pairs": reference_pairs,
        "query_cocluster_pairs": query_pairs,
        "shared_cocluster_pairs": shared,
        "compared_peptides": n,
    }


# --------------------------------------------------------------------------
# search-rule metrics
# --------------------------------------------------------------------------


def search_metrics(
    true_pairs: set[int], scored_pairs: set[int], all_possible_pairs: int
) -> dict[str, float]:
    found = len(true_pairs & scored_pairs)
    return {
        "true_pairs": len(true_pairs),
        "scored_unique_pairs": len(scored_pairs),
        "true_pairs_found": found,
        "missed_true_pairs": len(true_pairs) - found,
        "search_recall": found / len(true_pairs) if true_pairs else 1.0,
        "search_precision": found / len(scored_pairs) if scored_pairs else 1.0,
        "fraction_all_pairs_scored": len(scored_pairs) / all_possible_pairs,
    }


def cost_metrics(stats: dict) -> dict[str, float]:
    """Cost decomposition emitted by the tool. `index_candidate_hits` counts
    index hits with multiplicity; the other stages count distinct pairs."""
    return {
        "index_candidate_hits": stats.get("index_candidate_hits", 0),
        "anchor_bound_rejected": stats.get("anchor_bound_rejected", 0),
        "candidate_pairs_computed": stats.get("candidate_pairs_computed", 0),
        "alignment_evaluations": stats.get("alignment_evaluations", 0),
        "graph_edge_count": stats.get("graph_edge_count", 0),
        "elapsed_seconds": stats.get("elapsed_seconds", float("nan")),
    }


def peak_rss_kbytes(resource_file: Path) -> float:
    if not resource_file.exists():
        return float("nan")
    for line in resource_file.read_text().splitlines():
        if "Maximum resident set size" in line:
            return float(line.rsplit(":", 1)[1])
    return float("nan")


# --------------------------------------------------------------------------
# missed-pair audit
# --------------------------------------------------------------------------


class PeptideFeatures:
    """Terminal column-pair keys and anchor hypotheses for every node."""

    def __init__(self, sequences: list[str], tables: ScoreTables):
        self.tables = tables
        self.lengths = np.asarray([len(s) for s in sequences])
        self.anchors = np.asarray(
            [
                [
                    CODE[s[0]],
                    CODE[s[1]],
                    CODE[s[2]],
                    CODE[s[-3]],
                    CODE[s[-2]],
                    CODE[s[-1]],
                ]
                for s in sequences
            ],
            dtype=np.int32,
        )
        self.hypotheses = np.stack(
            [
                self.anchors[:, front] * 20 + self.anchors[:, 3 + back]
                for front, back in _ANCHOR_COMBINATIONS
            ],
            axis=1,
        )
        self.hypothesis_valid = np.stack(
            [
                (self.lengths - 3 + back) >= (front + 6)
                for front, back in _ANCHOR_COMBINATIONS
            ],
            axis=1,
        )

    def terminal_keys(self, ids: np.ndarray, terminus: int, geometry: str) -> np.ndarray:
        offset = 0 if terminus == 0 else 3
        pairs = (
            [(0, 1), (1, 2)]
            if geometry == "contiguous"
            else [(0, 1), (0, 2), (1, 2)]
        )
        return np.stack(
            [
                self.anchors[ids, offset + i] * 20 + self.anchors[ids, offset + j]
                for i, j in pairs
            ],
            axis=1,
        )

    def terminal_seed_hit(
        self, left: np.ndarray, right: np.ndarray, terminus: int, geometry: str,
        threshold_q: int,
    ) -> np.ndarray:
        keys_left = self.terminal_keys(left, terminus, geometry)
        keys_right = self.terminal_keys(right, terminus, geometry)
        hit = np.zeros(len(left), dtype=bool)
        for i in range(keys_left.shape[1]):
            for j in range(keys_right.shape[1]):
                hit |= self.tables.dimer[keys_left[:, i], keys_right[:, j]] >= threshold_q
        return hit

    def anchor_upper_bound(self, left: np.ndarray, right: np.ndarray) -> np.ndarray:
        """Row-max relaxation of the one-to-one hypothesis assignment, matching
        `Scorer::anchor_upper_bound_q`."""
        counts_left = self.hypothesis_valid[left].sum(1)
        counts_right = self.hypothesis_valid[right].sum(1)
        smaller_is_left = counts_left <= counts_right
        rows = np.where(smaller_is_left, counts_left, counts_right)
        best = np.zeros(len(left), dtype=np.int64)
        low, high = np.where(smaller_is_left, left, right), np.where(
            smaller_is_left, right, left
        )
        total = np.zeros(len(left), dtype=np.int64)
        for i in range(6):
            active_row = self.hypothesis_valid[low, i]
            if not active_row.any():
                continue
            row_best = np.full(len(left), -(10**9), dtype=np.int64)
            for j in range(6):
                active_col = self.hypothesis_valid[high, j]
                score = np.where(
                    active_col,
                    self.tables.pair_sum[self.hypotheses[low, i], self.hypotheses[high, j]],
                    -(10**9),
                )
                row_best = np.maximum(row_best, score)
            total += np.where(active_row, row_best, 0)
        best = np.clip(np.maximum(total, 0) // (2 * np.maximum(rows, 1)), 0, 1000)
        return best


def missed_pair_audit(
    missed_packed: np.ndarray,
    features: PeptideFeatures,
    geometry: str,
    seed_threshold_q: int,
    anchor_threshold_q: int,
) -> dict[str, float]:
    """Attribute each missed eligible pair to the conjunct that rejected it.

    A pair can only be missed because a terminal seed found no neighbouring
    column pair, or because the anchor upper bound rejected it. The latter must
    never happen: the bound is a mathematical upper bound on the exact anchor
    score, so a non-zero count means the implementation is wrong.
    """
    if len(missed_packed) == 0:
        return {
            "missed": 0,
            "missed_front_seed": 0,
            "missed_end_seed": 0,
            "missed_both_seeds": 0,
            "missed_seed_attributable": 0,
            "missed_retrievable_but_unexamined": 0,
            "missed_anchor_bound_unsound": 0,
        }
    left, right = unpack_pairs(missed_packed)
    front = features.terminal_seed_hit(left, right, 0, geometry, seed_threshold_q)
    end = features.terminal_seed_hit(left, right, 1, geometry, seed_threshold_q)
    bound = features.anchor_upper_bound(left, right)
    unsound = bound < anchor_threshold_q
    seed_attributable = (~front) | (~end)
    return {
        "missed": int(len(missed_packed)),
        "missed_front_seed": int((~front).sum()),
        "missed_end_seed": int((~end).sum()),
        "missed_both_seeds": int(((~front) & (~end)).sum()),
        "missed_seed_attributable": int(seed_attributable.sum()),
        # Pairs the index would have returned and the bound would have kept, but
        # which the clustering traversal never scored. The greedy paths only ever
        # score representative-to-unassigned pairs, so this is structural for
        # them rather than a candidate-search defect.
        "missed_retrievable_but_unexamined": int((~seed_attributable & ~unsound).sum()),
        "missed_anchor_bound_unsound": int(unsound.sum()),
    }
