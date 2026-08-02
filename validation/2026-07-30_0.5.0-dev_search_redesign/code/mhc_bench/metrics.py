#!/usr/bin/env python3
"""Clustering quality against MHC allele labels, robust to allele imbalance.

Alleles in this data differ in abundance by up to two orders of magnitude, so
every headline number is macro-averaged over alleles: each allele contributes
equally regardless of how many peptides it has. The micro-averaged counterpart
is reported alongside, since it answers a different question (how a randomly
chosen peptide fares) and is dominated by the largest allele.

Definitions, for peptide `i` with allele `y_i` in cluster `c_i`:

    same_i      peptides in c_i carrying allele y_i, including i
    precision_i same_i / |c_i|          "how pure is the cluster I landed in"
    recall_i    same_i / |{j: y_j=y_i}| "how much of my allele came with me"

Per-allele BCubed precision and recall average these over the allele's peptides.
Purity is precision; it is inflated by fragmentation, since singletons score 1,
so it is also reported chance-corrected against the allele's prior:

    adjusted_a = (precision_a - n_a/N) / (1 - n_a/N)

which is the objective. AMI and NMI are computed with scikit-learn and are
already chance-corrected (AMI) or normalised (NMI).
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score


def contingency(alleles: np.ndarray, clusters: np.ndarray):
    allele_ids, allele_index = np.unique(alleles, return_inverse=True)
    cluster_ids, cluster_index = np.unique(clusters, return_inverse=True)
    counts = sparse.coo_matrix(
        (np.ones(len(alleles), dtype=np.int64), (allele_index, cluster_index)),
        shape=(len(allele_ids), len(cluster_ids)),
    ).tocsr()
    return allele_ids, cluster_ids, counts


def evaluate(alleles: np.ndarray, clusters: np.ndarray) -> dict:
    n = len(alleles)
    allele_ids, cluster_ids, counts = contingency(alleles, clusters)
    allele_sizes = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    cluster_sizes = np.asarray(counts.sum(axis=0)).ravel().astype(np.float64)

    squared = counts.multiply(counts)
    # precision_a = sum_c n_ac^2 / |c|  / n_a ; recall_a = sum_c n_ac^2 / n_a^2
    per_cluster = squared.multiply(sparse.csr_matrix(1.0 / cluster_sizes))
    precision = np.asarray(per_cluster.sum(axis=1)).ravel() / allele_sizes
    recall = np.asarray(squared.sum(axis=1)).ravel() / (allele_sizes ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.where(precision + recall > 0,
                      2 * precision * recall / (precision + recall), 0.0)

    prior = allele_sizes / n
    adjusted = np.where(prior < 1.0, (precision - prior) / (1.0 - prior), 0.0)

    weights = allele_sizes / allele_sizes.sum()
    singletons = int((cluster_sizes == 1).sum())

    return {
        "peptides": int(n),
        "alleles": int(len(allele_ids)),
        "clusters": int(len(cluster_ids)),
        "singleton_clusters": singletons,
        "singleton_fraction_of_clusters": singletons / len(cluster_ids),
        "singleton_fraction_of_peptides": singletons / n,
        "largest_cluster": int(cluster_sizes.max()),
        "ami": float(adjusted_mutual_info_score(alleles, clusters)),
        "nmi": float(normalized_mutual_info_score(alleles, clusters)),
        "adjusted_purity_macro": float(adjusted.mean()),
        "adjusted_purity_micro": float((adjusted * weights).sum()),
        "bcubed_precision_macro": float(precision.mean()),
        "bcubed_precision_micro": float((precision * weights).sum()),
        "bcubed_recall_macro": float(recall.mean()),
        "bcubed_recall_micro": float((recall * weights).sum()),
        "bcubed_f1_macro": float(f1.mean()),
        "bcubed_f1_micro": float((f1 * weights).sum()),
        "min_allele_adjusted_purity": float(adjusted.min()),
    }


def objective(row: dict) -> float:
    """Selection score: the three quantities to maximise, equally weighted.

    They are on the same 0-1 scale and all chance-corrected or normalised. The
    singleton constraint is applied separately as a hard filter, not folded in
    here, so that a configuration cannot trade a violation against a better mean.
    """
    return float(np.mean([row["ami"], row["nmi"], row["adjusted_purity_macro"]]))
