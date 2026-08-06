#!/usr/bin/env python3
"""Shared pieces of the motif-layer analysis, mirroring `src/motif.rs`.

Extracted so the ablation scripts do not each re-implement (or worse, exec) the
prototype. The Rust implementation is the reference; this exists to run
experiments the binary does not expose, notably random-initialised EM at a fixed
component count, which is a control rather than a feature.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.special import gammaln

AA = "ACDEFGHIKLMNPQRSTVWY"
IDX = {a: i for i, a in enumerate(AA)}
ALPHABET = 20
NP = 9
MISSING = 20


def encode(peptides) -> np.ndarray:
    """(N, 9) frame indices; 20 marks an unobserved column.

    L >= 9 : columns 1-4 <- positions 1..4, columns 5-9 <- positions L-4..L
    L == 8 : columns 1-4 <- positions 1..4, column 5 unobserved,
             columns 6-9 <- positions 5..8
    """
    X = np.full((len(peptides), NP), MISSING, dtype=np.int64)
    for i, p in enumerate(peptides):
        L = len(p)
        if L >= NP:
            src = [0, 1, 2, 3, L - 5, L - 4, L - 3, L - 2, L - 1]
            cols = [0, 1, 2, 3, 4, 5, 6, 7, 8]
        elif L == 8:
            src = [0, 1, 2, 3, 4, 5, 6, 7]
            cols = [0, 1, 2, 3, 5, 6, 7, 8]
        else:
            continue
        for c, s in zip(cols, src):
            X[i, c] = IDX.get(p[s], MISSING)
    return X


def count_matrix(X: np.ndarray, labels: np.ndarray, k: int) -> np.ndarray:
    """(k, 9, 20) residue counts; unobserved columns contribute nothing."""
    C = np.zeros((k, NP, ALPHABET))
    for j in range(NP):
        ok = X[:, j] < ALPHABET
        if ok.any():
            np.add.at(C[:, j, :], (labels[ok], X[ok, j]), 1.0)
    return C


def background(X: np.ndarray) -> np.ndarray:
    counts = count_matrix(X, np.zeros(len(X), dtype=np.int64), 1)[0].sum(0)
    total = counts.sum()
    return counts / total if total > 0 else np.full(ALPHABET, 1.0 / ALPHABET)


def make_lml(alpha: np.ndarray):
    """Log marginal likelihood of labelled residues under one shared profile."""
    a0 = alpha.sum()
    ga = gammaln(alpha).sum()

    def lml(x):
        n = x.sum(-1)
        return (gammaln(a0) - gammaln(a0 + n) + gammaln(x + alpha).sum(-1) - ga).sum(-1)
    return lml


def em_from_counts(X: np.ndarray, C: np.ndarray, mixing: np.ndarray,
                   pseudo: np.ndarray, n_iter: int = 200, tol: float = 1e-6):
    """EM over a mixture of PWMs from explicit initial counts and weights.

    Returns (hard labels, final log likelihood). Peptides missing a column
    contribute nothing to it, so an 8-mer informs the other eight normally.
    """
    N = len(X)
    k = len(C)
    onehot = [sparse.csr_matrix((np.ones(N), (np.arange(N), X[:, j])), shape=(N, 21))
              for j in range(NP)]
    C = C.copy()
    w = np.maximum(mixing.astype(float), 1e-12)
    w = w / w.sum()
    previous = -np.inf
    R = None
    objective = -np.inf
    for _ in range(n_iter):
        theta = C + pseudo
        theta = theta / theta.sum(-1, keepdims=True)
        log_theta = np.zeros((k, NP, 21))
        log_theta[:, :, :ALPHABET] = np.log(theta)
        ll = np.zeros((N, k))
        for j in range(NP):
            ll += log_theta[:, j, :][:, X[:, j]].T
        ll += np.log(np.maximum(w, 1e-300))
        top = ll.max(1, keepdims=True)
        e = np.exp(ll - top)
        s = e.sum(1, keepdims=True)
        objective = float((np.log(s) + top).sum())
        R = e / s
        w = R.sum(0)
        w = np.maximum(w, 1e-12)
        w = w / w.sum()
        C = np.zeros((k, NP, ALPHABET))
        for j in range(NP):
            C[:, j, :] = (R.T @ onehot[j])[:, :ALPHABET]
        if abs(objective - previous) <= tol * abs(objective):
            break
        previous = objective
    return R.argmax(1), objective


def random_seeded_em(X: np.ndarray, k: int, pseudo: np.ndarray, seed: int,
                     restarts: int = 10, seed_peptides: int = 20):
    """Random-init EM at fixed k, keeping the maximum-likelihood restart.

    Components are seeded from small random subsets of peptides, so each starts
    peaked and distinct. Drawing responsibilities uniformly instead would start
    every component at the background distribution, which is a strawman rather
    than a fair random baseline.
    """
    best, best_objective = None, -np.inf
    n = len(X)
    for r in range(restarts):
        rng = np.random.default_rng(seed + r)
        m = min(seed_peptides, max(2, n // (2 * k)))
        C = np.zeros((k, NP, ALPHABET))
        for component in range(k):
            index = rng.choice(n, size=m, replace=False)
            C[component] = count_matrix(X[index], np.zeros(m, dtype=np.int64), 1)[0]
        labels, objective = em_from_counts(X, C, np.full(k, 1.0 / k), pseudo)
        if objective > best_objective:
            best, best_objective = labels, objective
    return best, best_objective
