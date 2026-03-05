from __future__ import annotations

from typing import Any

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.linalg import eigsh


def _ensure_csr_square(affinity: np.ndarray | csr_matrix) -> csr_matrix:
    """Validate affinity and convert it to CSR sparse matrix."""
    if issparse(affinity):
        A = affinity.tocsr().astype(np.float64, copy=False)
    else:
        arr = np.asarray(affinity, dtype=np.float64)
        if arr.ndim != 2:
            raise ValueError("affinity must be a 2D array or sparse matrix.")
        A = csr_matrix(arr)

    if A.shape[0] != A.shape[1]:
        raise ValueError("affinity matrix must be square (N x N).")
    return A


def _symmetrize_if_needed(A: csr_matrix) -> csr_matrix:
    """Ensure matrix symmetry for spectral decomposition stability."""
    asym = A - A.T
    if asym.nnz == 0:
        return A
    return ((A + A.T) * 0.5).tocsr()


def _kmeans_best_of_n(
    X: np.ndarray,
    *,
    n_clusters: int,
    n_init: int,
    max_iter: int,
    random_state: int,
) -> np.ndarray:
    """Run k-means multiple times and keep the best inertia solution."""
    best_labels: np.ndarray | None = None
    best_inertia = np.inf
    rng = np.random.default_rng(random_state)

    for _ in range(max(1, int(n_init))):
        trial_seed = int(rng.integers(0, 2**31 - 1))
        centroids, labels = kmeans2(
            data=X,
            k=n_clusters,
            iter=max(1, int(max_iter)),
            minit="++",
            seed=trial_seed,
        )

        inertia = float(np.sum((X - centroids[labels]) ** 2))
        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels.astype(np.int64, copy=False)

    if best_labels is None:
        raise RuntimeError("k-means failed to produce cluster labels.")
    return best_labels


def _randomized_top_eigvecs(
    S: csr_matrix,
    *,
    n_components: int,
    oversample: int,
    n_iter: int,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Approximate top eigenpairs of a symmetric sparse matrix.

    This randomized range-finder scales better than exact eigensolvers
    for large sparse graphs while keeping memory bounded by O(n * (k+p)).
    """
    n = S.shape[0]
    l = min(n, int(n_components) + max(0, int(oversample)))
    if l <= 0:
        raise ValueError("randomized eigensolver requires positive target dimension.")

    rng = np.random.default_rng(int(random_state))
    omega = rng.standard_normal(size=(n, l))

    Y = S @ omega
    for _ in range(max(0, int(n_iter))):
        Y = S @ (S @ Y)

    Q, _ = np.linalg.qr(Y, mode="reduced")
    B = Q.T @ (S @ Q)
    eigvals_small, eigvecs_small = np.linalg.eigh(B)
    order = np.argsort(eigvals_small)[::-1]
    keep = order[:n_components]

    eigvals = eigvals_small[keep]
    eigvecs = Q @ eigvecs_small[:, keep]
    return eigvals, eigvecs


def spectral_clustering_ncut(
    affinity: np.ndarray | csr_matrix,
    *,
    n_clusters: int,
    remove_isolated: bool = True,
    kmeans_nstart: int = 10,
    kmeans_itermax: int = 100,
    random_state: int = 42,
    enforce_nonnegative: bool = True,
    eig_solver: str = "auto",
    randomized_oversample: int = 20,
    randomized_n_iter: int = 2,
    randomized_switch_n_nodes: int = 15000,
    **extra_kwargs: Any,
) -> np.ndarray:
    """
    Spectral clustering (Ng-Jordan-Weiss style) on an affinity graph.

    Notes:
    - 26-neighborhood constraint is intentionally not applied here.
    - Isolated nodes get label -1 when `remove_isolated=True`.
    """
    A = _ensure_csr_square(affinity)
    A = _symmetrize_if_needed(A)

    if enforce_nonnegative:
        A = A.maximum(0.0)

    n_nodes = A.shape[0]
    if n_nodes == 0:
        return np.empty(0, dtype=np.int64)

    degree = np.asarray(A.sum(axis=1)).ravel()
    active_mask = degree > 0

    if not remove_isolated and not np.all(active_mask):
        raise ValueError(
            "Graph contains isolated nodes (degree=0). "
            "Set remove_isolated=True to ignore them."
        )

    active_idx = np.flatnonzero(active_mask)
    if active_idx.size == 0:
        return np.full(n_nodes, -1, dtype=np.int64)

    k = int(n_clusters)
    if k <= 0:
        raise ValueError("n_clusters must be a positive integer.")
    if active_idx.size < k:
        raise ValueError(
            f"Number of active nodes ({active_idx.size}) is smaller than n_clusters ({k})."
        )

    W = A[active_idx][:, active_idx].tocsr()
    d = np.asarray(W.sum(axis=1)).ravel()
    invsqrt = np.zeros_like(d)
    pos = d > 0
    invsqrt[pos] = 1.0 / np.sqrt(d[pos])
    Dmhalf = csr_matrix((invsqrt, (np.arange(d.size), np.arange(d.size))), shape=W.shape)
    S = (Dmhalf @ W @ Dmhalf).tocsr()
    S = _symmetrize_if_needed(S)

    if W.shape[0] == 1:
        labels = np.full(n_nodes, -1, dtype=np.int64)
        labels[active_idx[0]] = 0
        return labels

    if k >= W.shape[0]:
        raise ValueError(
            "n_clusters must be smaller than the number of active nodes for eigensolver."
        )

    solver = eig_solver.lower()
    if solver == "auto":
        solver = "randomized" if W.shape[0] >= int(randomized_switch_n_nodes) else "eigsh"

    if solver == "eigsh":
        eigvals, eigvecs = eigsh(S, k=k, which="LA")
        order = np.argsort(eigvals)[::-1]
        U = eigvecs[:, order]
    elif solver == "randomized":
        _, U = _randomized_top_eigvecs(
            S,
            n_components=k,
            oversample=randomized_oversample,
            n_iter=randomized_n_iter,
            random_state=random_state,
        )
    else:
        raise ValueError("eig_solver must be one of {'auto', 'eigsh', 'randomized'}.")

    row_norm = np.linalg.norm(U, axis=1, keepdims=True)
    row_norm[~np.isfinite(row_norm) | (row_norm == 0)] = 1.0
    Y = U / row_norm

    active_labels = _kmeans_best_of_n(
        Y,
        n_clusters=k,
        n_init=kmeans_nstart,
        max_iter=kmeans_itermax,
        random_state=random_state,
    )

    labels = np.full(n_nodes, -1, dtype=np.int64)
    labels[active_idx] = active_labels
    return labels


def resolve_spectral_hyper(cluster_hyper: dict[str, Any]) -> dict[str, Any]:
    """Normalize optional hyperparameters for spectral clustering."""
    resolved = dict(cluster_hyper)
    if "K" in resolved and "n_clusters" not in resolved:
        resolved["n_clusters"] = resolved.pop("K")

    if "seed" in resolved and "random_state" not in resolved:
        resolved["random_state"] = resolved.pop("seed")

    if "nstart" in resolved and "kmeans_nstart" not in resolved:
        resolved["kmeans_nstart"] = resolved.pop("nstart")

    if "itermax" in resolved and "kmeans_itermax" not in resolved:
        resolved["kmeans_itermax"] = resolved.pop("itermax")

    if "solver" in resolved and "eig_solver" not in resolved:
        resolved["eig_solver"] = resolved.pop("solver")

    # Metadata-only keys must not be forwarded to solver kwargs.
    resolved.pop("method", None)

    if "n_clusters" not in resolved:
        raise ValueError(
            "spectral clustering requires 'n_clusters' (or alias 'K') in cluster_hyper."
        )

    resolved.setdefault("remove_isolated", True)
    resolved.setdefault("kmeans_nstart", 10)
    resolved.setdefault("kmeans_itermax", 100)
    resolved.setdefault("random_state", 42)
    resolved.setdefault("enforce_nonnegative", True)
    resolved.setdefault("eig_solver", "auto")
    resolved.setdefault("randomized_oversample", 20)
    resolved.setdefault("randomized_n_iter", 2)
    resolved.setdefault("randomized_switch_n_nodes", 15000)
    return resolved


__all__ = ["spectral_clustering_ncut", "resolve_spectral_hyper"]
