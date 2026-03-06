from __future__ import annotations

from typing import Any, Sequence

import numpy as np
from scipy.sparse import csr_matrix, issparse

from data.clustering import ClusteringResult


def _labels_2d(label: np.ndarray) -> np.ndarray:
    arr = np.asarray(label)
    if arr.ndim == 1:
        return arr.reshape(1, -1)
    if arr.ndim == 2:
        return arr
    raise ValueError("label must be 1D or 2D array.")


def _normalize_similarity_matrix(similarity: np.ndarray | csr_matrix) -> np.ndarray | csr_matrix:
    """Force symmetry and zero diagonal, matching ROI_detection.R strict logic."""
    if issparse(similarity):
        S = similarity.tocsr().astype(np.float64, copy=False)
        S = ((S + S.T) * 0.5).tocsr()
        S.setdiag(0.0)
        S.eliminate_zeros()
        return S
    S = np.asarray(similarity, dtype=np.float64)
    if S.ndim != 2:
        raise ValueError("similarity matrix must be 2D.")
    S = 0.5 * (S + S.T)
    np.fill_diagonal(S, 0.0)
    return S


def _silhouette_cluster_strict_similarity(
    similarity: np.ndarray | csr_matrix,
    labels: np.ndarray,
    *,
    unassigned_labels: Sequence[int] = (-1, 0),
) -> dict[str, Any]:
    """
    Port of ROI_detection.R silhouette_cluster_strict.

    Uses all pairs in the similarity matrix (not edge-only), with:
    a_k = z_k^T S z_k / (n_k * (n_k - 1))
    b_k = (z_k^T S 1 - z_k^T S z_k) / (n_k * (N - n_k))
    s_k = (a_k - b_k) / max(a_k, b_k)
    """
    S = _normalize_similarity_matrix(similarity)
    L = np.asarray(labels, dtype=np.int64).copy()

    if S.shape[0] != S.shape[1]:
        raise ValueError("similarity matrix must be square.")
    n_voxels = int(S.shape[0])
    if L.shape[0] != n_voxels:
        raise ValueError("label length must match similarity matrix dimension.")

    for invalid in unassigned_labels:
        L[L == int(invalid)] = 0

    uniq = np.unique(L[L > 0])
    if uniq.size == 0:
        return {
            "overall": float("nan"),
            "per_cluster": {},
            "a": {},
            "b": {},
            "n_k": {},
        }

    one = np.ones(n_voxels, dtype=np.float64)
    if issparse(S):
        S_one = np.asarray(S @ one).ravel()
    else:
        S_one = S @ one

    per_cluster: dict[int, float] = {}
    a: dict[int, float] = {}
    b: dict[int, float] = {}
    n_k_out: dict[int, int] = {}

    for k in uniq:
        mask = L == int(k)
        n_k = int(np.count_nonzero(mask))
        n_k_out[int(k)] = n_k
        if n_k <= 1 or n_k >= n_voxels:
            a[int(k)] = float("nan")
            b[int(k)] = float("nan")
            per_cluster[int(k)] = float("nan")
            continue

        z = mask.astype(np.float64, copy=False)
        if issparse(S):
            Sz = np.asarray(S @ z).ravel()
        else:
            Sz = S @ z

        num_a = float(np.dot(z, Sz))
        num_b = float(np.dot(z, S_one) - num_a)
        den_a = float(n_k * (n_k - 1))
        den_b = float(n_k * (n_voxels - n_k))

        a_k = num_a / den_a
        b_k = num_b / den_b
        denom = max(a_k, b_k)
        s_k = (a_k - b_k) / denom if np.isfinite(denom) and denom != 0.0 else float("nan")

        a[int(k)] = a_k
        b[int(k)] = b_k
        per_cluster[int(k)] = s_k

    valid_scores = np.array(list(per_cluster.values()), dtype=np.float64)
    overall = float(np.nanmean(valid_scores)) if np.isfinite(valid_scores).any() else float("nan")

    return {
        "overall": overall,
        "per_cluster": per_cluster,
        "a": a,
        "b": b,
        "n_k": n_k_out,
    }


def build_silhouette_report(
    clustering_result: ClusteringResult,
    *,
    unassigned_labels: Sequence[int] = (-1, 0),
) -> dict[str, Any]:
    """Compute subject-wise and aggregate silhouette report from ClusteringResult."""
    if clustering_result.similarity_matrices is None:
        raise ValueError(
            "clustering_result.similarity_matrices is required for silhouette reporting."
        )

    labels_2d = _labels_2d(clustering_result.label).astype(np.int64, copy=False)
    n_subjects = labels_2d.shape[0]
    if len(clustering_result.similarity_matrices) != n_subjects:
        raise ValueError(
            "number of similarity matrices must match number of subjects in labels."
        )

    details: list[dict[str, Any]] = []
    per_subject = np.empty(n_subjects, dtype=np.float64)

    for subject_index in range(n_subjects):
        subject_detail = _silhouette_cluster_strict_similarity(
            clustering_result.similarity_matrices[subject_index],
            labels_2d[subject_index],
            unassigned_labels=unassigned_labels,
        )
        details.append(subject_detail)
        per_subject[subject_index] = float(subject_detail["overall"])

    mean = float(np.nanmean(per_subject)) if np.isfinite(per_subject).any() else float("nan")
    std = float(np.nanstd(per_subject, ddof=1)) if np.isfinite(per_subject).sum() >= 2 else float("nan")

    return {
        "per_subject": per_subject.tolist(),
        "mean": mean,
        "sd": std,
        "per_subject_details": details,
    }


__all__ = ["build_silhouette_report"]
