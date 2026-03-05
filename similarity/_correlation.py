from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


def _default_workers() -> int:
    cpu = os.cpu_count() or 1
    return max(1, cpu - 1)


def _ensure_voxel_time_matrix(
    data: np.ndarray,
) -> np.ndarray:
    """Validate and convert single-subject (V, T) matrix for similarity."""
    if not isinstance(data, np.ndarray) or data.ndim != 2:
        raise ValueError("data must be a 2D numpy array with shape (V, T).")
    return np.asarray(data, dtype=np.float64)


def _ensure_subject_tensor(data: np.ndarray) -> np.ndarray:
    """Validate multi-subject tensor with shape (S, V, T)."""
    if not isinstance(data, np.ndarray) or data.ndim != 3:
        raise ValueError("data must be a 3D numpy array with shape (S, V, T).")
    return np.asarray(data, dtype=np.float64)


def _row_standardize(X: np.ndarray, *, eps: float = 1e-8) -> np.ndarray:
    """Row-wise standardization for (V, T) matrix."""
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = eps
    return (X - mu) / sd


def _upper_block_pairs(n_voxels: int, block_size: int) -> list[tuple[int, int, int, int]]:
    """Generate upper-triangle block pairs: (i0, i1, j0, j1) with bi <= bj."""
    starts = list(range(0, n_voxels, block_size))
    blocks: list[tuple[int, int, int, int]] = []
    for bi, i0 in enumerate(starts):
        i1 = min(i0 + block_size, n_voxels)
        for bj in range(bi, len(starts)):
            j0 = starts[bj]
            j1 = min(j0 + block_size, n_voxels)
            blocks.append((i0, i1, j0, j1))
    return blocks


def _corr_block_pair(
    Z: np.ndarray,
    i0: int,
    i1: int,
    j0: int,
    j1: int,
    denom: float,
    r_min: float | None,
    clamp_negative: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute correlations for one upper-triangle block pair."""
    zi = Z[i0:i1]  # (Bi, T)
    zj = Z[j0:j1]  # (Bj, T)
    corr = (zi @ zj.T) / denom  # (Bi, Bj)

    if i0 == j0:
        rr, cc = np.triu_indices(i1 - i0, k=1)
        vals = corr[rr, cc]
    else:
        rr, cc = np.indices(corr.shape)
        rr = rr.ravel()
        cc = cc.ravel()
        vals = corr.ravel()

    if clamp_negative:
        vals = np.maximum(vals, 0.0)
    if r_min is not None:
        vals[vals < float(r_min)] = 0.0

    keep = np.isfinite(vals) & (vals != 0.0)
    if not np.any(keep):
        return (
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.float64),
        )

    rows = (rr[keep] + i0).astype(np.int64, copy=False)
    cols = (cc[keep] + j0).astype(np.int64, copy=False)
    weights = vals[keep].astype(np.float64, copy=False)
    return rows, cols, weights


def compute_correlation_adjacency(
    data: np.ndarray,
    voxel_coord: np.ndarray,
    *,
    prestandardized: bool = False,
    block_size: int = 512,
    chunk_size: int | None = None,
    n_jobs: int | None = None,
    r_min: float | None = None,
    clamp_negative: bool = False,
) -> csr_matrix:
    """
    Compute sparse weighted adjacency from all-pairs voxel correlations.

    Returns a symmetric `scipy.sparse.csr_matrix` of shape (V, V).
    """
    X = _ensure_voxel_time_matrix(data)
    n_voxels, t_len = X.shape
    if n_voxels != voxel_coord.shape[0]:
        raise ValueError("voxel_coord first dimension must match V in data.")
    if t_len < 2:
        raise ValueError("At least 2 time points are required to compute correlation.")

    if n_voxels < 2:
        return csr_matrix((n_voxels, n_voxels), dtype=np.float64)

    Z = X if prestandardized else _row_standardize(X)
    denom = float(t_len - 1)
    # Backward compatibility: allow chunk_size as alias of block_size.
    if chunk_size is not None:
        block_size = int(chunk_size)
    blk = max(1, int(block_size))
    block_pairs = _upper_block_pairs(n_voxels, blk)
    max_workers = _default_workers() if n_jobs is None else max(1, int(n_jobs))

    rows_parts: list[np.ndarray] = []
    cols_parts: list[np.ndarray] = []
    vals_parts: list[np.ndarray] = []

    if max_workers == 1 or len(block_pairs) == 1:
        for i0, i1, j0, j1 in block_pairs:
            r, c, v = _corr_block_pair(
                Z, i0, i1, j0, j1, denom, r_min, clamp_negative
            )
            if v.size > 0:
                rows_parts.append(r)
                cols_parts.append(c)
                vals_parts.append(v)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = [
                ex.submit(
                    _corr_block_pair,
                    Z,
                    i0,
                    i1,
                    j0,
                    j1,
                    denom,
                    r_min,
                    clamp_negative,
                )
                for i0, i1, j0, j1 in block_pairs
            ]
            for fut in futures:
                r, c, v = fut.result()
                if v.size > 0:
                    rows_parts.append(r)
                    cols_parts.append(c)
                    vals_parts.append(v)

    if not vals_parts:
        return csr_matrix((n_voxels, n_voxels), dtype=np.float64)

    upper_i = np.concatenate(rows_parts)
    upper_j = np.concatenate(cols_parts)
    upper_w = np.concatenate(vals_parts)

    rows = np.concatenate([upper_i, upper_j])
    cols = np.concatenate([upper_j, upper_i])
    vals = np.concatenate([upper_w, upper_w])

    adj = coo_matrix((vals, (rows, cols)), shape=(n_voxels, n_voxels))
    return adj.tocsr()


def compute_subjectwise_correlation_adjacency(
    data: np.ndarray,
    voxel_coord: np.ndarray,
    *,
    prestandardized: bool = False,
    block_size: int = 512,
    chunk_size: int | None = None,
    n_jobs: int | None = None,
    r_min: float | None = None,
    clamp_negative: bool = False,
) -> list[csr_matrix]:
    """
    Compute one sparse correlation adjacency per subject.

    Returns a list of length S, each entry shape (V, V).
    """
    X = _ensure_subject_tensor(data)
    n_subjects = X.shape[0]
    out: list[csr_matrix] = []
    for s in range(n_subjects):
        out.append(
            compute_correlation_adjacency(
                X[s],
                voxel_coord,
                prestandardized=prestandardized,
                block_size=block_size,
                chunk_size=chunk_size,
                n_jobs=n_jobs,
                r_min=r_min,
                clamp_negative=clamp_negative,
            )
        )
    return out


def resolve_correlation_hyper(similarity_hyper: dict[str, Any]) -> dict[str, Any]:
    """Normalize optional hyperparameters for subject-wise correlation."""
    resolved = dict(similarity_hyper)
    resolved.setdefault("prestandardized", False)
    if "chunk_size" in resolved and "block_size" not in resolved:
        resolved["block_size"] = resolved.pop("chunk_size")
    resolved.setdefault("block_size", 512)
    resolved.setdefault("n_jobs", _default_workers())
    resolved.setdefault("r_min", None)
    resolved.setdefault("clamp_negative", False)
    return resolved


__all__ = [
    "compute_correlation_adjacency",
    "compute_subjectwise_correlation_adjacency",
    "resolve_correlation_hyper",
]
