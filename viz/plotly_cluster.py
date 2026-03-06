from __future__ import annotations

from typing import Sequence

import numpy as np
import plotly.graph_objects as go
from plotly.colors import qualitative

from data.clustering import ClusteringResult


def _labels_2d(label: np.ndarray) -> np.ndarray:
    arr = np.asarray(label)
    if arr.ndim == 1:
        return arr.reshape(1, -1)
    if arr.ndim == 2:
        return arr
    raise ValueError("label must be 1D or 2D array.")


def _validate_coords(voxel_coord: np.ndarray) -> np.ndarray:
    coords = np.asarray(voxel_coord, dtype=np.float64)
    if coords.ndim != 2:
        raise ValueError("voxel_coord must be a 2D array with shape (V, C).")
    if coords.shape[1] < 3:
        raise ValueError("voxel_coord must provide at least 3 columns for x, y, z.")
    return coords[:, :3]


def _build_label_palette(labels: np.ndarray) -> dict[int, str]:
    palette = (
        qualitative.Plotly
        + qualitative.D3
        + qualitative.G10
        + qualitative.T10
        + qualitative.Alphabet
    )
    uniq = sorted(int(v) for v in np.unique(labels))
    return {lab: palette[idx % len(palette)] for idx, lab in enumerate(uniq)}


def plot_subject_clusters_3d(
    clustering_result: ClusteringResult,
    subject_index: int = 0,
    *,
    include_unassigned: bool = True,
    unassigned_labels: Sequence[int] = (-1, 0),
    marker_size: float = 2.2,
    marker_opacity: float = 0.9,
    title: str | None = None,
) -> go.Figure:
    """
    Create a 3D scatter plot for one subject's clustering labels.

    One trace is created per cluster label so labels get distinct colors.
    """
    labels_2d = _labels_2d(clustering_result.label)
    n_subjects, n_voxels = labels_2d.shape
    if not (0 <= subject_index < n_subjects):
        raise IndexError(f"subject_index={subject_index} is out of range for {n_subjects} subjects.")

    coords = _validate_coords(clustering_result.voxel_coord)
    if coords.shape[0] != n_voxels:
        raise ValueError("voxel_coord row count must match number of voxel labels.")

    labels = labels_2d[subject_index].astype(np.int64, copy=False)
    if not include_unassigned:
        mask = np.ones(labels.shape[0], dtype=bool)
        for invalid in unassigned_labels:
            mask &= labels != int(invalid)
        coords = coords[mask]
        labels = labels[mask]

    label_colors = _build_label_palette(labels)
    fig = go.Figure()
    for lab in sorted(label_colors):
        m = labels == lab
        if not np.any(m):
            continue
        fig.add_trace(
            go.Scatter3d(
                x=coords[m, 0],
                y=coords[m, 1],
                z=coords[m, 2],
                mode="markers",
                name=f"cluster {lab}",
                marker={
                    "size": marker_size,
                    "opacity": marker_opacity,
                    "color": label_colors[lab],
                },
            )
        )

    fig.update_layout(
        title=title or f"Subject {subject_index} voxel clusters",
        scene={"xaxis_title": "x", "yaxis_title": "y", "zaxis_title": "z"},
        legend={"itemsizing": "constant"},
        margin={"l": 0, "r": 0, "b": 0, "t": 40},
    )
    return fig


def plot_all_subject_clusters_3d(
    clustering_result: ClusteringResult,
    *,
    include_unassigned: bool = True,
    unassigned_labels: Sequence[int] = (-1, 0),
    marker_size: float = 2.2,
    marker_opacity: float = 0.9,
    title_prefix: str = "Subject",
) -> list[go.Figure]:
    """Create one 3D scatter figure per subject."""
    labels_2d = _labels_2d(clustering_result.label)
    figures: list[go.Figure] = []
    for subject_index in range(labels_2d.shape[0]):
        figures.append(
            plot_subject_clusters_3d(
                clustering_result,
                subject_index=subject_index,
                include_unassigned=include_unassigned,
                unassigned_labels=unassigned_labels,
                marker_size=marker_size,
                marker_opacity=marker_opacity,
                title=f"{title_prefix} {subject_index} voxel clusters",
            )
        )
    return figures


__all__ = ["plot_subject_clusters_3d", "plot_all_subject_clusters_3d"]
