from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(slots=True)
class ClusteringResult:
    """Layer 4 state object for clustering outputs."""

    data: np.ndarray
    voxel_coord: np.ndarray
    preprocess_method: str
    preprocess_hyper: dict[str, Any] = field(default_factory=dict)
    similarity_method: str = ""
    similarity_hyper: dict[str, Any] = field(default_factory=dict)
    label: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    cluster_hyper: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> dict[str, Any]:
        """Return minimal metadata snapshot for quick checks."""
        n_subjects, n_voxels, n_timepoints = self.data.shape
        if self.label.ndim == 2:
            assigned_mask = self.label >= 0
            n_assigned_per_subject = np.count_nonzero(assigned_mask, axis=1).astype(int)
            n_clusters_per_subject = np.array(
                [
                    np.unique(self.label[s][assigned_mask[s]]).size
                    if n_assigned_per_subject[s] > 0
                    else 0
                    for s in range(self.label.shape[0])
                ],
                dtype=int,
            )
            n_assigned = int(n_assigned_per_subject.sum())
        else:
            assigned_mask = self.label >= 0
            n_assigned = int(np.count_nonzero(assigned_mask))
            n_clusters_per_subject = np.array(
                [int(np.unique(self.label[assigned_mask]).size) if n_assigned > 0 else 0],
                dtype=int,
            )
            n_assigned_per_subject = np.array([n_assigned], dtype=int)

        return {
            "n_subjects": n_subjects,
            "n_voxels": n_voxels,
            "n_timepoints": n_timepoints,
            "preprocess_method": self.preprocess_method,
            "similarity_method": self.similarity_method,
            "label_shape": tuple(self.label.shape),
            "n_labels_total": n_assigned,
            "n_labels_per_subject": n_assigned_per_subject.tolist(),
            "n_clusters_per_subject": n_clusters_per_subject.tolist(),
        }


__all__ = ["ClusteringResult"]
