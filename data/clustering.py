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
        n_assigned = int(self.label.size)
        n_clusters = int(np.unique(self.label).size) if n_assigned > 0 else 0

        return {
            "n_subjects": n_subjects,
            "n_voxels": n_voxels,
            "n_timepoints": n_timepoints,
            "preprocess_method": self.preprocess_method,
            "similarity_method": self.similarity_method,
            "n_labels": n_assigned,
            "n_clusters": n_clusters,
        }


__all__ = ["ClusteringResult"]
