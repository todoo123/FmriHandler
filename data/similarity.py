from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


@dataclass(slots=True)
class SimilarityResult:
    """Layer 3 state object for similarity computation outputs."""

    data: np.ndarray
    voxel_coord: np.ndarray
    preprocess_method: str
    preprocess_hyper: dict[str, Any] = field(default_factory=dict)
    similarity_method: str = ""
    similarity_hyper: dict[str, Any] = field(default_factory=dict)

    def cluster(
        self,
        cluster_method: str,
        cluster_hyper: Mapping[str, Any] | None = None,
    ) -> "ClusteringResult":
        """Dispatch to layer-4 object with placeholder labels for MVP."""
        from .clustering import ClusteringResult

        resolved_cluster_hyper = dict(cluster_hyper or {})
        resolved_cluster_hyper.setdefault("method", cluster_method)
        labels = np.full(self.data.shape[1], -1, dtype=int)

        return ClusteringResult(
            data=self.data,
            voxel_coord=self.voxel_coord,
            preprocess_method=self.preprocess_method,
            preprocess_hyper=dict(self.preprocess_hyper),
            similarity_method=self.similarity_method,
            similarity_hyper=dict(self.similarity_hyper),
            label=labels,
            cluster_hyper=resolved_cluster_hyper,
        )


__all__ = ["SimilarityResult"]
