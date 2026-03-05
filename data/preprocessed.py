from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


@dataclass(slots=True)
class PreprocessedData:
    """Layer 2 state object for preprocessed fMRI data."""

    data: np.ndarray
    voxel_coord: np.ndarray
    preprocess_method: str
    preprocess_hyper: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Apply the specified preprocessing method upon initialization."""
        method = self.preprocess_method.lower()
        if method == "gaussian_smoothing":
            from transform import gaussian_smoothing_1d
            self.data = gaussian_smoothing_1d(self.data, **self.preprocess_hyper)
        elif method in ("none", "pass"):
            pass
        else:
            raise ValueError(f"Unknown or unimplemented preprocess method: {self.preprocess_method}")

    def compute_similarity(
        self,
        similarity_method: str,
        similarity_hyper: Mapping[str, Any] | None = None,
    ) -> "SimilarityResult":
        """Dispatch to layer-3 object with similarity metadata."""
        from .similarity import SimilarityResult

        return SimilarityResult(
            data=self.data,
            voxel_coord=self.voxel_coord,
            preprocess_method=self.preprocess_method,
            preprocess_hyper=dict(self.preprocess_hyper),
            similarity_method=similarity_method,
            similarity_hyper=dict(similarity_hyper or {}),
        )


__all__ = ["PreprocessedData"]
