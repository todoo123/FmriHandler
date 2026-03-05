from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np


@dataclass(slots=True)
class FmriVoxelData:
    """Layer 1 state object for raw fMRI voxel time-series."""

    data: np.ndarray
    voxel_coord: np.ndarray

    def validate(self) -> "FmriVoxelData":
        """Validate minimal shape constraints for MVP pipeline."""
        if not isinstance(self.data, np.ndarray) or self.data.ndim != 3:
            raise ValueError("data must be a 3D numpy array with shape (S, V, T).")
        if not isinstance(self.voxel_coord, np.ndarray) or self.voxel_coord.ndim != 2:
            raise ValueError("voxel_coord must be a 2D numpy array with shape (V, C).")
        if self.voxel_coord.shape[0] != self.data.shape[1]:
            raise ValueError("voxel_coord first dimension must match V in data (S, V, T).")
        return self

    def preprocess(
        self,
        preprocess_method: str,
        preprocess_hyper: Mapping[str, Any] | None = None,
    ) -> "PreprocessedData":
        """Dispatch to layer-2 object with method metadata."""
        from .preprocessed import PreprocessedData

        return PreprocessedData(
            data=self.data,
            voxel_coord=self.voxel_coord,
            preprocess_method=preprocess_method,
            preprocess_hyper=dict(preprocess_hyper or {}),
        )


__all__ = ["FmriVoxelData"]
