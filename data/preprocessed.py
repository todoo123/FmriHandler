from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Tuple, Union
import importlib

import numpy as np

from .similarity import SimilarityResult


ArrayLike = np.ndarray


def _get_package_root() -> str:
    name = __name__.split(".")[0]
    return name or "FmriHandler"


def _resolve_callable(dotted_path: str):
    package_root = _get_package_root()

    if ":" in dotted_path:
        module_path, func_name = dotted_path.split(":", 1)
    else:
        module_path, func_name = dotted_path.rsplit(".", 1)

    if not module_path.startswith(package_root):
        module_path = f"{package_root}.{module_path.lstrip('.')}"

    module = importlib.import_module(module_path)
    try:
        func = getattr(module, func_name)
    except AttributeError as exc:  # pragma: no cover
        raise RuntimeError(
            f"Engine function '{func_name}' not found in module '{module_path}'. "
            "Make sure the corresponding similarity function is implemented."
        ) from exc

    if not callable(func):  # pragma: no cover
        raise TypeError(f"Resolved object '{dotted_path}' is not callable.")

    return func


# Mapping from similarity method name → engine callable path in `similarity/`.
# Each engine function is expected to follow the convention
#     def func(features: ArrayLike, **kwargs) -> ArrayLike:
#         ...
#
# where:
#   - features: (V, D) matrix of voxel-wise features
#   - return:  (V, V) similarity or affinity matrix (dense or sparse-like)
_SIMILARITY_DISPATCH: Dict[str, str] = {
    "correlation": "similarity._correlation:correlation_matrix",
    "partial_correlation": "similarity._partial_corr:partial_correlation_matrix",
    "cosine": "similarity._cosine:cosine_similarity_matrix",
    "l2": "similarity._l2distance:l2_distance_matrix",
}


@dataclass(frozen=True)
class PreprocessedData:
    """
    Output of the transform layer for a single fMRI dataset.

    Attributes
    ----------
    features:
        Array of shape (V, D) where V is the number of voxels and
        D is the feature dimension produced by the chosen transform.
    coords:
        Voxel coordinates with shape (V, 3).
    fmri_metadata:
        Metadata propagated from the original `FMRIData` plus
        transform-specific information.
    """

    features: ArrayLike
    coords: ArrayLike
    fmri_metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        features = np.asarray(self.features)
        coords = np.asarray(self.coords)

        if features.ndim != 2:
            raise ValueError(
                "`features` must be 2D with shape (V, D); "
                f"got shape {features.shape!r}"
            )
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(
                "`coords` must be 2D with shape (V, 3); "
                f"got shape {coords.shape!r}"
            )
        if features.shape[0] != coords.shape[0]:
            raise ValueError(
                "Mismatch between number of voxels in `features` and `coords`: "
                f"{features.shape[0]} vs {coords.shape[0]}"
            )

        object.__setattr__(self, "features", features)
        object.__setattr__(self, "coords", coords)

    # ------------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------------
    @property
    def n_voxels(self) -> int:
        return int(self.features.shape[0])

    @property
    def n_features(self) -> int:
        return int(self.features.shape[1])

    # ------------------------------------------------------------------
    # Similarity layer dispatch
    # ------------------------------------------------------------------
    def compute_similarity(
        self,
        method: str,
        /,
        *,
        to_affinity: bool = True,
        **kwargs: Any,
    ) -> SimilarityResult:
        """
        Compute voxel-wise similarity or affinity matrix.

        Parameters
        ----------
        method:
            Name of the similarity measure, e.g. ``\"correlation\"``,
            ``\"partial_correlation\"``, ``\"cosine\"``, ``\"l2\"``.
        to_affinity:
            If ``True``, the engine is expected to return an affinity
            matrix suitable for clustering. If ``False``, the matrix
            may be a raw similarity or distance matrix.
        **kwargs:
            Passed directly to the underlying similarity engine function.

        Returns
        -------
        SimilarityResult
            Wrapper around the resulting matrix and metadata.
        """
        key = method.lower()
        if key not in _SIMILARITY_DISPATCH:
            raise ValueError(
                f"Unknown similarity method {method!r}. "
                f"Known methods: {sorted(_SIMILARITY_DISPATCH)}"
            )

        engine_path = _SIMILARITY_DISPATCH[key]
        engine = _resolve_callable(engine_path)

        matrix = engine(self.features, **kwargs)

        return SimilarityResult(
            matrix=matrix,
            coords=self.coords,
            preprocessed_metadata=dict(self.fmri_metadata),
            is_affinity=to_affinity,
            similarity_method=key,
            similarity_params=dict(kwargs),
        )

