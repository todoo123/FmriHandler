from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Union
import importlib

import numpy as np

from .clustering import ClusteringResult

try:  # optional dependency
    import scipy.sparse as sp  # type: ignore
except Exception:  # pragma: no cover
    sp = None  # type: ignore


ArrayLike = Union[np.ndarray, "sp.spmatrix"]  # type: ignore[name-defined]


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
            "Make sure the corresponding clustering function is implemented."
        ) from exc

    if not callable(func):  # pragma: no cover
        raise TypeError(f"Resolved object '{dotted_path}' is not callable.")

    return func


# Mapping from clustering method name → engine callable path in `cluster/`.
# Each engine function is expected to have a signature roughly like:
#     def func(matrix: ArrayLike, n_clusters: int, **kwargs) -> np.ndarray:
#         ...
#
# where:
#   - matrix:  (V, V) affinity or similarity matrix
#   - return:  integer labels of length V
_CLUSTER_DISPATCH: Dict[str, str] = {
    "spectral_craddock26": "cluster._spectral:spectral_craddock26",
    "spectral": "cluster._spectral:spectral_clustering",
    "hierarchical": "cluster._hierarchical:hierarchical_clustering",
    "kmeans": "cluster._kmeans:kmeans_from_similarity",
}


@dataclass(frozen=True)
class SimilarityResult:
    """
    Result of voxel-wise similarity / affinity computation.

    Attributes
    ----------
    matrix:
        (V, V) similarity or affinity matrix. May be a dense NumPy array
        or a SciPy sparse matrix.
    coords:
        Voxel coordinates with shape (V, 3).
    preprocessed_metadata:
        Metadata propagated from `PreprocessedData`.
    is_affinity:
        Whether the matrix is already an affinity matrix suitable for
        clustering (as opposed to a raw distance/similarity).
    similarity_method, similarity_params:
        Method name and arguments used to construct the matrix.
    """

    matrix: ArrayLike
    coords: np.ndarray
    preprocessed_metadata: Mapping[str, Any] = field(default_factory=dict)
    is_affinity: bool = True
    similarity_method: Optional[str] = None
    similarity_params: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        coords = np.asarray(self.coords)

        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(
                "`coords` must be 2D with shape (V, 3); "
                f"got shape {coords.shape!r}"
            )

        # Basic shape check for matrix (supports dense or sparse)
        mat = self.matrix
        if sp is not None and sp.issparse(mat):  # type: ignore[attr-defined]
            v0, v1 = mat.shape
        else:
            arr = np.asarray(mat)
            if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
                raise ValueError(
                    "`matrix` must be a square (V, V) array or sparse matrix; "
                    f"got shape {arr.shape!r}"
                )
            v0 = v1 = arr.shape[0]

        if v0 != coords.shape[0] or v1 != coords.shape[0]:
            raise ValueError(
                "Mismatch between matrix size and number of voxels in `coords`."
            )

        object.__setattr__(self, "coords", coords)

    # ------------------------------------------------------------------
    # Clustering layer dispatch
    # ------------------------------------------------------------------
    def cluster(
        self,
        method: str,
        /,
        *,
        n_clusters: int,
        **kwargs: Any,
    ) -> ClusteringResult:
        """
        Run a clustering algorithm on the similarity / affinity matrix.

        Parameters
        ----------
        method:
            Clustering method name, e.g. ``\"spectral_craddock26\"``,
            ``\"spectral\"``, ``\"hierarchical\"``, ``\"kmeans\"``.
        n_clusters:
            Desired number of clusters (ROIs).
        **kwargs:
            Passed directly to the underlying clustering engine function.
        """
        if n_clusters <= 0:
            raise ValueError("`n_clusters` must be a positive integer.")

        key = method.lower()
        if key not in _CLUSTER_DISPATCH:
            raise ValueError(
                f"Unknown clustering method {method!r}. "
                f"Known methods: {sorted(_CLUSTER_DISPATCH)}"
            )

        engine_path = _CLUSTER_DISPATCH[key]
        engine = _resolve_callable(engine_path)

        labels = engine(self.matrix, n_clusters=n_clusters, **kwargs)

        labels_arr = np.asarray(labels, dtype=int)
        if labels_arr.ndim != 1 or labels_arr.shape[0] != self.coords.shape[0]:
            raise ValueError(
                "Clustering engine must return a 1D array of length V with "
                "cluster labels for each voxel."
            )

        meta: Dict[str, Any] = dict(self.preprocessed_metadata)
        meta.update(
            {
                "similarity_method": self.similarity_method,
                "similarity_params": dict(self.similarity_params),
                "is_affinity": self.is_affinity,
                "cluster_method": key,
                "cluster_params": dict(kwargs),
                "n_clusters": int(n_clusters),
            }
        )

        return ClusteringResult(
            labels=labels_arr,
            coords=self.coords,
            similarity_matrix=self.matrix,
            metadata=meta,
        )

