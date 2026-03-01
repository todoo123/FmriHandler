from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Tuple
import importlib

import numpy as np

from .preprocessed import PreprocessedData


ArrayLike = np.ndarray


def _get_package_root() -> str:
    """
    Return the top-level package name, e.g. 'FmriHandler'.

    This makes import paths robust even if the package is renamed.
    """
    name = __name__.split(".")[0]
    return name or "FmriHandler"


def _resolve_callable(dotted_path: str):
    """
    Resolve a 'pkg.module:function' style path to a Python callable.

    The path is interpreted as relative to the top-level package if it
    does not already start with it.
    """
    package_root = _get_package_root()

    if ":" in dotted_path:
        module_path, func_name = dotted_path.split(":", 1)
    else:
        # Fallback: last component after dot is treated as function name
        module_path, func_name = dotted_path.rsplit(".", 1)

    if not module_path.startswith(package_root):
        module_path = f"{package_root}.{module_path.lstrip('.')}"

    module = importlib.import_module(module_path)
    try:
        func = getattr(module, func_name)
    except AttributeError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            f"Engine function '{func_name}' not found in module '{module_path}'. "
            "Make sure the corresponding transform function is implemented."
        ) from exc

    if not callable(func):  # pragma: no cover - defensive
        raise TypeError(f"Resolved object '{dotted_path}' is not callable.")

    return func


# Mapping from high-level preprocessing method name → engine callable path.
# The actual functions are expected to live in the `transform/` layer and
# must not import from `data/` to avoid circular dependencies.
#
# Each engine function should have a signature compatible with:
#     def func(X: ArrayLike, coords: ArrayLike, **kwargs) -> Tuple[ArrayLike, Dict[str, Any]]:
#         ...
#
# where:
#   - X:       (T, V) fMRI time-by-voxel matrix
#   - coords:  (V, 3) voxel (i, j, k) coordinates
#   - return:  (features, metadata)
#
# You can extend or override this mapping as you implement the transform
# layer, but keep method names stable because they form the public API.
_PREPROCESS_DISPATCH: Dict[str, str] = {
    # Time-domain smoothing
    "gaussian": "transform._smoothing:gaussian_smooth_timeseries",
    # Basis expansion (e.g. B-spline)
    "bspline": "transform._basis_expansion:bspline_expand_timeseries",
    # Functional PCA
    "fpca": "transform._fpca:fpca_project_timeseries",
}


@dataclass(frozen=True)
class FMRIData:
    """
    Container for raw fMRI time series and voxel metadata.

    This is the entry point of the pipeline:

        FMRIData  →  PreprocessedData  →  SimilarityResult  →  ClusteringResult

    The class itself performs only light validation and dispatches heavy
    numerical work to engine functions in `transform/`.
    """

    data: ArrayLike
    coords: ArrayLike
    tr: Optional[float] = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        data = np.asarray(self.data)
        coords = np.asarray(self.coords)

        if data.ndim != 2:
            raise ValueError(f"`data` must be 2D (T, V); got shape {data.shape!r}")
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(
                "`coords` must be 2D with shape (V, 3) containing voxel (i, j, k) "
                f"coordinates; got shape {coords.shape!r}"
            )
        if data.shape[1] != coords.shape[0]:
            raise ValueError(
                "Mismatch between number of voxels in `data` and `coords`: "
                f"{data.shape[1]} vs {coords.shape[0]}"
            )

        # Normalize internal representation
        object.__setattr__(self, "data", data)
        object.__setattr__(self, "coords", coords)

    # ------------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------------
    @property
    def n_timepoints(self) -> int:
        return int(self.data.shape[0])

    @property
    def n_voxels(self) -> int:
        return int(self.data.shape[1])

    # ------------------------------------------------------------------
    # Validation & utilities
    # ------------------------------------------------------------------
    def validate(self) -> "FMRIData":
        """
        Run inexpensive sanity checks and return self for chaining.

        This method intentionally does not mutate data; it raises if
        basic assumptions are violated.
        """
        if not np.isfinite(self.data).all():
            raise ValueError("`data` contains non-finite values (NaN or inf).")

        if not np.isfinite(self.coords).all():
            raise ValueError("`coords` contains non-finite values (NaN or inf).")

        return self

    # ------------------------------------------------------------------
    # Transform layer dispatch
    # ------------------------------------------------------------------
    def preprocess(
        self,
        method: str,
        /,
        **kwargs: Any,
    ) -> PreprocessedData:
        """
        Apply a preprocessing / transform method and return `PreprocessedData`.

        Parameters
        ----------
        method:
            High-level preprocessing method name, e.g.:

            - ``\"gaussian\"``  – Gaussian temporal smoothing
            - ``\"bspline\"``   – B-spline basis expansion
            - ``\"fpca\"``      – functional PCA projection

        **kwargs:
            Passed through to the underlying engine function registered
            in ``_PREPROCESS_DISPATCH``.

        Returns
        -------
        PreprocessedData
            New object containing transformed features and metadata.
        """
        key = method.lower()
        if key not in _PREPROCESS_DISPATCH:
            raise ValueError(
                f"Unknown preprocess method {method!r}. "
                f"Known methods: {sorted(_PREPROCESS_DISPATCH)}"
            )

        engine_path = _PREPROCESS_DISPATCH[key]
        engine = _resolve_callable(engine_path)

        features, info = engine(self.data, self.coords, **kwargs)
        features = np.asarray(features)

        if features.ndim != 2:
            raise ValueError(
                "Preprocess engine must return a 2D feature array of shape (V, D); "
                f"got shape {features.shape!r}"
            )
        if features.shape[0] != self.n_voxels:
            raise ValueError(
                "Preprocess engine returned a feature matrix whose first dimension "
                "does not match the number of voxels. "
                f"Expected {self.n_voxels}, got {features.shape[0]}"
            )

        meta: Dict[str, Any] = dict(self.metadata)
        meta.update(
            {
                "preprocess_method": key,
                "preprocess_params": dict(kwargs),
                "preprocess_info": dict(info or {}),
            }
        )

        return PreprocessedData(
            features=features,
            coords=self.coords,
            fmri_metadata=meta,
        )

