from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np

try:  # optional visualization / metrics dependencies
    import matplotlib.pyplot as plt  # type: ignore
except Exception:  # pragma: no cover
    plt = None  # type: ignore


ArrayLike = np.ndarray


@dataclass(frozen=True)
class ClusteringResult:
    """
    Final clustering result over voxels.

    Attributes
    ----------
    labels:
        Integer cluster labels of shape (V,). Label semantics are up to
        the clustering engine (e.g. 0 can denote isolated / background).
    coords:
        Voxel coordinates with shape (V, 3).
    similarity_matrix:
        Optional similarity / affinity matrix used to derive the labels.
        Stored for diagnostics and metric computations.
    metadata:
        Aggregated metadata from previous layers plus clustering info.
    """

    labels: ArrayLike
    coords: ArrayLike
    similarity_matrix: Optional[Any] = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        labels = np.asarray(self.labels, dtype=int)
        coords = np.asarray(self.coords)

        if labels.ndim != 1:
            raise ValueError("`labels` must be a 1D array of length V.")
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(
                "`coords` must be 2D with shape (V, 3); "
                f"got shape {coords.shape!r}"
            )
        if labels.shape[0] != coords.shape[0]:
            raise ValueError(
                "Mismatch between number of labels and voxel coordinates: "
                f"{labels.shape[0]} vs {coords.shape[0]}"
            )

        object.__setattr__(self, "labels", labels)
        object.__setattr__(self, "coords", coords)

    # ------------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------------
    @property
    def n_voxels(self) -> int:
        return int(self.labels.shape[0])

    @property
    def n_clusters(self) -> int:
        return int(len(np.unique(self.labels)))

    # ------------------------------------------------------------------
    # Lightweight reporting helpers
    # ------------------------------------------------------------------
    def to_dataframe(self):
        """
        Return a pandas DataFrame combining coordinates and labels.

        This is a convenience method; it requires pandas to be installed.
        """
        try:
            import pandas as pd  # type: ignore
        except Exception as exc:  # pragma: no cover
            raise RuntimeError(
                "pandas is required to use `ClusteringResult.to_dataframe()`."
            ) from exc

        df = pd.DataFrame(
            {
                "i": self.coords[:, 0],
                "j": self.coords[:, 1],
                "k": self.coords[:, 2],
                "label": self.labels,
            }
        )
        return df

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------
    def plot_voxels_3d(
        self,
        *,
        elev: float = 30.0,
        azim: float = 30.0,
        figsize: tuple[float, float] = (6.0, 6.0),
        cmap: str = "tab20",
        alpha: float = 1.0,
        show: bool = True,
    ) -> Optional[Any]:
        """
        Minimal 3D scatter plot of voxel clusters using matplotlib.

        This method is intentionally simple so that the core package
        remains lightweight. For richer interactive plots, you can build
        separate helpers under `viz/` and call them from here.
        """
        if plt is None:  # pragma: no cover
            raise RuntimeError(
                "matplotlib is required for `plot_voxels_3d`, "
                "but it is not installed."
            )

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="3d")

        sc = ax.scatter(
            self.coords[:, 0],
            self.coords[:, 1],
            self.coords[:, 2],
            c=self.labels,
            cmap=cmap,
            s=4,
            alpha=alpha,
        )
        ax.set_xlabel("i")
        ax.set_ylabel("j")
        ax.set_zlabel("k")
        ax.view_init(elev=elev, azim=azim)
        fig.colorbar(sc, ax=ax, shrink=0.6, label="Cluster")

        if show:
            plt.show()
            return None
        return fig

    # ------------------------------------------------------------------
    # Metrics
    # ------------------------------------------------------------------
    def report_silhouette(self) -> Optional[float]:
        """
        Compute a simple silhouette-like score if a similarity matrix is available.

        This is a placeholder that can be replaced with a more faithful
        implementation in `metrics/`. For now it returns ``None`` if the
        required information is missing.
        """
        if self.similarity_matrix is None:
            return None

        # Lazy import to avoid hard dependency if metrics are unused.
        try:
            from ..metrics import silhouette_from_sparseW_labels  # type: ignore
        except Exception:
            return None

        try:
            res = silhouette_from_sparseW_labels(self.similarity_matrix, self.labels)
        except Exception:  # pragma: no cover
            return None

        overall = getattr(res, "overall", None)
        if isinstance(overall, dict):
            return overall.get("overall")  # type: ignore[no-any-return]
        return overall  # type: ignore[no-any-return]

