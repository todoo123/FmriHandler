from __future__ import annotations

from dataclasses import dataclass, field
import inspect
from typing import TYPE_CHECKING, Any, Mapping

import numpy as np
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    from .clustering import ClusteringResult


@dataclass(slots=True)
class SimilarityResult:
    """Layer 3 state object for similarity computation outputs."""

    data: np.ndarray
    voxel_coord: np.ndarray
    preprocess_method: str
    preprocess_hyper: dict[str, Any] = field(default_factory=dict)
    similarity_method: str = ""
    similarity_hyper: dict[str, Any] = field(default_factory=dict)
    similarity_matrices: list[csr_matrix] | None = None

    def __post_init__(self) -> None:
        """Run selected similarity computation at initialization time."""
        method = self.similarity_method.lower()
        if method == "correlation":
            from similarity._correlation import (
                compute_subjectwise_correlation_adjacency,
                resolve_correlation_hyper,
            )

            resolved_hyper = resolve_correlation_hyper(self.similarity_hyper)
            self.similarity_hyper = resolved_hyper
            self.similarity_matrices = compute_subjectwise_correlation_adjacency(
                self.data,
                self.voxel_coord,
                **resolved_hyper,
            )
        elif method in ("none", "pass", ""):
            self.similarity_matrices = None
        else:
            raise ValueError(
                f"Unknown or unimplemented similarity method: {self.similarity_method}"
            )

    def cluster(
        self,
        cluster_method: str,
        cluster_hyper: Mapping[str, Any] | None = None,
    ) -> "ClusteringResult":
        """Dispatch clustering engine and return layer-4 state object."""
        from .clustering import ClusteringResult

        resolved_cluster_hyper = dict(cluster_hyper or {})
        method = cluster_method.lower()
        resolved_cluster_hyper.setdefault("method", method)

        if method in ("spectral", "spectral_ncut", "spectral_unconstrained"):
            if not self.similarity_matrices:
                raise ValueError(
                    "subject-wise similarity matrices are required for spectral clustering. "
                    "Run compute_similarity(...) with a supported method first."
                )

            from cluster._spectral import resolve_spectral_hyper, spectral_clustering_ncut

            resolved_cluster_hyper = resolve_spectral_hyper(resolved_cluster_hyper)
            spectral_call_hyper = {
                k: v for k, v in resolved_cluster_hyper.items() if k != "method"
            }
            # Runtime compatibility guard: older/newer spectral implementations
            # may expose different keyword signatures.
            sig = inspect.signature(spectral_clustering_ncut)
            accepted = {name for name, p in sig.parameters.items() if p.kind != inspect.Parameter.VAR_KEYWORD}
            has_var_kwargs = any(
                p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()
            )
            if not has_var_kwargs:
                spectral_call_hyper = {
                    k: v for k, v in spectral_call_hyper.items() if k in accepted
                }
            labels = np.stack(
                [
                    spectral_clustering_ncut(
                        subject_similarity,
                        **spectral_call_hyper,
                    )
                    for subject_similarity in self.similarity_matrices
                ],
                axis=0,
            )
        else:
            raise ValueError(
                f"Unknown or unimplemented clustering method: {cluster_method}"
            )

        return ClusteringResult(
            data=self.data,
            voxel_coord=self.voxel_coord,
            preprocess_method=self.preprocess_method,
            preprocess_hyper=dict(self.preprocess_hyper),
            similarity_method=self.similarity_method,
            similarity_hyper=dict(self.similarity_hyper),
            similarity_matrices=(
                None
                if self.similarity_matrices is None
                else list(self.similarity_matrices)
            ),
            label=labels,
            cluster_hyper=resolved_cluster_hyper,
        )

__all__ = ["SimilarityResult"]
