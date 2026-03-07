from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

Stage = Literal["preprocess", "similarity", "clustering", "report", "visualization"]
InputKind = Literal["auto", "text", "number", "checkbox", "dropdown"]
ValueType = Literal["str", "int", "float", "bool", "csv_int", "json"]


@dataclass(slots=True)
class HyperParamSpec:
    """UI-facing specification for a single hyperparameter."""

    name: str
    label: str
    value_type: ValueType
    default: Any = None
    required: bool = False
    input_kind: InputKind = "auto"
    choices: list[str] = field(default_factory=list)
    advanced: bool = False
    help_text: str = ""


@dataclass(slots=True)
class MethodSpec:
    """UI-facing specification for one executable method."""

    key: str
    label: str
    description: str
    hyperparams: list[HyperParamSpec] = field(default_factory=list)


PREPROCESS_METHODS: dict[str, MethodSpec] = {
    "none": MethodSpec(
        key="none",
        label="No Preprocess",
        description="Use raw validated input as-is.",
        hyperparams=[],
    ),
    "gaussian_smoothing": MethodSpec(
        key="gaussian_smoothing",
        label="Gaussian Smoothing",
        description="Apply temporal Gaussian smoothing.",
        hyperparams=[
            HyperParamSpec(
                name="fwhm_samples",
                label="FWHM (samples)",
                value_type="float",
                default=3.0,
                required=True,
            ),
            HyperParamSpec(
                name="n_jobs",
                label="Parallel Jobs (-1: all cores)",
                value_type="int",
                default=-1,
                required=True,
            ),
            HyperParamSpec(
                name="truncate",
                label="Kernel Truncate (sigma units)",
                value_type="float",
                default=3.0,
                advanced=True,
            ),
            HyperParamSpec(
                name="chunk_size",
                label="Chunk Size",
                value_type="int",
                default=20000,
                advanced=True,
            ),
        ],
    ),
}

SIMILARITY_METHODS: dict[str, MethodSpec] = {
    "correlation": MethodSpec(
        key="correlation",
        label="Correlation",
        description="Subject-wise sparse correlation adjacency.",
        hyperparams=[
            HyperParamSpec(
                name="block_size",
                label="Block Size",
                value_type="int",
                default=512,
            ),
            HyperParamSpec(
                name="n_jobs",
                label="Parallel Jobs (empty: auto)",
                value_type="int",
                default=None,
            ),
            HyperParamSpec(
                name="r_min",
                label="Minimum Correlation Threshold (empty: none)",
                value_type="float",
                default=None,
            ),
            HyperParamSpec(
                name="clamp_negative",
                label="Clamp Negative Correlations",
                value_type="bool",
                default=False,
                input_kind="checkbox",
            ),
            HyperParamSpec(
                name="prestandardized",
                label="Input Already Standardized",
                value_type="bool",
                default=False,
                input_kind="checkbox",
                advanced=True,
            ),
        ],
    )
}

CLUSTERING_METHODS: dict[str, MethodSpec] = {
    "spectral": MethodSpec(
        key="spectral",
        label="Spectral Clustering",
        description="Ng-Jordan-Weiss style spectral clustering.",
        hyperparams=[
            HyperParamSpec(
                name="n_clusters",
                label="Number of Clusters",
                value_type="int",
                default=20,
                required=True,
            ),
            HyperParamSpec(
                name="eig_solver",
                label="Eigen Solver",
                value_type="str",
                default="auto",
                input_kind="dropdown",
                choices=["auto", "eigsh", "randomized"],
            ),
            HyperParamSpec(
                name="kmeans_nstart",
                label="KMeans Restarts",
                value_type="int",
                default=10,
            ),
            HyperParamSpec(
                name="kmeans_itermax",
                label="KMeans Max Iter",
                value_type="int",
                default=100,
            ),
            HyperParamSpec(
                name="random_state",
                label="Random Seed",
                value_type="int",
                default=42,
            ),
            HyperParamSpec(
                name="remove_isolated",
                label="Remove Isolated Nodes",
                value_type="bool",
                default=True,
                input_kind="checkbox",
                advanced=True,
            ),
            HyperParamSpec(
                name="enforce_nonnegative",
                label="Enforce Nonnegative Affinity",
                value_type="bool",
                default=True,
                input_kind="checkbox",
                advanced=True,
            ),
            HyperParamSpec(
                name="randomized_oversample",
                label="Randomized Oversample",
                value_type="int",
                default=20,
                advanced=True,
            ),
            HyperParamSpec(
                name="randomized_n_iter",
                label="Randomized Power Iterations",
                value_type="int",
                default=2,
                advanced=True,
            ),
            HyperParamSpec(
                name="randomized_switch_n_nodes",
                label="Randomized Switch Node Threshold",
                value_type="int",
                default=15000,
                advanced=True,
            ),
        ],
    )
}

REPORT_METHODS: dict[str, MethodSpec] = {
    "silhouette": MethodSpec(
        key="silhouette",
        label="Silhouette Report",
        description="Per-subject and aggregate silhouette metrics.",
        hyperparams=[
            HyperParamSpec(
                name="unassigned_labels",
                label="Unassigned Labels (CSV int)",
                value_type="csv_int",
                default="-1,0",
            )
        ],
    )
}

VISUALIZATION_METHODS: dict[str, MethodSpec] = {
    "subject_clusters_3d": MethodSpec(
        key="subject_clusters_3d",
        label="3D Subject Cluster Plot",
        description="3D scatter by subject cluster labels.",
        hyperparams=[
            HyperParamSpec(
                name="subject_index",
                label="Subject Index",
                value_type="int",
                default=0,
                required=True,
            ),
            HyperParamSpec(
                name="include_unassigned",
                label="Include Unassigned Labels",
                value_type="bool",
                default=True,
                input_kind="checkbox",
            ),
            HyperParamSpec(
                name="marker_size",
                label="Marker Size",
                value_type="float",
                default=2.2,
            ),
            HyperParamSpec(
                name="marker_opacity",
                label="Marker Opacity",
                value_type="float",
                default=0.9,
            ),
            HyperParamSpec(
                name="unassigned_labels",
                label="Unassigned Labels (CSV int)",
                value_type="csv_int",
                default="-1,0",
                advanced=True,
            ),
            HyperParamSpec(
                name="title",
                label="Plot Title (optional)",
                value_type="str",
                default="",
                advanced=True,
            ),
        ],
    )
}


METHOD_REGISTRY: dict[Stage, dict[str, MethodSpec]] = {
    "preprocess": PREPROCESS_METHODS,
    "similarity": SIMILARITY_METHODS,
    "clustering": CLUSTERING_METHODS,
    "report": REPORT_METHODS,
    "visualization": VISUALIZATION_METHODS,
}


def list_method_choices(stage: Stage) -> list[tuple[str, str]]:
    """Return Gradio dropdown-friendly choices as (label, key) tuples."""
    methods = METHOD_REGISTRY[stage]
    return [(spec.label, spec.key) for spec in methods.values()]


def get_method(stage: Stage, method_key: str) -> MethodSpec:
    """
    Resolve a method specification by stage and key.

    Extension guide:
    1) Add a new MethodSpec entry to the stage map above.
    2) Ensure the backend executor in ui/pipeline.py supports the same method key.
    3) The Gradio UI updates automatically from this registry.
    """
    methods = METHOD_REGISTRY[stage]
    if method_key not in methods:
        raise ValueError(f"Unknown {stage} method: {method_key}")
    return methods[method_key]
