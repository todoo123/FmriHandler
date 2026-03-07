from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from data.fmri import FmriVoxelData
from viz import build_silhouette_report, plot_subject_clusters_3d


def initial_pipeline_state() -> dict[str, Any]:
    return {
        "raw": None,
        "prep": None,
        "sim": None,
        "clu": None,
        "reports": {},
        "data_meta": {},
        "history": [],
    }


def run_validation(
    state: dict[str, Any],
    data_file: str,
    voxel_coord_file: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    data = _load_numpy_array(data_file, preferred_keys=("data", "arr_0"))
    voxel_coord = _load_numpy_array(
        voxel_coord_file,
        preferred_keys=("voxel_coord", "coords", "arr_0"),
    )

    raw = FmriVoxelData(data=data, voxel_coord=voxel_coord).validate()

    state["raw"] = raw
    state["prep"] = None
    state["sim"] = None
    state["clu"] = None
    state["reports"] = {}
    state["data_meta"] = {
        "data_shape": tuple(int(x) for x in data.shape),
        "voxel_coord_shape": tuple(int(x) for x in voxel_coord.shape),
        "data_dtype": str(data.dtype),
        "coord_dtype": str(voxel_coord.dtype),
    }
    _add_history(state, "validate", {"ok": True})
    return state, state["data_meta"]


def run_preprocess_similarity(
    state: dict[str, Any],
    preprocess_method: str,
    preprocess_hyper: dict[str, Any],
    similarity_method: str,
    similarity_hyper: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw = state.get("raw")
    if raw is None:
        raise ValueError("Validation step is required before preprocess/similarity.")

    prep = raw.preprocess(preprocess_method=preprocess_method, preprocess_hyper=preprocess_hyper)
    sim = prep.compute_similarity(
        similarity_method=similarity_method,
        similarity_hyper=similarity_hyper,
    )
    state["prep"] = prep
    state["sim"] = sim
    state["clu"] = None
    state["reports"] = {}
    summary = {
        "preprocess_method": preprocess_method,
        "preprocess_hyper": preprocess_hyper,
        "similarity_method": similarity_method,
        "similarity_hyper_resolved": sim.similarity_hyper,
        "n_similarity_matrices": 0 if sim.similarity_matrices is None else len(sim.similarity_matrices),
    }
    _add_history(state, "preprocess_similarity", summary)
    return state, summary


def run_clustering(
    state: dict[str, Any],
    clustering_method: str,
    cluster_hyper: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    sim = state.get("sim")
    if sim is None:
        raise ValueError("Similarity step is required before clustering.")

    clu = sim.cluster(cluster_method=clustering_method, cluster_hyper=cluster_hyper)
    state["clu"] = clu
    state["reports"] = {}
    summary = clu.summary()
    summary["cluster_hyper_resolved"] = clu.cluster_hyper
    _add_history(state, "cluster", summary)
    return state, summary


def run_report(
    state: dict[str, Any],
    report_method: str,
    report_hyper: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    clu = state.get("clu")
    if clu is None:
        raise ValueError("Clustering step is required before reporting.")

    if report_method != "silhouette":
        raise ValueError(f"Unknown report method: {report_method}")

    unassigned = report_hyper.get("unassigned_labels", (-1, 0))
    report = build_silhouette_report(clu, unassigned_labels=tuple(unassigned))
    state["reports"][report_method] = report
    _add_history(state, "report", {"method": report_method})
    return state, report


def run_visualization(
    state: dict[str, Any],
    viz_method: str,
    viz_hyper: dict[str, Any],
):
    clu = state.get("clu")
    if clu is None:
        raise ValueError("Clustering step is required before visualization.")
    if viz_method != "subject_clusters_3d":
        raise ValueError(f"Unknown visualization method: {viz_method}")

    figure = plot_subject_clusters_3d(
        clu,
        subject_index=int(viz_hyper.get("subject_index", 0)),
        include_unassigned=bool(viz_hyper.get("include_unassigned", True)),
        unassigned_labels=tuple(viz_hyper.get("unassigned_labels", (-1, 0))),
        marker_size=float(viz_hyper.get("marker_size", 2.2)),
        marker_opacity=float(viz_hyper.get("marker_opacity", 0.9)),
        title=(viz_hyper.get("title") or None),
    )
    _add_history(state, "visualization", {"method": viz_method})
    return state, figure


def pipeline_snapshot(state: dict[str, Any]) -> dict[str, Any]:
    raw = state.get("raw")
    prep = state.get("prep")
    sim = state.get("sim")
    clu = state.get("clu")
    return {
        "raw_ready": raw is not None,
        "prep_ready": prep is not None,
        "sim_ready": sim is not None,
        "clu_ready": clu is not None,
        "reports": list((state.get("reports") or {}).keys()),
        "data_meta": state.get("data_meta", {}),
        "history_count": len(state.get("history", [])),
    }


def _load_numpy_array(file_path: str, preferred_keys: tuple[str, ...]) -> np.ndarray:
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    suffix = path.suffix.lower()
    if suffix == ".npy":
        return np.load(path)

    if suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            for key in preferred_keys:
                if key in archive:
                    return np.asarray(archive[key])
            if len(archive.files) == 1:
                return np.asarray(archive[archive.files[0]])
            raise ValueError(
                f"Could not select array from npz: {path}. "
                f"Provide one of keys={preferred_keys} or a single-array archive."
            )

    raise ValueError(f"Unsupported file format: {path}. Use .npy or .npz")


def _add_history(state: dict[str, Any], step: str, payload: dict[str, Any]) -> None:
    history = state.setdefault("history", [])
    history.append(
        {
            "step": step,
            "timestamp": datetime.utcnow().isoformat(timespec="seconds") + "Z",
            "payload": payload,
        }
    )
