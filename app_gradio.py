from __future__ import annotations

import traceback
from typing import Any

import gradio as gr

from ui import (
    build_widget_updates,
    create_hyperparam_widgets,
    get_method,
    initial_pipeline_state,
    list_method_choices,
    parse_hyperparam_values,
    pipeline_snapshot,
    run_clustering,
    run_preprocess_similarity,
    run_report,
    run_validation,
    run_visualization,
    split_widget_values,
)

MAX_BASIC_PARAMS = 8
MAX_ADV_PARAMS = 8
WIDGET_WIDTH_BASIC = MAX_BASIC_PARAMS * 5
WIDGET_WIDTH_ADV = MAX_ADV_PARAMS * 5


def _filter_specs(method_key: str, stage: str, *, advanced: bool) -> list[Any]:
    method = get_method(stage, method_key)
    return [spec for spec in method.hyperparams if spec.advanced is advanced]


def _update_hyper_widgets(stage: str, method_key: str) -> list[dict[str, Any]]:
    basic_specs = _filter_specs(method_key, stage, advanced=False)
    adv_specs = _filter_specs(method_key, stage, advanced=True)
    return [
        *build_widget_updates(basic_specs, MAX_BASIC_PARAMS),
        *build_widget_updates(adv_specs, MAX_ADV_PARAMS),
    ]


def _parse_hyper_values(
    stage: str,
    method_key: str,
    basic_flat: list[Any],
    adv_flat: list[Any],
) -> dict[str, Any]:
    basic_specs = _filter_specs(method_key, stage, advanced=False)
    adv_specs = _filter_specs(method_key, stage, advanced=True)

    b_keys, b_texts, b_numbers, b_bools, b_choices = split_widget_values(
        basic_flat, MAX_BASIC_PARAMS
    )
    a_keys, a_texts, a_numbers, a_bools, a_choices = split_widget_values(
        adv_flat, MAX_ADV_PARAMS
    )

    parsed_basic = parse_hyperparam_values(
        basic_specs,
        b_keys,
        b_texts,
        b_numbers,
        b_bools,
        b_choices,
    )
    parsed_adv = parse_hyperparam_values(
        adv_specs,
        a_keys,
        a_texts,
        a_numbers,
        a_bools,
        a_choices,
    )
    return {**parsed_basic, **parsed_adv}


def _ok(message: str) -> str:
    return f"OK: {message}"


def _error(message: str) -> str:
    return f"ERROR: {message}"


def _format_exception(exc: Exception) -> str:
    return _error(f"{exc}\n\n{traceback.format_exc()}")


def _on_validate(
    state: dict[str, Any],
    data_file: str | None,
    voxel_coord_file: str | None,
) -> tuple[dict[str, Any], str, dict[str, Any], dict[str, Any]]:
    if not data_file or not voxel_coord_file:
        return state, _error("Both data and voxel_coord files are required."), {}, pipeline_snapshot(state)
    try:
        state, summary = run_validation(state, data_file, voxel_coord_file)
        return state, _ok("Validation complete."), summary, pipeline_snapshot(state)
    except Exception as exc:  # pragma: no cover - UI error path
        return state, _format_exception(exc), {}, pipeline_snapshot(state)


def _on_run_preprocess_similarity(
    state: dict[str, Any],
    preprocess_method: str,
    similarity_method: str,
    *flat_values: Any,
) -> tuple[dict[str, Any], str, dict[str, Any], dict[str, Any]]:
    try:
        p_basic = list(flat_values[0:WIDGET_WIDTH_BASIC])
        p_adv = list(flat_values[WIDGET_WIDTH_BASIC : WIDGET_WIDTH_BASIC + WIDGET_WIDTH_ADV])
        s_start = WIDGET_WIDTH_BASIC + WIDGET_WIDTH_ADV
        s_basic = list(flat_values[s_start : s_start + WIDGET_WIDTH_BASIC])
        s_adv = list(flat_values[s_start + WIDGET_WIDTH_BASIC : s_start + WIDGET_WIDTH_BASIC + WIDGET_WIDTH_ADV])

        preprocess_hyper = _parse_hyper_values("preprocess", preprocess_method, p_basic, p_adv)
        similarity_hyper = _parse_hyper_values("similarity", similarity_method, s_basic, s_adv)
        state, summary = run_preprocess_similarity(
            state,
            preprocess_method=preprocess_method,
            preprocess_hyper=preprocess_hyper,
            similarity_method=similarity_method,
            similarity_hyper=similarity_hyper,
        )
        return state, _ok("Preprocess and similarity complete."), summary, pipeline_snapshot(state)
    except Exception as exc:  # pragma: no cover - UI error path
        return state, _format_exception(exc), {}, pipeline_snapshot(state)


def _on_run_clustering(
    state: dict[str, Any],
    clustering_method: str,
    *flat_values: Any,
) -> tuple[dict[str, Any], str, dict[str, Any], dict[str, Any]]:
    try:
        c_basic = list(flat_values[0:WIDGET_WIDTH_BASIC])
        c_adv = list(flat_values[WIDGET_WIDTH_BASIC : WIDGET_WIDTH_BASIC + WIDGET_WIDTH_ADV])
        cluster_hyper = _parse_hyper_values("clustering", clustering_method, c_basic, c_adv)
        state, summary = run_clustering(state, clustering_method=clustering_method, cluster_hyper=cluster_hyper)
        return state, _ok("Clustering complete."), summary, pipeline_snapshot(state)
    except Exception as exc:  # pragma: no cover - UI error path
        return state, _format_exception(exc), {}, pipeline_snapshot(state)


def _on_run_report(
    state: dict[str, Any],
    report_method: str,
    *flat_values: Any,
) -> tuple[dict[str, Any], str, dict[str, Any], dict[str, Any]]:
    try:
        r_basic = list(flat_values[0:WIDGET_WIDTH_BASIC])
        r_adv = list(flat_values[WIDGET_WIDTH_BASIC : WIDGET_WIDTH_BASIC + WIDGET_WIDTH_ADV])
        report_hyper = _parse_hyper_values("report", report_method, r_basic, r_adv)
        state, report = run_report(state, report_method=report_method, report_hyper=report_hyper)
        return state, _ok("Report complete."), report, pipeline_snapshot(state)
    except Exception as exc:  # pragma: no cover - UI error path
        return state, _format_exception(exc), {}, pipeline_snapshot(state)


def _on_run_visualization(
    state: dict[str, Any],
    viz_method: str,
    *flat_values: Any,
):
    try:
        v_basic = list(flat_values[0:WIDGET_WIDTH_BASIC])
        v_adv = list(flat_values[WIDGET_WIDTH_BASIC : WIDGET_WIDTH_BASIC + WIDGET_WIDTH_ADV])
        viz_hyper = _parse_hyper_values("visualization", viz_method, v_basic, v_adv)
        state, figure = run_visualization(state, viz_method=viz_method, viz_hyper=viz_hyper)
        return state, figure, _ok("Visualization rendered."), pipeline_snapshot(state)
    except Exception as exc:  # pragma: no cover - UI error path
        return state, None, _format_exception(exc), pipeline_snapshot(state)


def build_demo() -> gr.Blocks:
    initial_state = initial_pipeline_state()
    preprocess_default = next(iter(list_method_choices("preprocess")))[1]
    similarity_default = next(iter(list_method_choices("similarity")))[1]
    clustering_default = next(iter(list_method_choices("clustering")))[1]
    report_default = next(iter(list_method_choices("report")))[1]
    viz_default = next(iter(list_method_choices("visualization")))[1]

    with gr.Blocks(title="fMRI Handler - Gradio App") as demo:
        gr.Markdown(
            "## fMRI Handler Pipeline UI\n"
            "Upload -> Validate -> Preprocess/Similarity -> Cluster/Report -> Visualize"
        )

        state = gr.State(initial_state)
        snapshot_box = gr.JSON(label="Pipeline Snapshot", value=pipeline_snapshot(initial_state))

        with gr.Tab("1) Upload & Validation"):
            data_file = gr.File(
                label="Data (.npy or .npz)",
                file_types=[".npy", ".npz"],
                type="filepath",
            )
            coord_file = gr.File(
                label="Voxel Coordinates (.npy or .npz)",
                file_types=[".npy", ".npz"],
                type="filepath",
            )
            validate_button = gr.Button("Run Validation", variant="primary")
            validation_status = gr.Textbox(label="Status", lines=6, interactive=False)
            validation_summary = gr.JSON(label="Validation Summary")

        with gr.Tab("2) Preprocess & Similarity"):
            preprocess_method = gr.Dropdown(
                choices=list_method_choices("preprocess"),
                value=preprocess_default,
                label="Preprocess Method",
            )
            gr.Markdown("### Preprocess Hyperparameters")
            preprocess_basic = create_hyperparam_widgets(
                max_params=MAX_BASIC_PARAMS,
                initial_specs=_filter_specs(preprocess_default, "preprocess", advanced=False),
            )
            with gr.Accordion("Preprocess Advanced Hyperparameters", open=False):
                preprocess_adv = create_hyperparam_widgets(
                    max_params=MAX_ADV_PARAMS,
                    initial_specs=_filter_specs(preprocess_default, "preprocess", advanced=True),
                )

            similarity_method = gr.Dropdown(
                choices=list_method_choices("similarity"),
                value=similarity_default,
                label="Similarity Method",
            )
            gr.Markdown("### Similarity Hyperparameters")
            similarity_basic = create_hyperparam_widgets(
                max_params=MAX_BASIC_PARAMS,
                initial_specs=_filter_specs(similarity_default, "similarity", advanced=False),
            )
            with gr.Accordion("Similarity Advanced Hyperparameters", open=False):
                similarity_adv = create_hyperparam_widgets(
                    max_params=MAX_ADV_PARAMS,
                    initial_specs=_filter_specs(similarity_default, "similarity", advanced=True),
                )

            run_pre_sim_button = gr.Button("Run Preprocess + Similarity", variant="primary")
            pre_sim_status = gr.Textbox(label="Status", lines=6, interactive=False)
            pre_sim_summary = gr.JSON(label="Preprocess/Similarity Summary")

        with gr.Tab("3) Clustering & Metrics"):
            clustering_method = gr.Dropdown(
                choices=list_method_choices("clustering"),
                value=clustering_default,
                label="Clustering Method",
            )
            gr.Markdown("### Clustering Hyperparameters")
            clustering_basic = create_hyperparam_widgets(
                max_params=MAX_BASIC_PARAMS,
                initial_specs=_filter_specs(clustering_default, "clustering", advanced=False),
            )
            with gr.Accordion("Clustering Advanced Hyperparameters", open=False):
                clustering_adv = create_hyperparam_widgets(
                    max_params=MAX_ADV_PARAMS,
                    initial_specs=_filter_specs(clustering_default, "clustering", advanced=True),
                )
            run_clustering_button = gr.Button("Run Clustering", variant="primary")
            clustering_status = gr.Textbox(label="Status", lines=6, interactive=False)
            clustering_summary = gr.JSON(label="Clustering Summary")

            report_method = gr.Dropdown(
                choices=list_method_choices("report"),
                value=report_default,
                label="Report Method",
            )
            gr.Markdown("### Report Hyperparameters")
            report_basic = create_hyperparam_widgets(
                max_params=MAX_BASIC_PARAMS,
                initial_specs=_filter_specs(report_default, "report", advanced=False),
            )
            with gr.Accordion("Report Advanced Hyperparameters", open=False):
                report_adv = create_hyperparam_widgets(
                    max_params=MAX_ADV_PARAMS,
                    initial_specs=_filter_specs(report_default, "report", advanced=True),
                )
            run_report_button = gr.Button("Run Report", variant="secondary")
            report_status = gr.Textbox(label="Status", lines=6, interactive=False)
            report_json = gr.JSON(label="Report Output")

        with gr.Tab("4) Visualization"):
            viz_method = gr.Dropdown(
                choices=list_method_choices("visualization"),
                value=viz_default,
                label="Visualization Method",
            )
            gr.Markdown("### Visualization Hyperparameters")
            viz_basic = create_hyperparam_widgets(
                max_params=MAX_BASIC_PARAMS,
                initial_specs=_filter_specs(viz_default, "visualization", advanced=False),
            )
            with gr.Accordion("Visualization Advanced Hyperparameters", open=False):
                viz_adv = create_hyperparam_widgets(
                    max_params=MAX_ADV_PARAMS,
                    initial_specs=_filter_specs(viz_default, "visualization", advanced=True),
                )
            run_viz_button = gr.Button("Render Figure", variant="primary")
            viz_status = gr.Textbox(label="Status", lines=6, interactive=False)
            viz_figure = gr.Plot(label="Figure")

        validate_button.click(
            fn=_on_validate,
            inputs=[state, data_file, coord_file],
            outputs=[state, validation_status, validation_summary, snapshot_box],
        )

        preprocess_method.change(
            fn=lambda method: _update_hyper_widgets("preprocess", method),
            inputs=[preprocess_method],
            outputs=[*preprocess_basic.as_outputs(), *preprocess_adv.as_outputs()],
        )
        similarity_method.change(
            fn=lambda method: _update_hyper_widgets("similarity", method),
            inputs=[similarity_method],
            outputs=[*similarity_basic.as_outputs(), *similarity_adv.as_outputs()],
        )

        run_pre_sim_button.click(
            fn=_on_run_preprocess_similarity,
            inputs=[
                state,
                preprocess_method,
                similarity_method,
                *preprocess_basic.as_inputs(),
                *preprocess_adv.as_inputs(),
                *similarity_basic.as_inputs(),
                *similarity_adv.as_inputs(),
            ],
            outputs=[state, pre_sim_status, pre_sim_summary, snapshot_box],
        )

        clustering_method.change(
            fn=lambda method: _update_hyper_widgets("clustering", method),
            inputs=[clustering_method],
            outputs=[*clustering_basic.as_outputs(), *clustering_adv.as_outputs()],
        )
        report_method.change(
            fn=lambda method: _update_hyper_widgets("report", method),
            inputs=[report_method],
            outputs=[*report_basic.as_outputs(), *report_adv.as_outputs()],
        )

        run_clustering_button.click(
            fn=_on_run_clustering,
            inputs=[
                state,
                clustering_method,
                *clustering_basic.as_inputs(),
                *clustering_adv.as_inputs(),
            ],
            outputs=[state, clustering_status, clustering_summary, snapshot_box],
        )
        run_report_button.click(
            fn=_on_run_report,
            inputs=[state, report_method, *report_basic.as_inputs(), *report_adv.as_inputs()],
            outputs=[state, report_status, report_json, snapshot_box],
        )

        viz_method.change(
            fn=lambda method: _update_hyper_widgets("visualization", method),
            inputs=[viz_method],
            outputs=[*viz_basic.as_outputs(), *viz_adv.as_outputs()],
        )

        run_viz_button.click(
            fn=_on_run_visualization,
            inputs=[state, viz_method, *viz_basic.as_inputs(), *viz_adv.as_inputs()],
            outputs=[state, viz_figure, viz_status, snapshot_box],
        )

    demo.queue()
    return demo


if __name__ == "__main__":
    build_demo().launch()
