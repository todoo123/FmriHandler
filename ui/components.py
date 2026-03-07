from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

import gradio as gr

from .registry import HyperParamSpec


@dataclass(slots=True)
class HyperParamWidgets:
    max_params: int
    keys: list[gr.Textbox]
    texts: list[gr.Textbox]
    numbers: list[gr.Number]
    bools: list[gr.Checkbox]
    choices: list[gr.Dropdown]

    def as_inputs(self) -> list[gr.Component]:
        return [*self.keys, *self.texts, *self.numbers, *self.bools, *self.choices]

    def as_outputs(self) -> list[gr.Component]:
        return self.as_inputs()


def create_hyperparam_widgets(max_params: int = 8) -> HyperParamWidgets:
    keys: list[gr.Textbox] = []
    texts: list[gr.Textbox] = []
    numbers: list[gr.Number] = []
    bools: list[gr.Checkbox] = []
    choices: list[gr.Dropdown] = []

    for idx in range(max_params):
        keys.append(
            gr.Textbox(
                visible=False,
                value="",
                label=f"key_{idx}",
                interactive=False,
            )
        )
        texts.append(gr.Textbox(visible=False, label=f"text_{idx}", value=""))
        numbers.append(gr.Number(visible=False, label=f"num_{idx}", value=None))
        bools.append(gr.Checkbox(visible=False, label=f"bool_{idx}", value=False))
        choices.append(
            gr.Dropdown(
                visible=False,
                label=f"choice_{idx}",
                choices=[],
                value=None,
                allow_custom_value=False,
            )
        )

    return HyperParamWidgets(
        max_params=max_params,
        keys=keys,
        texts=texts,
        numbers=numbers,
        bools=bools,
        choices=choices,
    )


def _resolve_input_kind(spec: HyperParamSpec) -> str:
    if spec.input_kind != "auto":
        return spec.input_kind
    if spec.choices:
        return "dropdown"
    if spec.value_type == "bool":
        return "checkbox"
    if spec.value_type in ("int", "float"):
        return "number"
    return "text"


def _normalize_text_default(spec: HyperParamSpec) -> str:
    if spec.default is None:
        return ""
    if spec.value_type == "json":
        return json.dumps(spec.default)
    return str(spec.default)


def _spec_label(spec: HyperParamSpec) -> str:
    if spec.required:
        return f"{spec.label} *"
    return spec.label


def build_widget_updates(specs: list[HyperParamSpec], max_params: int) -> list[dict[str, Any]]:
    """Build gr.update payloads for key/text/number/bool/choice slots."""
    updates: list[dict[str, Any]] = []
    visible_specs = specs[:max_params]

    for idx in range(max_params):
        if idx >= len(visible_specs):
            updates.extend(
                [
                    gr.update(value="", visible=False),
                    gr.update(value="", visible=False, label=f"text_{idx}", info=None),
                    gr.update(value=None, visible=False, label=f"num_{idx}", info=None),
                    gr.update(value=False, visible=False, label=f"bool_{idx}", info=None),
                    gr.update(
                        value=None,
                        choices=[],
                        visible=False,
                        label=f"choice_{idx}",
                        info=None,
                    ),
                ]
            )
            continue

        spec = visible_specs[idx]
        kind = _resolve_input_kind(spec)
        info = spec.help_text if spec.help_text else None
        label = _spec_label(spec)

        updates.append(gr.update(value=spec.name, visible=False))
        updates.append(
            gr.update(
                value=_normalize_text_default(spec),
                visible=kind == "text",
                label=label,
                info=info,
            )
        )
        updates.append(
            gr.update(
                value=spec.default,
                visible=kind == "number",
                label=label,
                info=info,
            )
        )
        updates.append(
            gr.update(
                value=bool(spec.default) if spec.default is not None else False,
                visible=kind == "checkbox",
                label=label,
                info=info,
            )
        )
        updates.append(
            gr.update(
                value=spec.default,
                choices=spec.choices,
                visible=kind == "dropdown",
                label=label,
                info=info,
            )
        )

    return updates


def split_widget_values(
    flat_values: list[Any],
    max_params: int,
) -> tuple[list[str], list[Any], list[Any], list[Any], list[Any]]:
    n = max_params
    keys = [str(v or "") for v in flat_values[0:n]]
    texts = flat_values[n : 2 * n]
    numbers = flat_values[2 * n : 3 * n]
    bools = flat_values[3 * n : 4 * n]
    choices = flat_values[4 * n : 5 * n]
    return keys, texts, numbers, bools, choices


def parse_hyperparam_values(
    specs: list[HyperParamSpec],
    key_values: list[str],
    text_values: list[Any],
    number_values: list[Any],
    bool_values: list[Any],
    choice_values: list[Any],
) -> dict[str, Any]:
    parsed: dict[str, Any] = {}
    spec_by_name = {spec.name: spec for spec in specs}
    # Always seed defaults first, so missing UI inputs remain safe.
    for spec in specs:
        if spec.default is not None:
            parsed[spec.name] = _cast_value(spec, spec.default)

    # Prefer explicit hidden key mapping, but gracefully fall back to
    # positional slot mapping when hidden key values are not propagated
    # by the UI runtime.
    for idx in range(max(len(key_values), len(specs))):
        key = key_values[idx] if idx < len(key_values) else ""
        has_explicit_key = bool(key)
        spec = spec_by_name.get(key) if key else None
        if spec is None and idx < len(specs):
            spec = specs[idx]
        if spec is None:
            continue

        kind = _resolve_input_kind(spec)
        if kind == "text":
            raw = text_values[idx]
        elif kind == "number":
            raw = number_values[idx]
        elif kind == "checkbox":
            raw = bool_values[idx]
        elif kind == "dropdown":
            raw = choice_values[idx]
        else:
            raw = text_values[idx]

        # If no explicit key is attached and this is a bool field with
        # default=True, keep default unless True was explicitly sent.
        if (
            kind == "checkbox"
            and not has_explicit_key
            and bool(raw) is False
            and spec.default is True
        ):
            continue

        value = _cast_value(spec, raw)
        if value is None and not spec.required:
            continue
        if value is None and spec.required:
            raise ValueError(f"Required hyperparameter '{spec.name}' is missing.")
        parsed[spec.name] = value

    # Final required check after default + user merge.
    for spec in specs:
        if spec.required and spec.name not in parsed:
            raise ValueError(f"Required hyperparameter '{spec.name}' is missing.")

    return parsed


def _cast_value(spec: HyperParamSpec, raw: Any) -> Any:
    if raw is None or raw == "":
        if spec.default is not None:
            return spec.default
        return None

    if spec.value_type == "str":
        return str(raw)
    if spec.value_type == "int":
        return int(raw)
    if spec.value_type == "float":
        return float(raw)
    if spec.value_type == "bool":
        return bool(raw)
    if spec.value_type == "csv_int":
        if isinstance(raw, (list, tuple)):
            return tuple(int(v) for v in raw)
        return tuple(int(chunk.strip()) for chunk in str(raw).split(",") if chunk.strip())
    if spec.value_type == "json":
        if isinstance(raw, (dict, list)):
            return raw
        return json.loads(str(raw))
    return raw
