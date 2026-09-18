"""FastAPI backend for G2ELin — P7 (web UI) first slice, plus a network
editor / from-scratch builder layer.

Wraps the real ``g2elin_core`` compute functions (no mocking): static power
flow, linear modal analysis, time-series load flow, and nonlinear EMT
time-domain simulation. Two ways to get a ``Network`` to run these on:

- ``/api/presets/{id}/...`` — one of the fixed, named presets in
  ``presets.py``.
- ``/api/network/...`` (``network_routes.py``) — an arbitrary ``Network``
  JSON body, e.g. a preset the frontend cloned and edited, or one built
  from scratch in the drag-and-drop canvas. Nothing here is persisted
  server-side; the client holds the edited network and resends it whole
  each call.

Both sets of endpoints delegate to the same functions in ``analysis.py``
(each takes a ``Network`` and returns the exact response model its
endpoint returns) so the actual analysis logic exists in exactly one
place. This module's own handlers are just "resolve preset_id -> Network
(404 on unknown id), then call into analysis.py."

Also serves the plain-JS frontend in ``web/`` as static files alongside the
JSON API, and the built Sphinx documentation (if present) at ``/manual``
for the frontend's Documentation tab — see ``docs/sphinx/api_and_web.md``'s
"Mount layout" section for why ``/manual`` and not ``/docs`` (FastAPI's own
Swagger UI already lives there), and why mount order matters (specific
routes/mounts before the catch-all ``/``).
"""

from __future__ import annotations

import logging
import os
import threading
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request

from fastapi.responses import StreamingResponse
from fastapi.staticfiles import StaticFiles

from g2elin_core.network.schema import Network

from . import analysis
from .network_routes import router as network_router
from .presets import PRESETS, get_preset
from .schemas import (
    BatchPowerFlowRequest,
    BatchPowerFlowResponse,
    EmtRequest,
    EmtResponse,
    FreeResponseRequest,
    FreeResponseResponse,
    ModalResponse,
    ModeShapeRequest,
    ModeShapeResponse,
    PowerFlowRequest,
    PowerFlowResponse,
    PresetSummary,
    SensitivityRequest,
    SensitivityResponse,
    StatesResponse,
    StepResponseRequest,
    StepResponseResponse,
    TimeSeriesResponse,
    TopologyResponse,
)

logger = logging.getLogger("g2elin_api")

app = FastAPI(title="G2ELin API", description="Power-system analysis over G2ELin's ported presets.")


@app.middleware("http")
async def _revalidate_static(request: Request, call_next):
    """The web UI's HTML/CSS/JS are served without a cache lifetime, so
    browsers may reuse a stale copy after an update. ``no-cache`` makes them
    revalidate each load (a cheap 304 via StaticFiles' ETag when unchanged).
    """
    response = await call_next(request)
    if not request.url.path.startswith("/api/"):
        response.headers.setdefault("Cache-Control", "no-cache")
    return response


# Presets whose models are pre-built at startup when G2ELIN_WARMUP=1 (set in
# the Docker image): between them they cover every unit type (SM, GFM, GFL,
# infinite bus), so the one-off symbolic derivation per type happens before
# the first visitor asks for a modal analysis or an EMT run, not during it.
_WARMUP_PRESETS = ("wscc9_3sm", "wscc9_1sm_1gfm_1gfl", "sm_smib")


@app.on_event("startup")
def _warm_up_models() -> None:
    if os.environ.get("G2ELIN_WARMUP") != "1":
        return

    def run() -> None:
        for preset_id in _WARMUP_PRESETS:
            try:
                net = get_preset(preset_id).build()
                analysis.build_modal_from_network(net)
                analysis.build_nonlinear_model_from_network(net)
                logger.info("warm-up: built models for %s", preset_id)
            except Exception:  # noqa: BLE001 -- best effort; a real request reports the error properly
                logger.exception("warm-up failed for %s", preset_id)

    threading.Thread(target=run, name="g2elin-warmup", daemon=True).start()


def _resolve_preset_network(preset_id: str) -> Network:
    try:
        return get_preset(preset_id).build()
    except KeyError as e:
        raise HTTPException(status_code=404, detail=str(e)) from e


@app.get("/api/presets", response_model=list[PresetSummary])
def list_presets() -> list[PresetSummary]:
    summaries = []
    for info in sorted(PRESETS.values(), key=lambda p: p.id):
        net = info.build()
        summaries.append(
            PresetSummary(
                id=info.id,
                name=info.name,
                description=info.description,
                n_buses=len(net.buses),
                n_lines=len(net.lines),
                n_transformers=len(net.transformers),
                n_loads=len(net.loads),
                n_der_units=len(net.der_units),
                der_unit_types=sorted({d.unit_type.value for d in net.der_units}),
            )
        )
    return summaries


@app.get("/api/presets/{preset_id}/network", response_model=Network)
def get_preset_network(preset_id: str) -> Network:
    """The preset's own ``Network`` definition, verbatim -- how the web UI
    clones a preset into an editable copy (see ``/api/network/*`` below).
    """
    return _resolve_preset_network(preset_id)


@app.get("/api/presets/{preset_id}/topology", response_model=TopologyResponse)
def get_topology(preset_id: str) -> TopologyResponse:
    return analysis.topology_response(_resolve_preset_network(preset_id))


@app.get("/api/powerflow/algorithms")
def list_powerflow_algorithms() -> dict[str, str]:
    """Solver ids accepted by ``PowerFlowOptions.algorithm`` -> display label."""
    return analysis.POWERFLOW_ALGORITHMS


@app.post("/api/presets/{preset_id}/powerflow", response_model=PowerFlowResponse)
def run_powerflow(preset_id: str, req: PowerFlowRequest | None = None) -> PowerFlowResponse:
    # The body is optional: a bare POST solves with pandapower's defaults.
    return analysis.powerflow_response(_resolve_preset_network(preset_id), req.options if req else None)


@app.post("/api/presets/{preset_id}/powerflow/batch", response_model=BatchPowerFlowResponse)
def run_powerflow_batch(preset_id: str, req: BatchPowerFlowRequest) -> BatchPowerFlowResponse:
    return analysis.batch_powerflow_response(_resolve_preset_network(preset_id), req)


@app.post("/api/presets/{preset_id}/modal", response_model=ModalResponse)
def run_modal(preset_id: str) -> ModalResponse:
    return analysis.modal_response(_resolve_preset_network(preset_id))


@app.post("/api/presets/{preset_id}/modal/sensitivity", response_model=SensitivityResponse)
def run_modal_sensitivity(preset_id: str, req: SensitivityRequest) -> SensitivityResponse:
    return analysis.modal_sensitivity_response(_resolve_preset_network(preset_id), req)


@app.post("/api/presets/{preset_id}/modal/mode_shape", response_model=ModeShapeResponse)
def run_modal_mode_shape(preset_id: str, req: ModeShapeRequest) -> ModeShapeResponse:
    return analysis.modal_mode_shape_response(_resolve_preset_network(preset_id), req)


@app.post("/api/presets/{preset_id}/modal/free_response", response_model=FreeResponseResponse)
def run_modal_free_response(preset_id: str, req: FreeResponseRequest) -> FreeResponseResponse:
    return analysis.modal_free_response_response(_resolve_preset_network(preset_id), req)


@app.post("/api/presets/{preset_id}/modal/step_response", response_model=StepResponseResponse)
def run_modal_step_response(preset_id: str, req: StepResponseRequest) -> StepResponseResponse:
    return analysis.modal_step_response_response(_resolve_preset_network(preset_id), req)


@app.post("/api/presets/{preset_id}/timeseries", response_model=TimeSeriesResponse)
def run_timeseries(preset_id: str) -> TimeSeriesResponse:
    return analysis.timeseries_response(_resolve_preset_network(preset_id))


@app.get("/api/presets/{preset_id}/states", response_model=StatesResponse)
def list_states(preset_id: str) -> StatesResponse:
    """State names for the nonlinear (EMT) model — lets the frontend
    populate its state pickers without hard-coding per-preset names (which
    differ: e.g. only synchronous machines have ``dw_r_*`` states, and
    which DER id is the slack varies by preset).
    """
    return analysis.states_response(_resolve_preset_network(preset_id))


@app.post("/api/presets/{preset_id}/emt", response_model=EmtResponse)
def run_emt(preset_id: str, req: EmtRequest) -> EmtResponse:
    return analysis.emt_response(_resolve_preset_network(preset_id), req)


@app.post("/api/presets/{preset_id}/emt/live")
def run_emt_live(preset_id: str, req: EmtRequest, request: Request) -> StreamingResponse:
    """Newline-delimited JSON, one line per solver step -- see
    ``analysis.emt_live_stream``'s own docstring for the exact line
    shapes. Validation (``prepare_emt_live``) runs here, synchronously,
    *before* the streaming response is constructed, so a bad request
    still gets a normal 422 instead of a broken stream.
    """
    plan = analysis.prepare_emt_live(_resolve_preset_network(preset_id), req)
    return StreamingResponse(
        analysis.emt_live_stream(plan, req.perturb_kind, request), media_type="application/x-ndjson"
    )


app.include_router(network_router)

_repo_root = Path(__file__).resolve().parents[2]
_docs_dir = _repo_root / "docs" / "sphinx" / "_build" / "html"
if _docs_dir.is_dir():
    app.mount("/manual", StaticFiles(directory=str(_docs_dir), html=True), name="manual")

_web_dir = _repo_root / "web"
if _web_dir.is_dir():
    app.mount("/", StaticFiles(directory=str(_web_dir), html=True), name="web")
