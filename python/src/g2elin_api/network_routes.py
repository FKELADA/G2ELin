"""``/api/network/*`` — run any analysis on an arbitrary, client-provided
``Network`` instead of one of the fixed presets in ``presets.py``.

Every handler here is a thin wrapper delegating to the same functions in
``analysis.py`` the preset-id-keyed endpoints in ``main.py`` use — the
actual analysis logic lives in exactly one place; this module only differs
in *where the Network comes from* (a POSTed body, not a preset lookup).

Nothing here is persisted server-side: the client (the web UI's network
editor / drag-and-drop builder) holds the ``Network`` and resends it whole
on every call. A ``Network`` body can't ride a GET request, so every
endpoint here is POST, including the two (``/topology``, ``/states``) whose
preset-side counterparts are GET — an intentional, documented divergence,
not an oversight.

pydantic validates the incoming ``network`` field against the exact same
``Network`` model ``g2elin_core`` uses (bus-reference integrity,
exactly-one-slack, etc.) before any handler here runs, via FastAPI's
default request-body validation — an invalid network never reaches
``analysis.py`` at all, it 422s first with pydantic's own per-field error
list (see ``web/index.html``'s ``formatDetail()`` for how the frontend
renders that list as one readable message).
"""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import StreamingResponse

from . import analysis
from .schemas import (
    EmtResponse,
    FreeResponseResponse,
    ModalResponse,
    ModeShapeResponse,
    NetworkEmtRequest,
    NetworkFreeResponseRequest,
    NetworkModeShapeRequest,
    NetworkRequest,
    NetworkSensitivityRequest,
    NetworkStepResponseRequest,
    PowerFlowResponse,
    SensitivityResponse,
    StatesResponse,
    StepResponseResponse,
    TimeSeriesResponse,
    TopologyResponse,
    ValidateResponse,
)

router = APIRouter(prefix="/api/network", tags=["network"])


@router.post("/validate", response_model=ValidateResponse)
def network_validate(req: NetworkRequest) -> ValidateResponse:
    """Pre-flight structural check -- every issue found, not just the
    first, and cheap (no power flow run). Meant to be called before the
    other endpoints, e.g. by the web UI's network editor, so a hand-built
    network's problems can be shown all at once.
    """
    return analysis.validate_network_response(req.network)


@router.post("/topology", response_model=TopologyResponse)
def network_topology(req: NetworkRequest) -> TopologyResponse:
    return analysis.topology_response(req.network)


@router.post("/powerflow", response_model=PowerFlowResponse)
def network_powerflow(req: NetworkRequest) -> PowerFlowResponse:
    return analysis.powerflow_response(req.network)


@router.post("/modal", response_model=ModalResponse)
def network_modal(req: NetworkRequest) -> ModalResponse:
    return analysis.modal_response(req.network)


@router.post("/modal/sensitivity", response_model=SensitivityResponse)
def network_modal_sensitivity(req: NetworkSensitivityRequest) -> SensitivityResponse:
    return analysis.modal_sensitivity_response(req.network, req)


@router.post("/modal/mode_shape", response_model=ModeShapeResponse)
def network_modal_mode_shape(req: NetworkModeShapeRequest) -> ModeShapeResponse:
    return analysis.modal_mode_shape_response(req.network, req)


@router.post("/modal/free_response", response_model=FreeResponseResponse)
def network_modal_free_response(req: NetworkFreeResponseRequest) -> FreeResponseResponse:
    return analysis.modal_free_response_response(req.network, req)


@router.post("/modal/step_response", response_model=StepResponseResponse)
def network_modal_step_response(req: NetworkStepResponseRequest) -> StepResponseResponse:
    return analysis.modal_step_response_response(req.network, req)


@router.post("/timeseries", response_model=TimeSeriesResponse)
def network_timeseries(req: NetworkRequest) -> TimeSeriesResponse:
    return analysis.timeseries_response(req.network)


@router.post("/states", response_model=StatesResponse)
def network_states(req: NetworkRequest) -> StatesResponse:
    return analysis.states_response(req.network)


@router.post("/emt", response_model=EmtResponse)
def network_emt(req: NetworkEmtRequest) -> EmtResponse:
    return analysis.emt_response(req.network, req)


@router.post("/emt/live")
def network_emt_live(req: NetworkEmtRequest, request: Request) -> StreamingResponse:
    """Newline-delimited JSON, one line per solver step -- see
    ``main.py``'s identical preset-side route and
    ``analysis.emt_live_stream``'s own docstring for the line shapes.
    """
    plan = analysis.prepare_emt_live(req.network, req)
    return StreamingResponse(
        analysis.emt_live_stream(plan, req.perturb_kind, request), media_type="application/x-ndjson"
    )
