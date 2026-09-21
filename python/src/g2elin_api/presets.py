"""The set of networks the API exposes.

Mirrors ``g2elin_core.network.presets`` — the API layer doesn't add any new
networks, it just exposes them over HTTP. A network editor (accepting
arbitrary user-built ``Network`` instances) is added alongside this
(``/api/network/*``, see ``network_routes.py``); this module still owns the
fixed, named preset catalog those endpoints' "start from a preset" flow
clones from.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from g2elin_core.network.presets import (
    cigre_interconnected_1sm_1gfm_1gfl,
    cigre_islanded_1sm_1gfm_1gfl,
    cigre_islanded_1sm_2gfm_1gfl,
    cigre_islanded_1sm_3gfm_1gfl,
    cigre_islanded_2sm_2gfm_2gfl,
    gfl_smib,
    gfl_smsm,
    gfm_smib,
    gfm_smsm,
    sm_smib,
    sm_smsm,
    wscc9_1gfm_2gfl,
    wscc9_1sm_1gfm_1gfl,
    wscc9_1sm_2gfl,
    wscc9_1sm_2gfm,
    wscc9_2sm_1gfl,
    wscc9_2sm_1gfm,
    wscc9_3gfm,
    kundur_two_area,
    kundur_two_area_classic,
    wscc9_3sm,
)
from g2elin_core.network.schema import Network


@dataclass(frozen=True)
class PresetInfo:
    id: str
    name: str
    description: str
    build: Callable[[], Network]


# (id, name, description, build) -- one row per registered preset. Grouped
# by family (WSCC / CIGRE islanded / SMIB) to match the web UI's <optgroup>
# grouping in loadPresets().
_PRESET_SPECS: list[tuple[str, str, str, Callable[[], Network]]] = [
    ("wscc9_3sm", "WSCC 9-bus (3 SM)",
     "Classic WSCC 9-bus system, 3 synchronous machines. From WSCC_raw.m + script_WSCC.m.",
     wscc9_3sm),
    ("wscc9_2sm_1gfl", "WSCC 9-bus (2 SM + 1 GFL)",
     "WSCC 9-bus, 2 synchronous machines + 1 grid-following converter. script_WSCC.m case 'WSCC_2SM_1GFL'.",
     wscc9_2sm_1gfl),
    ("wscc9_1sm_2gfl", "WSCC 9-bus (1 SM + 2 GFL)",
     "WSCC 9-bus, 1 synchronous machine + 2 grid-following converters. script_WSCC.m case 'WSCC_1SM_2GFL'.",
     wscc9_1sm_2gfl),
    ("wscc9_1sm_1gfm_1gfl", "WSCC 9-bus (1 SM + 1 GFM + 1 GFL)",
     "WSCC 9-bus, 1 synchronous machine + 1 grid-forming (Droop) + 1 grid-following converter. "
     "script_WSCC.m case 'WSCC_1SM_1GFM_1GFL'.",
     wscc9_1sm_1gfm_1gfl),
    ("wscc9_1sm_2gfm", "WSCC 9-bus (1 SM + 2 GFM)",
     "WSCC 9-bus, 1 synchronous machine + 2 grid-forming (Droop) converters. script_WSCC.m case 'WSCC_1SM_2GFM'.",
     wscc9_1sm_2gfm),
    ("wscc9_2sm_1gfm", "WSCC 9-bus (2 SM + 1 GFM)",
     "WSCC 9-bus, 2 synchronous machines + 1 grid-forming (Droop) converter. script_WSCC.m case 'WSCC_2SM_1GFM'.",
     wscc9_2sm_1gfm),
    ("wscc9_1gfm_2gfl", "WSCC 9-bus (1 GFM + 2 GFL)",
     "WSCC 9-bus, 1 grid-forming (Droop, slack) + 2 grid-following converters. Power flow only -- a GFM "
     "slack isn't wired up for modal/EMT yet (see the preset's own docstring).",
     wscc9_1gfm_2gfl),
    ("wscc9_3gfm", "WSCC 9-bus (3 GFM)",
     "WSCC 9-bus, 3 grid-forming (Droop) converters, one as slack. Power flow only -- see wscc9_1gfm_2gfl.",
     wscc9_3gfm),
    ("cigre_islanded_1sm_1gfm_1gfl", "CIGRE islanded (1 SM + 1 GFM + 1 GFL)",
     "CIGRE MV benchmark feeder, islanded, 1 SM + 1 GFM + 1 GFL. "
     "CIGRE_raw.m + CIGRE_Islanded.m case 'CIGRE_Islanded_1SM_1GFM_1GFL'.",
     cigre_islanded_1sm_1gfm_1gfl),
    ("cigre_islanded", "CIGRE islanded (1 SM + 2 GFM + 1 GFL)",
     "CIGRE MV benchmark feeder, islanded, mixed SM/GFM/GFL. "
     "From CIGRE_raw.m + CIGRE/script_CIGRE_Islanded.m.",
     cigre_islanded_1sm_2gfm_1gfl),
    ("cigre_islanded_1sm_3gfm_1gfl", "CIGRE islanded (1 SM + 3 GFM + 1 GFL)",
     "CIGRE MV benchmark feeder, islanded, 1 SM + 3 GFM + 1 GFL. "
     "CIGRE_raw.m + CIGRE_Islanded.m case 'CIGRE_Islanded_1SM_3GFM_1GFL'.",
     cigre_islanded_1sm_3gfm_1gfl),
    ("cigre_islanded_2sm_2gfm_2gfl", "CIGRE islanded (2 SM + 2 GFM + 2 GFL)",
     "CIGRE MV benchmark feeder, islanded, 2 SM + 2 GFM + 2 GFL. "
     "CIGRE_raw.m + CIGRE_Islanded.m case 'CIGRE_Islanded_2SM_2GFM_2GFL'.",
     cigre_islanded_2sm_2gfm_2gfl),
    ("cigre_interconnected_1sm_1gfm_1gfl", "CIGRE interconnected (1 SM + 1 GFM + 1 GFL)",
     "CIGRE MV benchmark feeder connected to the upstream grid (an infinite bus behind a transformer "
     "at node 1, identical to the unit transformers), 1 SM + 1 GFM + 1 GFL. Topology/loads from CIGRE_raw.m (S0=1); "
     "grid attachment reconstructed -- preset_networks.m names the case but never defines it.",
     cigre_interconnected_1sm_1gfm_1gfl),
    ("sm_smib", "SMIB (SM)",
     "Single synchronous machine against an infinite bus. Topology/line/load from "
     "SMIB_raw.m; DER attachment reconstructed (the driving MATLAB script was never "
     "finished) -- see network.presets.sm_smib's docstring.",
     sm_smib),
    ("gfm_smib", "SMIB (GFM, Droop)",
     "Single grid-forming converter (Droop) against an infinite bus. See sm_smib for provenance.",
     gfm_smib),
    ("gfl_smib", "SMIB (GFL)",
     "Single grid-following converter against an infinite bus. See sm_smib for provenance.",
     gfl_smib),
    ("sm_smsm", "SMSM (SM vs SM)",
     "Single synchronous machine against a synchronous machine acting as the grid (slack), on the "
     "SMIB line and load. Reconstruction, like the SMIB presets.",
     sm_smsm),
    ("gfm_smsm", "SMSM (GFM vs SM)",
     "Single grid-forming converter (Droop) against a synchronous machine acting as the grid. See sm_smsm.",
     gfm_smsm),
    ("gfl_smsm", "SMSM (GFL vs SM)",
     "Single grid-following converter against a synchronous machine acting as the grid. The GFL absorbs "
     "0.9 MVAr of the line charging -- at Q = 0 the under-excited grid machine is unstable (see gfl_smsm).",
     gfl_smsm),
    ("kundur_two_area", "Kundur two-area (11-bus, 4 SM)",
     "Kundur's two-area system (Power System Stability and Control, Example 12.6): two areas of two "
     "900 MVA machines, a weak 220 km double-circuit tie carrying 400 MW, and the textbook case for "
     "inter-area oscillations. With this tool's own controls (AVR, PSS and governor).",
     kundur_two_area),
    ("kundur_two_area_classic", "Kundur two-area (book's controls)",
     "The same system under the book's own assumptions -- fast static exciter, no PSS, constant "
     "mechanical torque -- where the inter-area mode near 0.6 Hz comes out negatively damped. The case "
     "the example exists for.",
     kundur_two_area_classic),
]

PRESETS: dict[str, PresetInfo] = {
    id_: PresetInfo(id=id_, name=name, description=description, build=build)
    for id_, name, description, build in _PRESET_SPECS
}


def get_preset(preset_id: str) -> PresetInfo:
    try:
        return PRESETS[preset_id]
    except KeyError:
        raise KeyError(f"unknown preset {preset_id!r}; known presets: {sorted(PRESETS)}") from None
