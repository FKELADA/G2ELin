"""Regenerates the IEEE preset data files from pandapower's own test cases.

Run from ``python/``::

    python tools/build_ieee_presets.py

The conversion (``network.pandapower_import``) is not cheap: each generator's
voltage setpoint has to be retuned until the *grid* bus it feeds sits where the
published case put it, which takes tens of power-flow solves -- 6 s for the
14-bus case and 40 s for the 118-bus one. That is far too slow to do every time
a preset is asked for, so it is done once, here, and the result is written to
``src/g2elin_core/network/data/*.json`` for the preset functions to load.

Re-run it after changing the importer, and check the report it prints: it says
what the source case did not contain and the importer had to assume.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "src" / "g2elin_core" / "network" / "data"

# (file stem, display name, pandapower factory name, nominal frequency)
CASES = [
    ("ieee14", "IEEE_14_bus", "case14", 60.0),
    ("ieee39", "IEEE_39_bus_New_England", "case39", 60.0),
    ("ieee118", "IEEE_118_bus", "case118", 60.0),
]


def main() -> int:
    import pandapower as pp
    import pandapower.networks as nw

    from g2elin_core.network.pandapower_import import from_pandapower, regulate_grid_buses
    from g2elin_core.network.validation import validate_network
    from g2elin_core.powerflow import run_power_flow

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    failed = False
    for stem, name, factory, f_hz in CASES:
        source = getattr(nw, factory)()
        pp.runpp(source)
        network, report = from_pandapower(source, name=name, f_hz=f_hz)
        regulate_grid_buses(network, report, max_iter=60)

        result = run_power_flow(network)
        errors = [i.message for i in validate_network(network) if i.severity == "error"]
        worst_dv = max(
            abs(float(result.bus_table().set_index("bus")["vm_pu"][b]) - float(source.res_bus.vm_pu[b]))
            for b in source.bus.index
        ) if result.converged else float("nan")

        (DATA_DIR / f"{stem}.json").write_text(network.model_dump_json(indent=1), encoding="utf-8")
        size_kb = (DATA_DIR / f"{stem}.json").stat().st_size / 1024
        print(f"{stem:8s} {len(network.buses):4d} buses {len(network.der_units):3d} units "
              f"-> {size_kb:6.0f} kB | converged={result.converged} worst dV vs source {worst_dv:.2e} pu")
        print(f"         {report.summary()}")
        if errors:
            print(f"         VALIDATION ERRORS: {errors[:2]}")
        if errors or not result.converged:
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
