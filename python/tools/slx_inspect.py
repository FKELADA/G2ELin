"""Reads G2ELib_V1.slx (a Simulink *library*) directly, for P6 (full EMT)
scoping and block-by-block porting work.

A ``.slx`` file is a zip archive of XML — not an opaque binary — so the
actual block diagram (block types, parameters, and wiring) is readable
without MATLAB/Simulink. This is the tool that made that discovery useful
rather than academic: it walks the library's subsystem hierarchy (starting
from ``system_root``, following ``SubSystem`` blocks' nested ``<System
Ref="..."/>`` links) and can print either a summary inventory or one
subsystem's full block-and-wiring detail for hand tracing.

Usage:
    python tools/slx_inspect.py inventory   [path-to-slx]   # library-wide summary
    python tools/slx_inspect.py show <system_id> [path-to-slx]  # one subsystem's blocks + lines

The block-type histogram matters for scoping: G2ELib_V1.slx's 1135 blocks
across 55 subsystems are overwhelmingly ordinary signal-flow blocks (Gain,
Sum, Product, Integrator, TransferFcn, StateSpace — fully readable, exact
parameter values included) rather than opaque Simscape/SimPowerSystems
physical-network blocks (PMComponent/PMIOPort/Ground — only ~67 of the
1135, concentrated in the transformer/load/network subsystems, where the
underlying physics is standard and already known from the EMT port). See
docs/emt_inventory.md for the full catalog this tool generated, and its
worked example (Functions/symGFM_types.m's Droop Q-V control, independently
verified against this library's actual `system_1846` block diagram).
"""

from __future__ import annotations

import sys
import zipfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

DEFAULT_SLX_PATH = Path(__file__).resolve().parents[2] / "G2ELib_V1.slx"


@dataclass
class BlockInfo:
    block_type: str
    name: str
    sid: str
    params: dict[str, str]
    subsystem_ref: str | None  # for BlockType="SubSystem": the nested system_id


@dataclass
class LineInfo:
    src: str  # "<sid>#out:<port>"
    dst: str  # "<sid>#in:<port>"


@dataclass
class SystemInfo:
    system_id: str
    blocks: list[BlockInfo]
    lines: list[LineInfo]


def load_systems(slx_path: Path = DEFAULT_SLX_PATH) -> dict[str, SystemInfo]:
    """Parses every ``simulink/systems/system_*.xml`` entry in the archive."""
    systems: dict[str, SystemInfo] = {}
    with zipfile.ZipFile(slx_path) as z:
        names = [n for n in z.namelist() if n.startswith("simulink/systems/system_") and n.endswith(".xml")]
        for n in names:
            system_id = n.rsplit("/", 1)[-1][: -len(".xml")]
            root = ET.fromstring(z.read(n))
            blocks = []
            for b in root.findall("Block"):
                params = {p.get("Name"): (p.text or "") for p in b.findall("P")}
                sysref = b.find("System")
                blocks.append(
                    BlockInfo(
                        block_type=b.get("BlockType", ""),
                        name=b.get("Name", ""),
                        sid=b.get("SID", ""),
                        params=params,
                        subsystem_ref=sysref.get("Ref") if sysref is not None else None,
                    )
                )
            lines = []
            for line_el in root.findall("Line"):
                _collect_lines(line_el, lines)
            systems[system_id] = SystemInfo(system_id=system_id, blocks=blocks, lines=lines)
    return systems


def _collect_lines(line_el, out: list[LineInfo]) -> None:
    """A <Line> can have nested <Branch> elements (a fan-out signal) — each
    one is its own src/dst edge, sharing the parent's Src.
    """
    src = None
    for p in line_el.findall("P"):
        if p.get("Name") == "Src":
            src = p.text
    for p in line_el.findall("P"):
        if p.get("Name") == "Dst" and src:
            out.append(LineInfo(src=src, dst=p.text))
    for branch in line_el.findall("Branch"):
        _collect_lines_with_src(branch, src, out)


def _collect_lines_with_src(el, src: str | None, out: list[LineInfo]) -> None:
    for p in el.findall("P"):
        if p.get("Name") == "Dst" and src:
            out.append(LineInfo(src=src, dst=p.text))
    for branch in el.findall("Branch"):
        _collect_lines_with_src(branch, src, out)


def block_type_histogram(systems: dict[str, SystemInfo]) -> Counter:
    counts: Counter = Counter()
    for sysinfo in systems.values():
        for b in sysinfo.blocks:
            counts[b.block_type] += 1
    return counts


def hierarchy_lines(systems: dict[str, SystemInfo], root: str = "system_root") -> list[str]:
    out: list[str] = []

    def walk(system_id: str, depth: int, seen: set[str]) -> None:
        if system_id in seen or system_id not in systems:
            return
        seen.add(system_id)
        for b in systems[system_id].blocks:
            if b.block_type == "SubSystem" and b.subsystem_ref:
                n = len(systems.get(b.subsystem_ref, SystemInfo(b.subsystem_ref, [], [])).blocks)
                out.append("  " * depth + f"- {b.name}  [{b.subsystem_ref}, {n} blocks]")
                walk(b.subsystem_ref, depth + 1, seen)

    walk(root, 0, set())
    return out


def print_inventory(slx_path: Path = DEFAULT_SLX_PATH) -> None:
    systems = load_systems(slx_path)
    hist = block_type_histogram(systems)
    print("=== Block type histogram (whole library) ===")
    for btype, count in hist.most_common():
        print(f"  {count:4d}  {btype}")
    print(f"\ntotal blocks: {sum(hist.values())}  across {len(systems)} systems\n")
    print("=== Hierarchy from system_root ===")
    for line in hierarchy_lines(systems):
        print(line)


def print_system(system_id: str, slx_path: Path = DEFAULT_SLX_PATH) -> None:
    systems = load_systems(slx_path)
    if system_id not in systems:
        raise SystemExit(f"no such system: {system_id} (have {sorted(systems)[:5]}...)")
    info = systems[system_id]
    print(f"=== {system_id}: {len(info.blocks)} blocks, {len(info.lines)} lines ===\n")
    for b in info.blocks:
        interesting = {
            k: v for k, v in b.params.items()
            if k in ("Gain", "Inputs", "Operator", "A", "B", "C", "D", "U0", "Y0", "Numerator", "Denominator")
        }
        extra = f"  {interesting}" if interesting else ""
        print(f"  [{b.sid}] {b.block_type:<12} {b.name!r}{extra}")
    print()
    for line in info.lines:
        print(f"  {line.src}  ->  {line.dst}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)
    cmd = args[0]
    if cmd == "inventory":
        path = Path(args[1]) if len(args) > 1 else DEFAULT_SLX_PATH
        print_inventory(path)
    elif cmd == "show":
        if len(args) < 2:
            raise SystemExit("usage: slx_inspect.py show <system_id> [path-to-slx]")
        path = Path(args[2]) if len(args) > 2 else DEFAULT_SLX_PATH
        print_system(args[1], path)
    else:
        raise SystemExit(__doc__)
