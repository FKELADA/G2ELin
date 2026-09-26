"""Turn the recorded WebM into an MP4, speeding up the live-trace stretches.

Playwright records in real time. A live time-domain run is paced by the
solver, so it cannot be hurried in the browser -- but tour4.js timestamps
each one into segments.json, so exactly those stretches can be sped up here
and nothing else is touched. The result keeps the live trace (which is the
point: you see it fill) without keeping its real duration.
"""
from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).parent
def _ffmpeg() -> str:
    """ffmpeg, from wherever this machine has one.

    imageio-ffmpeg ships a binary and is a one-line install, which is why it
    is tried first -- it needs no system package and no PATH entry.
    """
    if os.environ.get("TOUR_FFMPEG"):
        return os.environ["TOUR_FFMPEG"]
    try:
        import imageio_ffmpeg

        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        pass
    found = shutil.which("ffmpeg")
    if not found:
        sys.exit("no ffmpeg: pip install imageio-ffmpeg, or set TOUR_FFMPEG")
    return found


FFMPEG = _ffmpeg()


def source() -> pathlib.Path:
    out_dir = pathlib.Path(os.environ.get("TOUR_OUT", HERE / "recording"))
    vids = sorted(out_dir.glob("*.webm"), key=lambda p: p.stat().st_size)
    if not vids:
        sys.exit(f"no recording in {out_dir} -- run record_tour.js first")
    return vids[-1]


def duration(path: pathlib.Path) -> float:
    out = subprocess.run(
        [str(FFMPEG), "-i", str(path)], capture_output=True, text=True, errors="replace"
    ).stderr
    for line in out.splitlines():
        if "Duration:" in line:
            h, m, s = line.split("Duration:")[1].split(",")[0].strip().split(":")
            return int(h) * 3600 + int(m) * 60 + float(s)
    sys.exit("could not read the recording's duration")


def plan(total: float, fast: list[dict], pad: float = 0.6) -> list[tuple[float, float, float]]:
    """(start, end, speed) covering the whole timeline exactly once.

    The fast stretches are trimmed slightly at each end so the caption that
    introduces a run, and the first moment of its result, stay at full speed
    -- what is sped up is the middle, where the trace is filling.
    """
    spans = []
    for f in sorted(fast, key=lambda f: f["start"]):
        a, b = f["start"] + pad, f["end"] - pad
        if b - a > 1.5:
            spans.append((a, b, float(f["speed"])))
    out, cursor = [], 0.0
    for a, b, speed in spans:
        if a > cursor:
            out.append((cursor, a, 1.0))
        out.append((a, b, speed))
        cursor = b
    if cursor < total:
        out.append((cursor, total, 1.0))
    return [(a, b, s) for a, b, s in out if b - a > 0.05]


def main() -> None:
    src = source()
    seg_file = HERE / "segments.json"
    total = duration(src)
    fast = json.loads(seg_file.read_text())["fast"] if seg_file.exists() else []
    parts = plan(total, fast)

    sped = sum(b - a for a, b, s in parts if s != 1.0)
    kept = sum((b - a) / s for a, b, s in parts)
    print(f"recording {total:.1f}s, {len(fast)} live stretches ({sped:.1f}s) -> {kept:.1f}s")

    # One trim+setpts per part, then concat. Video only; there is no audio.
    chains, labels = [], []
    for i, (a, b, s) in enumerate(parts):
        chains.append(f"[0:v]trim=start={a:.3f}:end={b:.3f},setpts=(PTS-STARTPTS)/{s}[v{i}]")
        labels.append(f"[v{i}]")
    graph = ";".join(chains) + ";" + "".join(labels) + f"concat=n={len(parts)}:v=1:a=0[out]"

    out = HERE / "G2ELin_web_ui_tour.mp4"
    cmd = [
        str(FFMPEG), "-y", "-i", str(src),
        "-filter_complex", graph, "-map", "[out]",
        "-c:v", "libx264", "-preset", "slow", "-crf", "23",
        "-pix_fmt", "yuv420p", "-movflags", "+faststart", "-r", "25",
        str(out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    if r.returncode != 0:
        print(r.stderr[-3000:])
        sys.exit("ffmpeg failed")
    print(f"wrote {out}  ({out.stat().st_size / 1e6:.1f} MB, {duration(out):.1f}s)")


if __name__ == "__main__":
    main()
