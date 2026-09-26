# Recording the web-interface tour

A scripted, captioned walk through the whole interface, recorded to video.
It drives real Chrome through Playwright — the tool answers for itself, so
what the video shows is what the tool does.

The result lives in `Documentation/G2ELin_web_UI_tour.mp4` (outside the
repository: it is tens of megabytes).

## What it needs

```powershell
npm install playwright-core            # drives the browser; Chrome itself is reused
pip install imageio-ffmpeg             # ships the ffmpeg binary the encoder uses
```

Neither goes in the project's dependencies: this is a thing you run once in
a while, not part of the tool.

## Running it

Start the interface on a known port, then record:

```powershell
.venv\Scripts\python.exe -m uvicorn g2elin_api.main:app --port 8765
node tools/ui_tour/record_tour.js          # ~25 minutes, writes recording/*.webm
.venv\Scripts\python.exe tools/ui_tour/encode.py   # -> G2ELin_web_ui_tour.mp4
```

Settings, all with working defaults: `TOUR_BASE` (default
`http://127.0.0.1:8765/`), `TOUR_CHROME` (where Chrome is installed),
`TOUR_OUT` (where the raw recording goes), `TOUR_FFMPEG`.

## How it is built

**Every scene is independent.** One that fails is logged and skipped, so a
recording always completes and the log says what was missed. That matters
because a run is half an hour: a tour that aborts at minute twenty has cost
twenty minutes and produced nothing.

**The captions are an overlay** (`overlay.js`), injected into the page
rather than added afterwards — a caption and the thing it describes are on
screen at the same moment, and the cursor is drawn too, so clicks read as
clicks. A caption is cleared when its scene ends; left standing it outlives
what it describes.

**Live tracing stays on for the time-domain runs, and is sped up
afterwards.** The pace of a live trace is the solver's, and nothing in the
browser can hurry it — so `record_tour.js` timestamps each live run into
`segments.json` and `encode.py` speeds exactly those stretches up (5x) when
it writes the MP4. Everything else plays at its real speed. That keeps the
trace filling in, which is the point of showing it, without keeping its
real duration: a 25-minute recording becomes about 13 minutes.

**What is out of service is stated, not assumed.** The breaker hit areas on
the diagram overlap, so a click meant for one can toggle its neighbour —
that is how a load once stayed open through a whole chapter, and the modal
scenes ran on 126 states instead of 128. The clicks stay, because they are
what the viewer sees; `setService()` then says outright what should be
open, so a stray toggle cannot follow the tour into the next scene.

## Changing it

Selectors drift as the interface changes. Before re-recording after a gap,
check the ones the script drives — there is a probe pattern in the git
history for this, and it is much cheaper than finding out twenty minutes
into a take. The log line for each scene gives its duration; a scene that
reports `0.0s` did nothing.
