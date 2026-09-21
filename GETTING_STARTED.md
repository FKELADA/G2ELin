# Getting started with G2ELin

Three steps, no coding, about ten minutes — most of which is waiting for the
download. Nothing is installed outside the folder you unzip.

## 1. Install Python (once per computer)

Download it from <https://www.python.org/downloads/> and run the installer.

> **Windows:** tick **"Add python.exe to PATH"** on the first screen of the
> installer, before clicking Install. That one checkbox is what makes the
> rest work.

> **macOS:** after installing, double-click **"Install
> Certificates.command"** in the Python folder that Finder opens.

## 2. Download G2ELin

On the GitHub page: the green **Code** button → **Download ZIP**, then unzip
it (on Windows: right-click → *Extract All*) somewhere like your Documents
folder.

## 3. Start it

Open the unzipped folder and:

| | |
| --- | --- |
| **Windows** | double-click **`start-windows.bat`** |
| **macOS / Linux** | run **`./start-macos-linux.sh`** in a Terminal opened in that folder |

On Windows you may get a blue *"Windows protected your PC"* box the first
time, because the file came from the internet: click **More info** → **Run
anyway**.

A black window opens and says what it is doing. **The first run takes a few
minutes** — it downloads the scientific libraries and builds the manual.
When it is ready, your browser opens on the tool. If it does not, go to
<http://127.0.0.1:8000>.

That black window *is* the program: keep it open while you work, and close
it (or press `Ctrl+C`) to stop. Next time, double-click the same file — it
starts in seconds.

## What to look at first

1. **Network** — choose a preset, e.g. *CIGRE MV interconnected*. The
   diagram appears immediately; click any element to see and edit its
   parameters.
2. **Power Flow** — *Run power flow*. The network is coloured by the
   results; hover an element to read its own.
3. **Modal Analysis → Eigenvalue map** — the same network linearised: every
   mode, with frequency and damping.
4. **EMT Simulation** — pick a disturbance and run it. This is the slow one:
   a large network takes a few minutes.

Every figure has **PNG / SVG / CSV** buttons in its corner, and every table a
**CSV** one.

**Keeping your work.** On the Network page, *My networks* saves whatever you
have built or modified under a name, and brings it back next time — even after
closing the browser. Unsaved edits are kept too and offered back on your next
visit. Saved networks live in your browser, so use the same address each time
(`http://127.0.0.1:8000`), and press **Export all (backup)** now and then: that
file can be re-imported here or on another computer.

There is a narrated video tour of the whole interface in
[`Documentation/G2ELin_web_UI_tour.mp4`](Documentation/).

## If it does not start

| What you see | What to do |
| --- | --- |
| `python is not recognised` | Python is not on the PATH: run its installer again, choose *Modify*, tick *"Add python.exe to PATH"* |
| The window flashes and closes | Open the folder, type `cmd` in the address bar, press Enter, then type `start-windows.bat` — the error stays on screen |
| `address already in use` | G2ELin is already running in another window; use that one, or close it first |
| Anything else | Send the last ten lines of the black window — they say what failed |

To remove everything: delete the folder. The only thing it created is
`python/.venv/` inside it.

---

Working from the command line, running the tests, the notebooks, or
deploying it on a server: see
[`python/docs/sphinx/installation.md`](python/docs/sphinx/installation.md),
which is also the **Documentation → Installation and setup** page inside the
tool itself.
