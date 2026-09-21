# Getting started (no coding needed)

This page is for trying G2ELin out without knowing Python. There are three
steps, and nothing is installed outside the folder you download.

## 1. Install Python (once)

Download it from [python.org/downloads](https://www.python.org/downloads/)
and run the installer.

**On Windows, tick "Add python.exe to PATH"** on the first screen of the
installer before clicking Install. That single checkbox is what lets
everything else work.

On macOS, run the installer and then also double-click
*"Install Certificates.command"* in the Python folder that Finder opens.

## 2. Download G2ELin

On the repository page on GitHub: the green **Code** button →
**Download ZIP**. Then unzip it — on Windows, right-click the file →
*Extract All* — somewhere like your Documents folder.

(If you use Git, `git clone` the repository instead; it is the same thing.)

## 3. Start it

Open the folder you just unzipped and:

- **Windows** — double-click **`start-windows.bat`**.
  Windows may show a blue *"Windows protected your PC"* box the first time,
  because the file was downloaded: click *More info* → *Run anyway*.
- **macOS / Linux** — open the Terminal in that folder and run
  `./start-macos-linux.sh`.

A black window opens and prints what it is doing. **The first time it takes
a few minutes**: it is downloading the scientific libraries it needs and
building the manual. When it is ready your browser opens on the tool.

If the browser does not open by itself, go to <http://127.0.0.1:8000>.

## While it is running

Keep that black window open — it *is* the program. To stop, close it or
press `Ctrl+C` in it. To use G2ELin again later, double-click the same file:
from then on it starts in seconds.

## Where to go first

The **Home** page introduces the tool and lists the preset networks. A good
first pass:

1. **Network** — pick a preset, for example *CIGRE MV interconnected*. The
   diagram is drawn straight away; click any element to see and edit its
   parameters.
2. **Power Flow** — press *Run power flow*. The network is coloured by the
   results, and hovering an element shows its own.
3. **Modal Analysis → Eigenvalue map** — the same network, linearised: every
   mode with its frequency and damping.
4. **EMT Simulation** — pick a disturbance and press *Run*. This is the slow
   one: a large network takes a few minutes.

**Keeping your work.** On the Network page, *My networks* saves whatever you
have built or modified under a name and brings it back next time, even after
closing the browser; unsaved edits are kept too and offered back on your next
visit. Saved networks live in your browser, so open the tool at the same
address each time (`http://127.0.0.1:8000`), and press **Export all (backup)**
now and then: that file can be re-imported here or on another computer.

There is also a narrated video tour of the whole interface in the
repository's `Documentation/` folder.

## If it does not start

- **"python is not recognised"** — Python is installed but not on the PATH.
  Run its installer again, choose *Modify*, and tick *"Add python.exe to
  PATH"*.
- **The window flashes and disappears** — start it from a terminal instead
  so the message stays on screen: open the folder in Windows Explorer, type
  `cmd` in the address bar, press Enter, then type `start-windows.bat`.
- **"address already in use"** — G2ELin is already running in another
  window; use that one, or close it and start again.
- **Anything else** — send the last ten lines of that black window; they say
  what failed.

{doc}`installation` has the same steps as plain commands, plus the test
suite, the notebooks and the server deployment.
