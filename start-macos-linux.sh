#!/usr/bin/env bash
# G2ELin - one-click setup and start (macOS / Linux).
# Creates a private Python environment under python/.venv, installs the tool
# into it, builds the manual once, and opens the web interface. Running it
# again just starts the interface.
set -e
cd "$(dirname "$0")/python"

echo
echo "  G2ELin"
echo "  ------"
echo "  The first run installs everything and takes a few minutes."
echo "  Later runs start in seconds. Keep this window open while you work."
echo

# --- 1. Find Python ---------------------------------------------------------
PY=""
for candidate in python3 python; do
  if command -v "$candidate" >/dev/null 2>&1; then PY="$candidate"; break; fi
done
if [ -z "$PY" ]; then
  echo "  [X] Python was not found on this computer."
  echo "      Install Python 3.10 or newer from https://www.python.org/downloads/"
  echo "      and run this file again."
  exit 1
fi

# --- 2. A private environment, created once ---------------------------------
if [ ! -x ".venv/bin/python" ]; then
  echo "  [1/4] Creating a private Python environment ..."
  "$PY" -m venv .venv
else
  echo "  [1/4] Python environment already set up."
fi
VENV=".venv/bin/python"

# --- 3. The tool and its dependencies ---------------------------------------
echo "  [2/4] Installing G2ELin and its dependencies ..."
"$VENV" -m pip install --upgrade pip --quiet --disable-pip-version-check
"$VENV" -m pip install -e ".[api,docs]" --quiet --disable-pip-version-check

# --- 4. The manual, built once ----------------------------------------------
if [ ! -f "docs/sphinx/_build/html/index.html" ]; then
  echo "  [3/4] Building the documentation (once) ..."
  "$VENV" tools/build_docs.py >/dev/null 2>&1 || echo "  [!] The manual did not build - everything else still works."
else
  echo "  [3/4] Documentation already built."
fi

# --- 5. Serve, on the first free port, and open a browser -------------------
# Another copy of G2ELin (or anything else) may already hold 8000.
PORT=$("$VENV" -c "import socket;print(next(p for p in range(8000,8100) if socket.socket().connect_ex(('127.0.0.1',p))!=0))" 2>/dev/null || echo 8000)

echo "  [4/4] Starting the web interface on http://127.0.0.1:$PORT"
echo
echo "  Press Ctrl+C to stop G2ELin."
echo
( sleep 6
  if command -v open >/dev/null 2>&1; then open "http://127.0.0.1:$PORT"
  elif command -v xdg-open >/dev/null 2>&1; then xdg-open "http://127.0.0.1:$PORT"
  fi ) >/dev/null 2>&1 &

exec "$VENV" -m uvicorn g2elin_api.main:app --port "$PORT"
