@echo off
rem G2ELin - one-click setup and start (Windows).
rem Creates a private Python environment under python\.venv, installs the
rem tool into it, builds the manual once, and opens the web interface.
rem Running it again just starts the interface.
setlocal
title G2ELin
cd /d "%~dp0python"

echo.
echo   G2ELin
echo   ------
echo   The first run installs everything and takes a few minutes.
echo   Later runs start in seconds. Keep this window open while you work.
echo.

rem --- 1. Find Python -------------------------------------------------------
set "PY="
py -3 --version >nul 2>&1 && set "PY=py -3"
if not defined PY (
  python --version >nul 2>&1 && set "PY=python"
)
if not defined PY (
  echo   [X] Python was not found on this computer.
  echo.
  echo       Install Python 3.10 or newer from https://www.python.org/downloads/
  echo       and TICK "Add python.exe to PATH" on the first screen of the installer.
  echo       Then run this file again.
  echo.
  pause
  exit /b 1
)

rem --- 2. A private environment, created once -------------------------------
if not exist ".venv\Scripts\python.exe" (
  echo   [1/4] Creating a private Python environment ...
  %PY% -m venv .venv
  if errorlevel 1 (
    echo   [X] Could not create the environment.
    pause
    exit /b 1
  )
) else (
  echo   [1/4] Python environment already set up.
)
set "VENV=.venv\Scripts\python.exe"

rem --- 3. The tool and its dependencies -------------------------------------
echo   [2/4] Installing G2ELin and its dependencies ...
"%VENV%" -m pip install --upgrade pip --quiet --disable-pip-version-check
"%VENV%" -m pip install -e ".[api,docs]" --quiet --disable-pip-version-check
if errorlevel 1 (
  echo   [X] Installation failed. The lines above say why.
  pause
  exit /b 1
)

rem --- 4. The manual, rebuilt when its sources have changed -----------------
rem Built only when missing before, so updating the tool left the interface
rem serving the old manual with nothing to say so. Compare the newest source
rem against what was built (through a file: quoting python -c inside for /f
rem is a minefield).
set "DOCS=stale"
"%VENV%" -c "import pathlib,sys;b=pathlib.Path('docs/sphinx/_build/html/index.html');src=[p for p in pathlib.Path('docs/sphinx').rglob('*') if p.is_file() and '_build' not in p.parts];print('fresh' if b.exists() and src and b.stat().st_mtime>=max(p.stat().st_mtime for p in src) else 'stale')" > "%TEMP%\g2elin_docs.txt" 2>nul
if exist "%TEMP%\g2elin_docs.txt" (
  set /p DOCS=<"%TEMP%\g2elin_docs.txt"
  del "%TEMP%\g2elin_docs.txt" >nul 2>&1
)
if /i "%DOCS%"=="fresh" (
  echo   [3/4] Documentation already up to date.
) else (
  echo   [3/4] Building the documentation ...
  "%VENV%" tools\build_docs.py >nul 2>&1
  if errorlevel 1 echo   [!] The manual did not build - everything else still works.
)

rem --- 5. Serve, on the first free port, and open a browser ------------------
rem Another copy of G2ELin (or anything else) may already hold 8000.
rem (through a file: quoting a python -c call inside for /f is a minefield)
set "PORT=8000"
"%VENV%" -c "import socket;print(next(p for p in range(8000,8100) if socket.socket().connect_ex(('127.0.0.1',p))!=0))" > "%TEMP%\g2elin_port.txt" 2>nul
if exist "%TEMP%\g2elin_port.txt" (
  set /p PORT=<"%TEMP%\g2elin_port.txt"
  del "%TEMP%\g2elin_port.txt" >nul 2>&1
)

echo   [4/4] Starting the web interface on http://127.0.0.1:%PORT%
echo.
echo   Close this window (or press Ctrl+C) to stop G2ELin.
echo.
echo   G2ELin is free and open access. If it helps your work, you can
echo   support it at https://buymeacoffee.com/FadiKelada
echo.
start "" /b cmd /c "ping -n 7 127.0.0.1 >nul & start "" http://127.0.0.1:%PORT%"
"%VENV%" -m uvicorn g2elin_api.main:app --port %PORT%

echo.
echo   G2ELin has stopped.
pause
