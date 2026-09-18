# G2ELin web app: FastAPI backend + web UI + Sphinx documentation.
#
# Build context is this python/ directory:
#   docker build -t g2elin .
#   docker run -p 10000:10000 g2elin      -> http://localhost:10000
# Render builds it from the repo via ../render.yaml.

# ---- Stage 1: build the Sphinx documentation (served at /manual) ----------
FROM python:3.12-slim AS docs
WORKDIR /src
COPY pyproject.toml README.md ./
COPY src ./src
# api too: the reference pages autodoc g2elin_api, which imports FastAPI.
RUN pip install --no-cache-dir -e ".[docs,api]"
COPY docs ./docs
COPY notebooks ./notebooks
COPY tools ./tools
# build_docs.py runs with -W --keep-going: stale widget outputs saved in one
# notebook raise warnings that make it exit non-zero while still writing the
# full site, so check for the output instead of trusting the exit code.
RUN python tools/build_docs.py || true; test -f docs/sphinx/_build/html/index.html

# ---- Stage 2: runtime ------------------------------------------------------
FROM python:3.12-slim
RUN useradd -m -u 1000 app
USER app
ENV HOME=/home/app \
    PATH=/home/app/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    MPLCONFIGDIR=/tmp/mpl \
    G2ELIN_WARMUP=1 \
    PORT=10000
WORKDIR /home/app/app

COPY --chown=app pyproject.toml README.md ./
COPY --chown=app src ./src
# Editable install on purpose: g2elin_api.main locates web/ and the built docs
# relative to its own source file, which only holds when the package runs
# from this tree rather than from site-packages.
RUN pip install --no-cache-dir --user -e ".[api]"
COPY --chown=app web ./web
COPY --chown=app --from=docs /src/docs/sphinx/_build/html ./docs/sphinx/_build/html

EXPOSE 10000
# Shell form so $PORT (set by the host, e.g. Render) is expanded.
CMD uvicorn g2elin_api.main:app --host 0.0.0.0 --port ${PORT} --timeout-keep-alive 75 --proxy-headers --forwarded-allow-ips "*"
