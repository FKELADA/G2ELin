# G2ELin web app as a Hugging Face Docker Space (also runs anywhere Docker
# does: `docker build -t g2elin . && docker run -p 7860:7860 g2elin`).
#
# Built from the staging directory tools/deploy_hf_space.py assembles:
# pyproject.toml + README.md + src/ + web/ + the pre-built Sphinx HTML.
FROM python:3.12-slim

# Spaces run the container as uid 1000.
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    MPLCONFIGDIR=/tmp/mpl \
    G2ELIN_WARMUP=1
WORKDIR /home/user/app

COPY --chown=user pyproject.toml README.md ./
COPY --chown=user src ./src
# Editable install on purpose: g2elin_api.main locates web/ and the built docs
# relative to its own source file (python/src/g2elin_api/main.py -> python/),
# which only holds when the package runs from this tree, not site-packages.
RUN pip install --no-cache-dir --user -e ".[api]"

COPY --chown=user web ./web
COPY --chown=user docs/sphinx/_build/html ./docs/sphinx/_build/html

EXPOSE 7860
CMD ["uvicorn", "g2elin_api.main:app", "--host", "0.0.0.0", "--port", "7860", "--timeout-keep-alive", "75"]
