FROM mcr.microsoft.com/devcontainers/python:1-3.12-bookworm

# uv manages the project environment (.venv) from uv.lock
COPY --from=ghcr.io/astral-sh/uv:0.13.7 /uv /uvx /bin/

# Non-Python tools from bioconda, in their own prefix appended to PATH so
# they never shadow the project's Python. Channels use the prefix.dev
# mirrors, as anaconda.org is blocked on the institute network.
COPY --from=mambaorg/micromamba:2.10.0 /bin/micromamba /usr/local/bin/micromamba
RUN micromamba create --yes --prefix /opt/tools --override-channels \
        --channel https://prefix.dev/conda-forge \
        --channel https://prefix.dev/bioconda \
        rtk \
    && micromamba clean --all --yes
ENV PATH="${PATH}:/opt/tools/bin"
