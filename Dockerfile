FROM python:3.13-slim

# install uv
COPY --from=ghcr.io/astral-sh/uv:0.8.17 /uv /bin/uv
ENV UV_SYSTEM_PYTHON=1

# install git
RUN apt-get update && \
    apt-get install -y --no-install-recommends git && \
    rm -rf /var/lib/apt/lists/*

# copy code
WORKDIR /code
COPY .python-version /code/python-version.py
COPY pyproject.toml /code/pyproject.toml
COPY src /code/src

# install package
RUN uv pip install -e .
