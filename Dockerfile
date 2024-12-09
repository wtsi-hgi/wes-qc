FROM python:3.9.5-slim
WORKDIR /wes-qc
RUN pip install uv
RUN apt-get update && apt-get install -y libpq-dev build-essential wget
COPY pyproject.toml uv.lock /wes-qc
RUN uv sync
COPY . /wes-qc