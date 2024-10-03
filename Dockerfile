FROM python:3.12.7-bookworm

ARG DUCKDB_VERSION="1.1.1"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"
ARG BCFTOOLS_VERSION="1.21"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

RUN pip install --no-cache-dir "duckdb==${DUCKDB_VERSION}"

RUN curl -L -o duckdb_cli-linux-amd64.zip "${DUCKDB_URI}" \
  && unzip duckdb_cli-linux-amd64.zip \
  && cp duckdb /usr/bin/duckdb \
  && rm duckdb_cli-linux-amd64.zip duckdb

RUN curl -L -o bcftools.tar.bz2 "${BCFTOOLS_URI}" \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && make install \
  && cd .. \
  && rm -r "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2
