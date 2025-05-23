FROM python:3.13.2-bookworm AS base

ARG DUCKDB_VERSION="1.2.1"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"
ARG BCFTOOLS_VERSION="1.21"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
ARG HTSLIB_VERSION="1.21"
ARG HTSLIB_URI="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

ENV PIP_ROOT_USER_ACTION=ignore

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gawk \
    zstd \
  && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir "duckdb==${DUCKDB_VERSION}"

RUN curl -L -o duckdb_cli-linux-amd64.zip "${DUCKDB_URI}" \
  && unzip duckdb_cli-linux-amd64.zip -d /usr/local/bin \
  && rm duckdb_cli-linux-amd64.zip

RUN curl -L -o bcftools.tar.bz2 "${BCFTOOLS_URI}" \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && ./configure --prefix=/usr/local \
  && make \
  && make install \
  && cd .. \
  && rm -fr "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2

RUN curl -L -o htslib.tar.bz2 "${HTSLIB_URI}" \
  && tar -jxf htslib.tar.bz2 \
  && cd "htslib-${HTSLIB_VERSION}" \
  && ./configure --prefix=/usr/local \
    --enable-libcurl \
    --enable-gcs \
    --enable-s3 \
    --with-libdeflate \
  && make \
  && make install \
  && cd .. \
  && rm -rf "htslib-${HTSLIB_VERSION}" htslib.tar.bz2

CMD ["bash"]

FROM base AS pipeline

RUN mkdir -p /opt/outlier-exclusion/scripts

COPY scripts/*.py /opt/outlier-exclusion/scripts/
COPY scripts/*.awk /opt/outlier-exclusion/scripts/
