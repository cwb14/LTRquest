# syntax=docker/dockerfile:1
#
# LTRquest — hermetic image.
#
# The detector normally fetches and builds four helper repos into --tools-dir on
# first use. That is convenient on a workstation and wrong in a container: it
# makes the image's behaviour depend on what GitHub served that morning. Here
# they are cloned at pinned commits during the build, and LTRQUEST_TOOLS_DIR
# points every run at the result, so a container never reaches the network.
#
#   docker build -t ltrquest:1.0.1 .
#   docker run --rm -v "$PWD:/data" -w /data ltrquest:1.0.1 \
#       ltrquest --genome genome.fa --proteins prot.fa --threads 8

# ---------------------------------------------------------------- build stage
FROM mambaorg/micromamba:1.5.10-jammy AS build

USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential git zlib1g-dev ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Pinned so a rebuild six months from now produces the same tools.
ARG KMER2LTR_REF=02fe54a6c657dff9fd34c8877c7d5a562ebabf64
ARG TESORTER2_REF=4997e92f175f819b748866f4e66cc3c9f6116f50
ARG TRFMOD_REF=3e891db310124f7e5f7a630a1c006650be9d1f3a
ARG SDUST_REF=89c42cb41ba598e9cfa07c2ef99ae8c08f769b3e

ENV TOOLS=/opt/ltrquest/tools
RUN mkdir -p "$TOOLS" \
 && git clone https://github.com/cwb14/Kmer2LTR.git "$TOOLS/Kmer2LTR" \
 && git -C "$TOOLS/Kmer2LTR" checkout --quiet "$KMER2LTR_REF" \
 && git clone --branch feat/minimap2 https://github.com/cwb14/TEsorter2.git "$TOOLS/TEsorter2" \
 && git -C "$TOOLS/TEsorter2" checkout --quiet "$TESORTER2_REF" \
 && git clone https://github.com/lh3/TRF-mod "$TOOLS/TRF-mod" \
 && git -C "$TOOLS/TRF-mod" checkout --quiet "$TRFMOD_REF" \
 && make -C "$TOOLS/TRF-mod" -f compile.mak \
 && git clone https://github.com/lh3/sdust.git "$TOOLS/sdust" \
 && git -C "$TOOLS/sdust" checkout --quiet "$SDUST_REF" \
 && make -C "$TOOLS/sdust" \
 && find "$TOOLS" -name .git -type d -prune -exec rm -rf {} + \
 && test -x "$TOOLS/TRF-mod/trf-mod" && test -x "$TOOLS/sdust/sdust"

# --------------------------------------------------------------- runtime stage
FROM mambaorg/micromamba:1.5.10-jammy

LABEL org.opencontainers.image.title="LTRquest" \
      org.opencontainers.image.description="Iterative detection of nested LTR retrotransposons" \
      org.opencontainers.image.source="https://github.com/cwb14/LTRquest" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.version="1.0.1"

USER root
# procps: Nextflow's task tracer needs `ps`. awscli: AWS Batch stages S3 inputs
# with the CLI from *inside* the container (conf/awsbatch.config points at it).
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba install -y -n base -c conda-forge awscli \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH=/opt/conda/bin:$PATH \
    LTRQUEST_TOOLS_DIR=/opt/ltrquest/tools \
    MPLBACKEND=Agg \
    PYTHONNOUSERSITE=1

COPY --from=build /opt/ltrquest/tools /opt/ltrquest/tools

COPY . /opt/ltrquest/src
RUN pip install --no-cache-dir /opt/ltrquest/src \
 && rm -rf /opt/ltrquest/src

# Fail the build, not the first user, if the image is wired up wrong.
RUN ltrquest --help > /dev/null \
 && python -c "import ltrquest, matplotlib, numpy; print('ltrquest', ltrquest.__version__)" \
 && for t in gt ltr_finder mmseqs hmmsearch blastn minimap2 miniprot mafft; do \
        command -v "$t" > /dev/null || { echo "missing: $t" >&2; exit 1; }; \
    done \
 && test -f "$LTRQUEST_TOOLS_DIR/Kmer2LTR/Kmer2LTR.py" \
 && test -f "$LTRQUEST_TOOLS_DIR/Kmer2LTR/flag_fp_families.py" \
 && test -x "$LTRQUEST_TOOLS_DIR/TRF-mod/trf-mod"

USER $MAMBA_USER
WORKDIR /data
ENTRYPOINT []
CMD ["ltrquest", "--help"]
