# syntax=docker/dockerfile:1.7
# ============================================================================
# seekr-base: foundation image for SEEKR cloud compute.
#
# Contains: CUDA 12.6 runtime, micromamba, OpenMM 8.4.0 (CUDA build), and the
# seekr package (private repo, cloned at build time via a BuildKit secret).
#
# Build (from the seekrflow repo root):
#   DOCKER_BUILDKIT=1 docker build \
#     --network=hose \
#     -f docker/seekr-base.Dockerfile \
#     --secret id=github_token,src=$HOME/.config/seekr/github_token \
#     -t seekr-base:latest .
#
# The github_token file must contain ONLY a GitHub fine-grained PAT with
# read access to the private seekr repository. The token is mounted only
# during the clone step and never persisted in any image layer.
# ============================================================================
FROM nvidia/cuda:12.6.3-runtime-ubuntu22.04

# --- OS-level dependencies -------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        git \
        bzip2 \
    && rm -rf /var/lib/apt/lists/*

# --- micromamba ------------------------------------------------------------
# Lightweight, fast conda-compatible package manager.
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# --- Python / OpenMM environment ------------------------------------------
# Installed into the base conda prefix so it is on PATH by default.
# cuda-version=12.6 pins OpenMM's CUDA build to match the base image.
RUN micromamba install -y -n base -c conda-forge \
        python=3.11 \
        pip \
        numpy \
        attrs \
        cattrs \
        openmm=8.4.0 \
        cuda-version=12.6 \
        parmed \
        rdkit \
        mdtraj \
        networkx \
    && micromamba clean --all --yes

# --- seekr (private repo) --------------------------------------------------
# Cloned with a BuildKit secret so the token never lands in a layer.
#
# SEEKR_GIT_REF selects the branch/tag to clone (default: main).
# SEEKR_CACHE_BUST forces Docker to re-run this layer even when nothing else
# changed -- pass a changing value (e.g. the latest commit SHA or a timestamp)
# at build time so the clone always picks up the newest seekr code:
#   --build-arg SEEKR_CACHE_BUST="$(date +%s)"
ARG SEEKR_GIT_REF=main
ARG SEEKR_CACHE_BUST=0
RUN --mount=type=secret,id=github_token \
    echo "Cache bust: ${SEEKR_CACHE_BUST}" \
    && GITHUB_TOKEN="$(cat /run/secrets/github_token)" \
    && git clone --depth 1 --branch "${SEEKR_GIT_REF}" \
        "https://x-access-token:${GITHUB_TOKEN}@github.com/lvotapka/seekr.git" \
        /opt/seekr \
    && git -C /opt/seekr rev-parse HEAD \
    && python -m pip install --no-cache-dir /opt/seekr

# --- Sanity: verify OpenMM sees a CUDA platform at build time --------------
# (Build-time check only lists platforms; GPU presence is validated at run.)
RUN python -c "import openmm; print('OpenMM', openmm.version.version); \
print([openmm.Platform.getPlatform(i).getName() \
for i in range(openmm.Platform.getNumPlatforms())])"

CMD ["/bin/bash"]