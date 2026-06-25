# syntax=docker/dockerfile:1.7
# ============================================================================
# seekr-engines-bd: seekr-base + Browndye2 (built from source).
#
# Use this image for runs that require Brownian dynamics (b-surface / BD
# milestoning). Pure-MD/SEEKR runs can use seekr-base directly and skip the
# extra Browndye build weight.
#
# Build:
#   DOCKER_BUILDKIT=1 docker build \
#     -f docker/seekr-engines-bd.Dockerfile \
#     -t seekr-engines-bd:latest .
# ============================================================================
FROM seekr-base:latest

# Browndye2 build + runtime dependencies (Ubuntu 22.04 set from the
# official instructions: make, gcc, g++, ocaml, libexpat-dev, liblapack-dev,
# apbs). C++17-capable compilers are required.

# TODO: use Opam to install Ocaml instead of this?
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        make \
        cmake \
        gcc \
        g++ \
        ocaml \
        ocaml-findlib \
        libexpat1-dev \
        liblapack-dev \
        apbs \
    && rm -rf /var/lib/apt/lists/*

# Download and compile Browndye2 into /opt/browndye2.
# A single RUN keeps the working directory consistent (each RUN otherwise
# starts a fresh shell back at /, so a bare `cd` would not persist).
WORKDIR /opt/browndye2
RUN curl -Ls https://browndye.ucsd.edu/downloads/browndye2.tar.gz \
        -o /tmp/browndye2.tar.gz \
    && tar xzf /tmp/browndye2.tar.gz -C /opt \
    && cmake -B build \
    && cmake --build build \
    && cmake --install build --prefix /opt/browndye2 \
    && rm -f /tmp/browndye2.tar.gz

# Put Browndye binaries on PATH.
ENV PATH=/opt/browndye2/bin:$PATH

# Sanity check.
RUN which bd_top || echo "WARNING: bd_top not found on PATH"

CMD ["/bin/bash"]