# Rebuilds the thesis solvers (thesis_code/) in the environment the published results were
# produced and re-verified in:
#
#   Ubuntu 26.04, g++ 15.2.0, CMake 4.2.3, CGAL 6.1.1, Boost 1.90, GMP 6.3.0, MPFR 4.2.2,
#   Python 3.13.14, nanobind 3.1.0, scikit-build-core 1.1.2, uv 0.12.2
#
# System packages come from the Ubuntu archive snapshot of 2026-11-07, so these exact
# versions stay installable after the live archive moves on.
#
#   docker build -t potatoes-research .
#   docker run --rm -it potatoes-research
#   # inside: uv run reproduce/fetch_sources.py && uv run reproduce/build_corpus.py
#   #         uv run python reproduce/verify_thesis_code.py all --budget 60

FROM ubuntu:26.04@sha256:6a06111c3b0a8600feadb270aec771dc239bf407b25625f93f44ab00e11f4d3c

ARG SNAPSHOT=20261107T000000Z
# The snapshot service is HTTPS-only, so CA certificates come from the live archive first.
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && sed -i '/^Signed-By:/a Snapshot: enable' /etc/apt/sources.list.d/ubuntu.sources \
    && apt-get update --snapshot "$SNAPSHOT" \
    && apt-get install -y --no-install-recommends --snapshot "$SNAPSHOT" \
        g++=4:15.2.0-5ubuntu1 \
        cmake=4.2.3-2ubuntu2 \
        ninja-build=1.13.2-1 \
        libcgal-dev=6.1.1-2 \
        libboost-dev=1.90.0.1ubuntu3 \
        libgmp-dev=2:6.3.0+dfsg-5ubuntu2 \
        libmpfr-dev=4.2.2-3 \
        git \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.12.2@sha256:069a51314a7bb6031777a9273205fe1b0b19e914ef418207d1338b268df641dd /uv /uvx /usr/local/bin/

WORKDIR /work
COPY . .
RUN uv python install 3.13.14 && uv sync --frozen --python 3.13.14

CMD ["bash"]
