# potatoes-research

Results, corpus pipeline and implementation code for the thesis: *Peeling Potatoes Exactly and Almost Exactly: Computing the Convex Skull of a Polygon* (Zane Dmitri Pereira, The University of Melbourne, 2026).

This repository is a record of the thesis: the published results, the scripts that produced them, and the implementation code that ran, archived and checked to reproduce them exactly. The code in this repository has been trimmed to only include what was needed to reproduce published results and is licensed under [GPL-3.0-or-later](#licence). To use the implementations in your own work, use the [LGPL-3.0-or-later](https://github.com/AnotherZane/potato-peeling/blob/main/LICENSE) licensed [potato-peeling](https://github.com/AnotherZane/potato-peeling) library instead, which can also be used in proprietary and non-GPL software.

## Contents
| | |
|---|---|
| `thesis_code/` | The C++ implementations and their Python bindings, built as the `potato_peeling` module: Chang & Yap (exact), Aronov et al. (meshed), Cabello et al. (randomised) and refinement |
| `reproduce/` | Scripts to fetch the sources, build the corpus, solve, analyse, and check the archive |
| `reproduce/published/` | The results the thesis reports, unchanged, with `NOTES.md` |
| `db_schema.sql` | Schema of the polygon corpus |
| `Dockerfile` | The build environment, pinned |

## Setup

### With Docker or Podman

```
docker build -t potatoes-research .
docker run --rm -it potatoes-research
```

The image pins the environment the published results were produced and checked in: Ubuntu 26.04, g++ 15.2.0, CMake 4.2.3, CGAL 6.1.1, Boost 1.90, GMP 6.3.0, MPFR 4.2.2 and Python 3.13.14. Run the commands below inside it.

### Without

Install a C++17 compiler, CMake, and the CGAL, Boost, GMP and MPFR development packages, and [uv](https://docs.astral.sh/uv/). `uv sync` then builds the implementations (Python 3.13 or later).

## Building the corpus

```
uv run reproduce/fetch_sources.py
uv run reproduce/build_corpus.py
```

`fetch_sources.py` downloads Natural Earth 1:110m countries (v5.1.1) and ten Google Fonts typefaces into `data/`, checking each file against a pinned SHA-256. If a source has moved, it falls back to the copies archived on Zenodo ([10.5281/zenodo.23938660](https://doi.org/10.5281/zenodo.23938660)). `build_corpus.py` writes the 824 non-convex polygons to `poly_db.sqlite`.

## Reproducing the results

Each `solve_*.py` runs one algorithm over the corpus and each `analyse_*.py` derives one result from stored solves. Every script's docstring says what it produces and how to run it. Each run writes a new timestamped database to `reproduce/results/`, with a `.json` file recording the machine, versions and command line. Analyses read the newest matching result in `reproduce/results/`, falling back to `reproduce/published/`, so they can be run against the published results without solving anything.

Some solves are long. Chang & Yap is limited by memory: by default it skips polygons above n = 70 and gives each 600 s, and the published run solved 545 polygons in about 7.5 hours. Refinement takes about 12 hours over both starting implementations and all four configurations.

## Checking the archive

```
uv run python reproduce/verify_thesis_code.py all --budget 60
```

This re-solves published rows with `thesis_code/` and compares the exact answers, cheapest rows first within the budget (seconds of published solve time per check). Without `--budget` it re-runs everything, which takes about 20 hours. Any mismatch fails.

`reproduce/published/NOTES.md` covers how the archived code and the published results relate: the commits the results record, the timings, and the corpus hash.

Refinement's search uses `cos`, `sin` and `atan2`. Its answers are always exact and certified, but on a different maths library they can differ from the published ones.

## Licence

GPL-3.0-or-later; see `LICENSE`. The implementations use CGAL packages licensed under the GPL.

The corpus sources keep their own licences. Natural Earth is in the public domain. The fonts are under the SIL Open Font License, except Fontdiner Swanky (Apache 2.0). Their licence files are fetched into `data/fonts/licenses/`.

## Acknowledgement

Claude Opus 5.5 (Anthropic) assisted in creating this repository, including cleaning the archived code and comments, verifying correctness against the published results, updating the reproduce scripts, and creating the build environment. Any decisions, reviews and final edits are the author's.
