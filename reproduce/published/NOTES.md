# Notes on the published results

These are the results the thesis reports, copied unchanged from the thesis code repository.

## Seeded rows

`analyse_quality_*.sqlite` is unchanged and still contains `seeded` rows (824) from a seeded Cabello variant the thesis does not use; its code is not in this repository. Re-running `analyse_quality.py` produces the same table without them.

## git_commit and git_dirty

`git_commit` refers to the thesis code repository, which is private. This repository was started with fresh history, so those commits do not appear in it.

Every run records `git_dirty: True`. The uncommitted changes were to the `reproduce/`
scripts. The solver code was the committed core. Re-solving with the code in `thesis_code/` reproduces the published answers bit for bit: Aronov and Cabello (baseline) over all 824 polygons, Chang & Yap over its 298 cheapest optima, and refinement (`solve_polish`) over all 824 polygons for `base` and 439-554 for each other config, from both starts. `reproduce/verify_thesis_code.py` reruns this. Every `analyse_*` table regenerates identically from the published solve results.

## Chang & Yap timings

The thesis build of Chang & Yap recorded a trace of every solve and wrote a msgpack dump of its tables, from which `solve_chang_and_yap.py` read `n_xi`, `a_entries` and `m_entries`. The code in `thesis_code/` does neither: the solver returns the three counts directly, and the answers and counts are unchanged. The published `elapsed` times include the tracing and dump overhead, so a fresh run measures lower times than the thesis reports.

## Corpus hash

Each database's `meta.corpus_sha1_16` hashes the corpus file of the thesis code repository. `reproduce/build_corpus.py` produces the same polygons (ids, geometry, metadata) in a file with a different schema that removes redundant fields, so a fresh run records a different hash.
