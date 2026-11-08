"""Analyses for "Circuit-aware protection against burst errors in lattice surgery".

Each study module exposes `analyze(data_dir) -> dict`, which recomputes that study's stored analysis from its
per-event failure counts. `lsburst.paper` compares the recomputed values with the numbers printed in the paper.

Modules:
    count_sweep: refined protection-count sweep at d = 9 (Table I, Fig. 2b).
    p1_archive: original cohorts for Fig. 2a (precision_d9), Fig. 5a (precision_d13) and Fig. 5b,c (regime).
    r9_estimated_decoder: estimated-burst and other realistic decoders (break-even burst probability).
    r1_response_models: three stochastic response models at d = 13 (SM S5).
    paper: the paper's reported values and the checks against recomputed results.
"""

from lsburst._io import DATA_DIR, REPO_ROOT

__all__ = ["DATA_DIR", "REPO_ROOT"]
