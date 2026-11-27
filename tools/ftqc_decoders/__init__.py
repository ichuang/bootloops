"""`ftqc_decoders`: Certified zero-float QEC syndrome decoders, effective-distance /
trapping-set verifiers, Stim `.dem` & `sinter` interchange, and audited external
decoder backends (`Tesseract`, `PyMatching`, `ldpc`, `MWPF`, `Chromobius`).
"""

from __future__ import annotations

from .bp_osd import DecodeResult, decode_bp_osd
from .decimal_bounds import RefusalError, fraction_to_outward_strs, parse_exact_probability
from .dem import (
    DecodingProblem,
    bb72_problem,
    bb144_gross_problem,
    bit_positions,
    repetition_problem,
    toric_x_problem,
)
from .effective_distance import (
    certify_operational_radius,
    run_decoder,
    search_trapping_sets,
    verify_failure_witness,
)
from .mwpm import MWPMResult, decode_mwpm
from .relay_bp import decode_relay_bp
from .sat_mwd import (
    SyndromeFormula,
    build_weight0_lrat_proof,
    decode_sat_mwd,
    verify_lrat_optimality,
)
from .ssf_metacheck import decode_small_set_flip, decode_with_metachecks
from .union_find import decode_union_find

VERSION = "0.1.0"

__all__ = [
    "VERSION",
    "DecodeResult",
    "DecodingProblem",
    "MWPMResult",
    "RefusalError",
    "SyndromeFormula",
    "bb72_problem",
    "bb144_gross_problem",
    "bit_positions",
    "build_weight0_lrat_proof",
    "certify_operational_radius",
    "decode_bp_osd",
    "decode_mwpm",
    "decode_relay_bp",
    "decode_sat_mwd",
    "decode_small_set_flip",
    "decode_union_find",
    "decode_with_metachecks",
    "fraction_to_outward_strs",
    "parse_exact_probability",
    "repetition_problem",
    "run_decoder",
    "search_trapping_sets",
    "toric_x_problem",
    "verify_failure_witness",
    "verify_lrat_optimality",
]
