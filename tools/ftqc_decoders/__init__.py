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
from .mwpm import MWPMResult, decode_mwpm
from .relay_bp import decode_relay_bp

VERSION = "0.1.0"

__all__ = [
    "VERSION",
    "DecodeResult",
    "DecodingProblem",
    "MWPMResult",
    "RefusalError",
    "bb72_problem",
    "bb144_gross_problem",
    "bit_positions",
    "decode_bp_osd",
    "decode_mwpm",
    "decode_relay_bp",
    "fraction_to_outward_strs",
    "parse_exact_probability",
    "repetition_problem",
    "toric_x_problem",
]
