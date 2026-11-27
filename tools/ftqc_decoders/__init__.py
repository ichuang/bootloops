"""`ftqc_decoders`: Certified zero-float QEC syndrome decoders, effective-distance /
trapping-set verifiers, Stim `.dem` & `sinter` interchange, Quantum Tanner local
single-shot decoding, and audited external decoder backends.
"""

from __future__ import annotations

from .adapters.external_decoders import (
    build_check_matrix_for_external,
    decode_external,
    probe_external_backends,
)
from .adapters.sinter_bridge import BootLoopsSinterDecoder, CompiledBootLoopsDecoder
from .bp_osd import DecodeResult, decode_bp_osd
from .bridges import (
    GadgetDecoderRuleTable,
    decode_open_record_probe,
    estimate_reaction_decoder_latency,
    with_procedural_decoders,
)
from .certify import SCHEMA, build_receipt, evaluate_study, load_problem_from_spec, verify_receipt
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
from .enclosures import (
    enclose_exact_shells,
    enclose_peierls_cluster,
    enclose_stratified_shells,
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
from .tanner_local import (
    coset_reduced_weight,
    decode_tanner_local,
    extract_tanner_neighborhoods,
    simulate_tanner_single_shot_rounds,
)
from .union_find import decode_union_find

VERSION = "0.1.0"

__all__ = [
    "SCHEMA",
    "VERSION",
    "BootLoopsSinterDecoder",
    "CompiledBootLoopsDecoder",
    "DecodeResult",
    "DecodingProblem",
    "GadgetDecoderRuleTable",
    "MWPMResult",
    "RefusalError",
    "SyndromeFormula",
    "bb72_problem",
    "bb144_gross_problem",
    "bit_positions",
    "build_check_matrix_for_external",
    "build_receipt",
    "build_weight0_lrat_proof",
    "certify_operational_radius",
    "coset_reduced_weight",
    "decode_bp_osd",
    "decode_external",
    "decode_mwpm",
    "decode_open_record_probe",
    "decode_relay_bp",
    "decode_sat_mwd",
    "decode_small_set_flip",
    "decode_tanner_local",
    "decode_union_find",
    "decode_with_metachecks",
    "enclose_exact_shells",
    "enclose_peierls_cluster",
    "enclose_stratified_shells",
    "estimate_reaction_decoder_latency",
    "evaluate_study",
    "extract_tanner_neighborhoods",
    "fraction_to_outward_strs",
    "load_problem_from_spec",
    "parse_exact_probability",
    "probe_external_backends",
    "repetition_problem",
    "run_decoder",
    "search_trapping_sets",
    "simulate_tanner_single_shot_rounds",
    "toric_x_problem",
    "verify_failure_witness",
    "verify_lrat_optimality",
    "verify_receipt",
    "with_procedural_decoders",
]
