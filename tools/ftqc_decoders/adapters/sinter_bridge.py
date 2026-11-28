"""`sinter.Decoder` / `CompiledDecoder` compatibility wrapper for `ftqc_decoders`.

Allows any deterministic zero-float `ftqc_decoders` algorithm (`bp_osd`, `bp_lsd`,
`relay_bp`, `union_find`, `mwpm`, `ssf`) to plug directly into `sinter.collect(...)`
via `compile_decoder_for_dem(dem)` without temporary disk I/O (Gidney, `sinter` /
`stim` ecosystem standard).

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

import importlib.util
from typing import Any

from ..decimal_bounds import RefusalError
from ..dem import DecodingProblem
from ..effective_distance import run_decoder

if importlib.util.find_spec("sinter") is not None:
    import sinter as _sinter

    _DecoderBase: type = _sinter.Decoder
    _CompiledBase: type = _sinter.CompiledDecoder
else:
    _DecoderBase = object
    _CompiledBase = object


class CompiledBootLoopsDecoder(_CompiledBase):  # type: ignore[misc]
    """In-memory bit-packed shot decoder implementing Sinter's `CompiledDecoder` protocol."""

    def __init__(self, problem: DecodingProblem, decoder_spec: str | dict[str, Any]) -> None:
        self.problem = problem
        self.decoder_spec = decoder_spec

    def decode_shots_bit_packed(self, *, bit_packed_detection_event_data: Any) -> Any:
        import numpy as np

        shots = np.asarray(bit_packed_detection_event_data, dtype=np.uint8)
        n_shots = shots.shape[0]
        n_obs_bytes = (self.problem.num_observables + 7) // 8
        out = np.zeros((n_shots, n_obs_bytes), dtype=np.uint8)

        for row_idx in range(n_shots):
            row_bytes = shots[row_idx].tobytes()
            syn_mask = int.from_bytes(row_bytes, "little") & ((1 << self.problem.num_checks) - 1)
            res = run_decoder(self.problem, syn_mask, self.decoder_spec)
            obs_m = res.observable_mask
            for b_i in range(n_obs_bytes):
                out[row_idx, b_i] = (obs_m >> (8 * b_i)) & 0xFF
        return out


class BootLoopsSinterDecoder(_DecoderBase):  # type: ignore[misc]
    """`sinter.Decoder`-compatible factory compiling a Stim `DetectorErrorModel` into
    `CompiledBootLoopsDecoder`. Note: Sinter serializes custom decoders across worker
    processes when `num_workers > 1`; standard top-level class definition supports pickle!
    """

    def __init__(self, decoder_spec: str | dict[str, Any] = "bp_lsd") -> None:
        self.decoder_spec = decoder_spec

    def compile_decoder_for_dem(self, *, dem: Any) -> CompiledBootLoopsDecoder:
        dem_text = str(dem)
        if not dem_text.strip():
            raise RefusalError("EMPTY-SINTER-DEM", "Received empty DetectorErrorModel in compile_decoder_for_dem")
        problem = DecodingProblem.from_stim_dem_text(dem_text, name="sinter_dem")
        return CompiledBootLoopsDecoder(problem, self.decoder_spec)
