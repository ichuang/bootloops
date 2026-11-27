"""Live bridges from `ftqc_decoders` to sibling `ftqc_gadgets`, `ftqc_audit`,
`ftqc_extractors`, `ftqc_compose`, and `ftqc_uncompute` packages.

Provides:
  1. `GadgetDecoderRuleTable` & `with_procedural_decoders`: plugs `ftqc_decoders`
     directly into `ftqc_gadgets/effects.py` (`Decoder.rule = {'type': 'ftqc_decoders', ...}`),
     eliminating the `2^(n-k)` coset-leader lookup-table wall on `qr:47`, `bch`,
     `hgp`, and `bb72`/`bb144`!
  2. `decode_open_record_probe`: operational (unpostselected) spacetime decoding of
     `ftqc_audit.RecordModel` / `ftqc_extractors` (`gross_circuit`, `reading_audit`,
     HGP) probes with forward replay verification.
  3. `estimate_reaction_decoder_latency`: deterministic classical decoder work &
     backlog latency estimator for `ftqc_uncompute` `HMR` / spooky-pebbling
     `reaction_depth` schedules (cf. B. M. Terhal, Rev. Mod. Phys. 87, 307, 2015, §VII).

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any

from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions
from .effective_distance import run_decoder


def _ensure_tools_on_path() -> None:
    tools_dir = str(Path(__file__).resolve().parents[1])
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)


class GadgetDecoderRuleTable:
    """On-demand cached syndrome decoder table matching `ftqc_gadgets.effects._HierTable`.

    Invoked by `ftqc_gadgets.effects.Analysis` whenever a `Decoder` carries
    `rule = {"type": "ftqc_decoders", "algorithm": "bp_osd", ...}`.
    """

    def __init__(
        self,
        rule: dict[str, Any],
        checks: list[int],
        npos: int,
        per: list[int],
    ) -> None:
        prior = rule.get("prior", "1/1000")
        self.problem = DecodingProblem.from_matrices(
            checks,
            (),
            npos,
            name=str(rule.get("name", "gadget_procedural_decoder")),
            default_prior=prior,
        )
        self.spec: dict[str, Any] = {
            k: v for k, v in rule.items() if k not in ("type", "prior", "name")
        }
        if "algorithm" not in self.spec:
            self.spec["algorithm"] = "bp_osd"
        if self.spec["algorithm"] in ("bp_osd", "bp_lsd") and "osd_order" not in self.spec:
            self.spec["osd_order"] = 2
        self.per = per
        self.cache_effect: dict[int, int | None] = {0: 0}
        self.cache_pattern: dict[int, int | None] = {0: 0}

    def pattern(self, values: int) -> int | None:
        if values in self.cache_pattern:
            return self.cache_pattern[values]
        res = run_decoder(self.problem, values, self.spec)
        pat = res.error_mask if (res.converged and self.problem.syndrome_of(res.error_mask) == values) else None
        if len(self.cache_pattern) < 2_000_000:
            self.cache_pattern[values] = pat
        return pat

    def get(self, values: int) -> int | None:
        if values in self.cache_effect:
            return self.cache_effect[values]
        pat = self.pattern(values)
        if pat is None:
            eff = None
        else:
            eff = 0
            for p_idx in bit_positions(pat):
                eff ^= self.per[p_idx]
        if len(self.cache_effect) < 2_000_000:
            self.cache_effect[values] = eff
        return eff


def with_procedural_decoders(
    gadget: Any,
    *,
    algorithm: str = "bp_osd",
    osd_order: int = 2,
    prior: str = "1/1000",
    extra_options: dict[str, Any] | None = None,
) -> Any:
    """Return a copy of an `ftqc_gadgets.Gadget` with every `Decoder` backed on-demand
    by `ftqc_decoders` (`rule={'type': 'ftqc_decoders', ...}`), bypassing `leader_table`.
    """
    _ensure_tools_on_path()
    from ftqc_gadgets.gadget import Decoder, Gadget

    rule_dict: dict[str, Any] = {
        "type": "ftqc_decoders",
        "algorithm": algorithm,
        "osd_order": osd_order,
        "prior": prior,
        **(extra_options or {}),
    }
    new_decs = [
        Decoder(d.name, d.checks, d.bits, d.kind, d.targets, dict(rule_dict))
        for d in gadget.decoders
    ]
    return Gadget(
        gadget.circuit,
        gadget.inputs,
        gadget.outputs,
        new_decs,
        gadget.body_start,
        f"{gadget.name}+{algorithm}",
        gadget.parts,
    )


def decode_open_record_probe(
    record_model: Any,
    physical_faults: list[tuple[int, str]],
    *,
    decoder_spec: str | dict[str, Any] = "bp_lsd",
    location_prior: str = "1/1000",
) -> dict[str, Any]:
    """Operationally decode a non-zero detector history on an `ftqc_audit.RecordModel`
    (such as an `ftqc_extractors` BB144 `gross_circuit` or `reading_audit` probe) without
    postselection!

    1. Replays `physical_faults` forward through `record_model.replay` to get the true
       detector history `det_bits` and un-decoded observable flip `raw_obs_flip`.
    2. Runs `decoder_spec` on `det_bits` over `DecodingProblem.from_record_model`.
    3. Applies the decoder's predicted observable correction `res.observable_mask`,
       checking whether all logical observables are operationally restored (`corrected_obs_flip == 0`).
    """
    det_bits, raw_obs_flip = record_model.replay(physical_faults)
    prob = DecodingProblem.from_record_model(record_model, location_prior=location_prior)
    res = run_decoder(prob, det_bits, decoder_spec)
    corrected_obs_flip = raw_obs_flip ^ res.observable_mask
    return {
        "circuit_sha256": record_model.circuit.sha256(),
        "decoder": res.decoder,
        "injected_fault_count": len(physical_faults),
        "detector_weight": det_bits.bit_count(),
        "raw_observable_flip_mask": raw_obs_flip,
        "decoder_converged": res.converged,
        "predicted_observable_flip_mask": res.observable_mask,
        "corrected_observable_flip_mask": corrected_obs_flip,
        "logical_success": res.converged and (corrected_obs_flip == 0),
        "decoder_stats": res.stats,
    }


def estimate_reaction_decoder_latency(
    problem: DecodingProblem,
    decoder_spec: str | dict[str, Any],
    reaction_depth: int,
    *,
    syndrome_rounds_per_step: int = 1,
    cycle_ticks: int = 14,
    bp_iter_ticks: int = 1,
    gf2_pivot_ticks: int = 2,
) -> dict[str, Any]:
    """Compute deterministic physical reaction-stall latency for `ftqc_uncompute`
    measurement-based uncomputation (`HMR` / `del_and`, parallel spooky pebbling).

    Measures single-fault worst-case deterministic decoder work (`max_bp_iterations`,
    `max_cluster_rank`, `max_eliminated_rank`) across all `n` single-mechanism syndromes:
      - `BP+LSD` / `Relay-BP` / `UF` incur `O(1)` local cluster rank (`max_cluster_rank <= 2`),
      - Global `BP+OSD` incurs `O(rank(H))` global Gaussian elimination pivots (`Terhal 2015` backlog!).
    """
    if reaction_depth < 0 or syndrome_rounds_per_step < 1 or cycle_ticks < 1:
        raise RefusalError("INVALID-LATENCY-PARAMS", "Invalid reaction latency parameters")

    max_bp_it = 0
    max_elim = 0
    max_clust = 0
    dec_label = str(decoder_spec)

    for j in range(problem.num_mechanisms):
        syn = problem.mech_detectors[j]
        res = run_decoder(problem, syn, decoder_spec)
        dec_label = res.decoder
        st = res.stats
        bp_i = int(st.get("bp_iterations", st.get("total_iterations", 0)))
        el_r = int(st.get("eliminated_rank", 0))
        cl_r = int(st.get("max_cluster_rank", st.get("max_cluster_checks", el_r)))
        if bp_i > max_bp_it:
            max_bp_it = bp_i
        if el_r > max_elim:
            max_elim = el_r
        if cl_r > max_clust:
            max_clust = cl_r

    # Parallel hardware executes independent clusters concurrently (`max_clust` depth vs `max_elim` serial)
    parallel_decoder_ticks = max_bp_it * bp_iter_ticks + max_clust * gf2_pivot_ticks
    serial_decoder_ticks = max_bp_it * bp_iter_ticks + max_elim * gf2_pivot_ticks
    qec_measurement_ticks = syndrome_rounds_per_step * cycle_ticks

    per_step_parallel_ticks = qec_measurement_ticks + parallel_decoder_ticks
    per_step_serial_ticks = qec_measurement_ticks + serial_decoder_ticks

    return {
        "problem_name": problem.name,
        "decoder": dec_label,
        "reaction_depth": reaction_depth,
        "single_shot_work_upper": {
            "max_bp_iterations": max_bp_it,
            "max_parallel_cluster_rank": max_clust,
            "max_global_eliminated_rank": max_elim,
        },
        "ticks_per_reaction_step_parallel": per_step_parallel_ticks,
        "ticks_per_reaction_step_serial": per_step_serial_ticks,
        "total_reaction_stall_ticks_parallel": reaction_depth * per_step_parallel_ticks,
        "total_reaction_stall_ticks_serial": reaction_depth * per_step_serial_ticks,
    }
