"""Deterministic receipt generation & source-bound verification (`ftqc_decoders-receipt-1`).

Every `ftqc_decoders` study (`decode`, `certify_radius`, `trapping_set`,
`enclose_exact_shells`, `enclose_peierls_cluster`, `enclose_stratified_shells`,
`metacheck_decode`, `reaction_latency`) is serialized into a zero-float, canonical-JSON
SHA-256-bound receipt that can be replayed deterministically from the receipt alone or
cross-verified against an explicit `--target` study file.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any

from .bridges import estimate_reaction_decoder_latency
from .decimal_bounds import RefusalError
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
)
from .enclosures import (
    enclose_exact_shells,
    enclose_peierls_cluster,
    enclose_stratified_shells,
)
from .ssf_metacheck import decode_with_metachecks

SCHEMA = "ftqc_decoders-receipt-1"


def _ensure_tools_on_path() -> None:
    tools_dir = str(Path(__file__).resolve().parents[1])
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)


def canonical_json(obj: Any) -> str:
    """Deterministic NaN-refusing compact JSON serialization."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def _check_no_floats_in_json_tree(obj: Any, path: str = "root") -> None:
    if not isinstance(obj, (str, int, bool, type(None), dict, list)):
        raise RefusalError("FLOAT-INPUT-REFUSED", f"Non-exact numeric type {type(obj).__name__} forbidden at {path}: {obj!r}")
    if isinstance(obj, dict):
        for k, v in obj.items():
            _check_no_floats_in_json_tree(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for idx, v in enumerate(obj):
            _check_no_floats_in_json_tree(v, f"{path}[{idx}]")


def load_problem_from_spec(prob_spec: dict[str, Any], *, base_dir: Path | None = None) -> DecodingProblem:
    """Construct a `DecodingProblem` from a declarative JSON dictionary."""
    _check_no_floats_in_json_tree(prob_spec, "problem")
    kind = str(prob_spec.get("kind", "")).strip()
    prior = prob_spec.get("prior", "1/1000")
    sector = str(prob_spec.get("sector", "X"))

    if kind == "css_code":
        return DecodingProblem.from_css_code(prob_spec["code"], sector=sector, prior=prior)
    if kind == "bb72":
        return bb72_problem(sector=sector, prior=prior)
    if kind == "bb144":
        return bb144_gross_problem(sector=sector, prior=prior)
    if kind == "repetition":
        return repetition_problem(int(prob_spec["distance"]), prior=prior)
    if kind == "toric":
        return toric_x_problem(int(prob_spec["length"]), prior=prior)
    if kind == "stim_dem":
        if "text" in prob_spec:
            dem_txt = str(prob_spec["text"])
        elif "path" in prob_spec:
            p = Path(prob_spec["path"])
            if not p.is_absolute() and base_dir is not None:
                p = base_dir / p
            dem_txt = p.read_text()
        else:
            raise RefusalError("MISSING-DEM-SOURCE", "stim_dem problem requires 'text' or 'path'")
        return DecodingProblem.from_stim_dem_text(
            dem_txt,
            name=str(prob_spec.get("name", "stim_dem")),
            merge_parallel=bool(prob_spec.get("merge_parallel", True)),
        )
    if kind == "matrices":
        return DecodingProblem.from_matrices(
            [int(r) for r in prob_spec["check_rows"]],
            [int(r) for r in prob_spec.get("observable_rows", [])],
            int(prob_spec["num_mechanisms"]),
            name=str(prob_spec.get("name", "matrix_problem")),
            default_prior=prior,
            stabilizer_generators=[int(g) for g in prob_spec.get("stabilizer_generators", [])],
            metacheck_rows=[int(m) for m in prob_spec.get("metacheck_rows", [])],
        )
    if kind == "extractor_gross_probe":
        _ensure_tools_on_path()
        from ftqc_faults.audit.record_model import RecordModel
        from ftqc_extractors.gross_probes import experiment

        rounds = int(prob_spec.get("rounds", 2))
        basis = str(prob_spec.get("basis", "Z"))
        circ = experiment(rounds, basis)
        rm = RecordModel(circ)
        return DecodingProblem.from_record_model(rm, name=f"gross_r{rounds}_{basis}", location_prior=prior)
    if kind == "surgery_cone":
        _ensure_tools_on_path()
        from ftqc_codes import gf2 as fc_gf2
        from ftqc_extractors import surgery_sequence

        p = Path(prob_spec["spec_path"])
        if not p.is_absolute() and base_dir is not None:
            p = base_dir / p
        seq_spec = json.loads(p.read_text())
        op_idx = int(prob_spec.get("operation_index", 0))
        code0 = surgery_sequence.memory(seq_spec["memory"])
        n0, hx0, hz0 = code0
        x_cone, z_cone = list(hx0), list(hz0)
        op = seq_spec["operations"][op_idx]
        cone_w, _imgs = surgery_sequence.append_cone(code0, op, n0, x_cone, z_cone, seq_spec["max_degree"])
        ech_z = fc_gf2.Echelon(z_cone)
        dual_logs = [v for v in fc_gf2.nullspace(x_cone, cone_w) if ech_z.add(v)]
        return DecodingProblem.from_matrices(
            z_cone,
            dual_logs,
            cone_w,
            name=f"surgery_cone_{op_idx}:X",
            default_prior=prior,
            stabilizer_generators=x_cone,
            metacheck_rows=op["metachecks"],
        )

    raise RefusalError("UNKNOWN-PROBLEM-KIND", f"Unsupported problem kind {kind!r}")


def evaluate_study(study: dict[str, Any], *, base_dir: Path | None = None) -> tuple[DecodingProblem, dict[str, Any]]:
    """Execute a declarative decoding study and return `(problem, result_dict)`."""
    _check_no_floats_in_json_tree(study, "study")
    if not isinstance(study, dict) or "study" not in study or "problem" not in study:
        raise RefusalError("INVALID-STUDY-SPEC", "Study JSON must contain 'study' and 'problem' fields")

    problem = load_problem_from_spec(study["problem"], base_dir=base_dir)
    kind = str(study["study"])
    decoder_spec = study.get("decoder", "bp_osd")

    if kind == "decode":
        if "syndrome_checks" in study:
            syn = sum(1 << int(c) for c in study["syndrome_checks"])
        elif "fault_mechanisms" in study:
            e_in = sum(1 << int(m) for m in study["fault_mechanisms"])
            syn = problem.syndrome_of(e_in)
        else:
            syn = int(study.get("syndrome", 0))
        res = run_decoder(problem, syn, decoder_spec)
        out: dict[str, Any] = {
            "status": "proved" if res.converged else "non_converged",
            "decoder": res.decoder,
            "syndrome_checks": bit_positions(syn),
            "converged": res.converged,
            "recovery_mechanisms": bit_positions(res.error_mask),
            "hamming_weight": res.hamming_weight,
            "soft_weight": res.soft_weight,
            "observable_mask": res.observable_mask,
            "stats": res.stats,
        }
        if "fault_mechanisms" in study:
            e_in = sum(1 << int(m) for m in study["fault_mechanisms"])
            obs_res = problem.observable_of(e_in ^ res.error_mask)
            out["logical_success"] = res.converged and (obs_res == 0)
            out["residual_observable_mask"] = obs_res
        return problem, out

    if kind == "certify_radius":
        res_r = certify_operational_radius(
            problem,
            decoder_spec,
            int(study.get("target_radius", 1)),
            budget_cap=int(study.get("budget_cap", 200_000)),
        )
        return problem, res_r

    if kind == "trapping_set":
        res_t = search_trapping_sets(
            problem,
            decoder_spec,
            max_exhaustive_weight=int(study.get("max_exhaustive_weight", 1)),
            max_splits_per_generator=int(study.get("max_splits_per_generator", 32)),
            compare_decoders=study.get("compare_decoders"),
        )
        return problem, res_t

    if kind == "enclose_exact_shells":
        res_e = enclose_exact_shells(
            problem,
            decoder_spec,
            study["physical_p"],
            verified_radius=int(study.get("verified_radius", 0)),
            max_exact_weight=int(study.get("max_exact_weight", 2)),
            budget_cap=int(study.get("budget_cap", 200_000)),
            digits=int(study.get("digits", 30)),
        )
        return problem, res_e

    if kind == "enclose_peierls_cluster":
        res_p = enclose_peierls_cluster(
            problem,
            study["physical_p"],
            verified_component_radius=int(study["verified_component_radius"]),
            cluster_degree_override=study.get("cluster_degree_override"),
            digits=int(study.get("digits", 30)),
        )
        return problem, res_p

    if kind == "enclose_linked_cluster":
        from .linked_cluster import enclose_linked_cluster

        res_l = enclose_linked_cluster(
            problem,
            study.get("base_decoder", "union_find"),
            study["physical_p"],
            detector_range=int(study.get("detector_range", 1)),
            max_cluster_size=int(study.get("max_cluster_size", 2)),
            kappa=int(study.get("observable_mask", 1)),
            budget_cap=int(study.get("budget_cap", 200_000)),
            digits=int(study.get("digits", 30)),
        )
        return problem, res_l

    if kind == "enclose_stratified_shells":
        res_s = enclose_stratified_shells(
            problem,
            decoder_spec,
            study["physical_p"],
            verified_radius=int(study.get("verified_radius", 0)),
            max_shell_weight=int(study.get("max_shell_weight", 4)),
            shots_per_shell=int(study.get("shots_per_shell", 100)),
            failure_delta=study.get("failure_delta", "1/10000"),
            seed=int(study.get("seed", 20261126)),
            digits=int(study.get("digits", 30)),
        )
        return problem, res_s

    if kind == "metacheck_decode":
        e_data = sum(1 << int(m) for m in study.get("data_faults", []))
        e_meas = sum(1 << int(c) for c in study.get("syndrome_flips", []))
        noisy_syn = problem.syndrome_of(e_data) ^ e_meas
        mc_out = decode_with_metachecks(
            problem,
            noisy_syn,
            syndrome_prior=study.get("syndrome_prior", "1/1000"),
            meta_decoder=str(study.get("meta_decoder", "bp_osd")),
            data_decoder=str(study.get("data_decoder", "bp_osd")),
        )
        obs_res = problem.observable_of(e_data ^ mc_out["data_error_mask"])
        mc_out["logical_success"] = bool(mc_out["meta_converged"] and mc_out["data_converged"] and obs_res == 0)
        mc_out["status"] = "proved" if mc_out["logical_success"] else "failed"
        return problem, mc_out

    if kind == "reaction_latency":
        lat_out = estimate_reaction_decoder_latency(
            problem,
            decoder_spec,
            int(study.get("reaction_depth", 10)),
            syndrome_rounds_per_step=int(study.get("syndrome_rounds_per_step", 1)),
            cycle_ticks=int(study.get("cycle_ticks", 14)),
        )
        lat_out["status"] = "proved"
        return problem, lat_out

    raise RefusalError("UNKNOWN-STUDY-KIND", f"Unsupported study kind {kind!r}")


def build_receipt(study: dict[str, Any], *, base_dir: Path | None = None) -> dict[str, Any]:
    """Run `study`, verify zero-float invariants, and emit an `ftqc_decoders-receipt-1` dict."""
    problem, result = evaluate_study(study, base_dir=base_dir)
    receipt = {
        "schema": SCHEMA,
        "study_sha256": sha256_json(study),
        "problem_name": problem.name,
        "problem_sha256": problem.canonical_sha256(),
        "problem_summary": {
            "num_checks": problem.num_checks,
            "num_mechanisms": problem.num_mechanisms,
            "num_observables": problem.num_observables,
            "is_graphlike": problem.is_graphlike(),
            "num_stabilizer_generators": len(problem.stabilizer_generators),
            "num_metachecks": len(problem.metacheck_rows),
        },
        "study_spec": study,
        "result": result,
    }
    _check_no_floats_in_json_tree(receipt, "receipt")
    return receipt


def verify_receipt(
    receipt: dict[str, Any],
    *,
    target_study: dict[str, Any] | None = None,
    base_dir: Path | None = None,
) -> dict[str, Any]:
    """Deterministically regenerate and verify every field of an `ftqc_decoders-receipt-1`."""
    _check_no_floats_in_json_tree(receipt, "receipt")
    if receipt.get("schema") != SCHEMA:
        raise RefusalError("INVALID-RECEIPT-SCHEMA", f"Expected schema {SCHEMA!r}, got {receipt.get('schema')!r}")
    if "study_spec" not in receipt:
        raise RefusalError("MISSING-STUDY-SPEC", "Receipt lacks embedded study_spec")

    spec = target_study if target_study is not None else receipt["study_spec"]
    if sha256_json(spec) != receipt.get("study_sha256"):
        raise RefusalError("STUDY-DIGEST-MISMATCH", "Receipt study_sha256 does not match supplied target study")
    if canonical_json(spec) != canonical_json(receipt["study_spec"]):
        raise RefusalError("STUDY-CONTENT-MISMATCH", "Target study differs from receipt's embedded study_spec")

    regenerated = build_receipt(spec, base_dir=base_dir)
    if canonical_json(regenerated) != canonical_json(receipt):
        raise RefusalError("RECEIPT-REPLAY-MISMATCH", "Recomputed study receipt differs from saved receipt")

    return {
        "verified": True,
        "schema": SCHEMA,
        "study": spec["study"],
        "problem_name": receipt["problem_name"],
        "problem_sha256": receipt["problem_sha256"],
        "status": receipt["result"].get("status", receipt["result"].get("mode", "ok")),
    }
