"""`ftqc_uncompute-receipt-1` certificate builder and standalone zero-float re-verifier.

Enforces BootLoops certification invariants:
  - Strict JSON `parse_float` rejection (`FLOAT-LITERAL-IN-RECEIPT` refusal).
  - Canonical SHA-256 input binding.
  - Independent re-execution of HMR symbolic checks, Fiat-Shamir trajectory bounds,
    disjoint Kraus verification (when requested / small width), spooky pebbling schedules,
    and superposition mask enclosures.

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
from typing import Any

from .hmr import verify_hmr_circuit
from .ir_dag import LOWERABLE_KINDS, ArithDAG
from .ir_kickmix import KickMixCircuit
from .kickmix_sim import simulate_trajectories
from .kraus_oracle import verify_channel_kraus
from .masking import LadderStep, verify_mask_spec, verify_streaming_fibonacci_ladder
from .pebbling import PebbleDAG, PebbleStep, validate_pebbling_schedule
from .suggest import suggest_for_dag, suggest_for_kmx

RECEIPT_SCHEMA = "ftqc_uncompute-receipt-1"


class ReceiptRefusalError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def parse_strict_json(text: str) -> dict[str, Any]:
    """Parse JSON while strictly refusing IEEE-754 floating-point numbers (`N6` gate)."""

    def _reject_float(token: str) -> None:
        raise ReceiptRefusalError(
            "FLOAT-LITERAL-IN-RECEIPT",
            f"Floating-point literal {token!r} forbidden in certified JSON; encode rationals/Decimals as strings",
        )

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            if key in result:
                raise ReceiptRefusalError("DUPLICATE-KEY", key)
            result[key] = value
        return result

    obj = json.loads(text, parse_float=_reject_float, parse_constant=_reject_float,
                     object_pairs_hook=unique_object)
    if not isinstance(obj, dict):
        raise ReceiptRefusalError("INVALID-RECEIPT-ROOT", "Receipt root must be a JSON object")
    return obj


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def build_receipt_for_kmx(
    circuit: KickMixCircuit,
    *,
    source_name: str = "<memory>",
    shots: int = 512,
    run_kraus: bool | None = None,
    salt: str = "ftqc_uncompute-fs-v1",
) -> dict[str, Any]:
    """Certify a `KickMixCircuit` and return an `ftqc_uncompute-receipt-1` dict."""
    hmr_res = verify_hmr_circuit(circuit)
    sim_res = simulate_trajectories(
        circuit,
        shots=shots,
        salt=salt,
        branch_injective=hmr_res.branch_injective,
    )

    do_kraus = (circuit.num_qubits <= 8) if run_kraus is None else bool(run_kraus)
    kraus_dict: dict[str, Any] | None = None
    kraus_ok = True
    if do_kraus:
        # Note: if circuit implements a non-identity map without an external oracle, we only assert
        # trace preservation & HMR symbolic validity unless identity is expected (`sim_res.passed`)!
        if sim_res.passed:
            kr = verify_channel_kraus(circuit)
            kraus_dict = kr.to_dict()
            kraus_ok = kr.exact_match
        else:
            # Circuit transforms basis states (e.g. adder `a, b -> a, a+b`): check trace preservation on |+>
            from .kraus_oracle import CTX, ONE_D, ZERO_D, _vec_norm2, evolve_kraus_ensemble

            dim = 1 << circuit.num_qubits
            inv_n = CTX.divide(ONE_D, CTX.sqrt(Decimal(dim)))
            psi_plus = [(inv_n, ZERO_D)] * dim
            ens = evolve_kraus_ensemble(circuit, psi_plus)
            tr_val = ZERO_D
            for wt, bvec in ens:
                tr_val = CTX.add(tr_val, CTX.multiply(wt, _vec_norm2(bvec)))
            tp_err = CTX.subtract(ONE_D, tr_val).copy_abs()
            kraus_ok = tp_err <= Decimal("1e-24")
            kraus_dict = {
                "mode": "trace_preservation_only",
                "trace_preservation_error": "0" if tp_err <= Decimal("1e-36") else str(tp_err),
                "branch_count": len(ens),
                "num_qubits": circuit.num_qubits,
            }

    mask_reports: list[dict[str, Any]] = []
    masks_ok = True
    for md in circuit.masks:
        m_res = verify_mask_spec(
            md.kind,
            modulus=int(md.params["modulus"]) if "modulus" in md.params else None,
            width=int(md.params["width"]) if "width" in md.params else None,
            runway_bits=int(md.params["runway_bits"]) if "runway_bits" in md.params else None,
        )
        if not m_res.valid:
            masks_ok = False
        mask_reports.append({"register": md.register, **m_res.to_dict()})

    suggestions = [s.to_dict() for s in suggest_for_kmx(circuit)]
    stats = circuit.operation_statistics()

    overall_valid = hmr_res.valid and (sim_res.dirty_resets == 0) and kraus_ok and masks_ok
    verdict = "CERTIFIED" if overall_valid else "REFUSED"

    return {
        "schema": RECEIPT_SCHEMA,
        "target_kind": "kmx_circuit",
        "verification_options": {"shots": shots, "salt": salt, "run_kraus": do_kraus},
        "source_name": source_name,
        "target_sha256": circuit.sha256(),
        "verdict": verdict,
        "summary": {
            "num_qubits": circuit.num_qubits,
            "num_bits": circuit.num_bits,
            "reaction_depth": circuit.reaction_depth(),
            "operation_statistics": stats.to_dict(),
        },
        "hmr_verification": hmr_res.to_dict(),
        "fiat_shamir_trajectories": sim_res.to_dict(),
        "kraus_verification": kraus_dict,
        "mask_verifications": mask_reports,
        "suggestions": suggestions,
    }


def build_receipt_for_spec_dict(
    spec: dict[str, Any],
    *,
    raw_bytes: bytes,
    source_name: str = "<spec.json>",
) -> dict[str, Any]:
    """Certify a JSON specification (`ArithDAG`, `pebbling_game`, or `masked_ladder`)."""
    kind = spec.get("kind") or spec.get("schema")
    target_hash = _sha256_bytes(raw_bytes)

    if kind in ("ftqc_uncompute-dag-1", "arith_dag") or ("registers" in spec and "nodes" in spec):
        dag = ArithDAG.from_dict(spec)
        costs = dag.estimate_costs()
        sugs = [s.to_dict() for s in suggest_for_dag(dag)]
        unlowered = [n.node_id for n in dag.nodes if n.op_kind not in LOWERABLE_KINDS]
        lowered_hmr: dict[str, Any] | None = None
        valid = not unlowered
        if valid:
            kmx = dag.lower_to_kmx()
            h_res = verify_hmr_circuit(kmx)
            lowered_hmr = h_res.to_dict()
            valid = valid and h_res.valid
        for r in dag.registers:
            if r.kind == "mask" or r.mask_kind:
                mk = r.mask_kind or "uniform_mod_n"
                mr = verify_mask_spec(mk, modulus=r.modulus, width=r.width)
                valid = valid and mr.valid
        return {
            "schema": RECEIPT_SCHEMA,
            "target_kind": "arith_dag",
            "source_name": source_name,
            "target_sha256": target_hash,
            "verdict": "CERTIFIED" if valid else "REFUSED",
            **({"unlowered_nodes": unlowered} if unlowered else {}),
            "dag_name": dag.name,
            "estimated_costs": costs,
            "lowered_hmr_verification": lowered_hmr,
            "suggestions": sugs,
        }

    if kind == "pebbling_schedule":
        if "length" in spec:
            dag_p = PebbleDAG.line(
                int(spec["length"]),
                window=int(spec["window"]) if spec.get("window") is not None else None,
                window_ancilla_weight=int(spec.get("window_ancilla_weight", 1)),
            )
        else:
            nodes = tuple(int(v) for v in spec["nodes"])
            preds = {int(k): tuple(int(u) for u in v) for k, v in spec["preds"].items()}
            targets = frozenset(int(t) for t in spec["targets"])
            dag_p = PebbleDAG(nodes=nodes, preds=preds, targets=targets)
        steps = [PebbleStep.from_dict(s) for s in spec["steps"]]
        p_res = validate_pebbling_schedule(
            dag_p,
            steps,
            max_pebbles=int(spec["max_pebbles"]) if spec.get("max_pebbles") is not None else None,
            max_weighted_space=int(spec["max_weighted_space"]) if spec.get("max_weighted_space") is not None else None,
        )
        return {
            "schema": RECEIPT_SCHEMA,
            "target_kind": "pebbling_schedule",
            "source_name": source_name,
            "target_sha256": target_hash,
            "verdict": "CERTIFIED" if p_res.valid else "REFUSED",
            "pebbling_verification": p_res.to_dict(),
        }

    if kind == "masked_ladder":
        steps_l = [LadderStep.from_dict(s) for s in spec["steps"]]
        exp_w = (
            {str(k): int(v) for k, v in spec["expected_control_weights"].items()}
            if spec.get("expected_control_weights") is not None
            else None
        )
        l_res = verify_streaming_fibonacci_ladder(
            int(spec["modulus"]),
            int(spec["eigen_base"]),
            steps_l,
            mask_kind=str(spec.get("mask_kind", "mult_mod_n")),
            expected_control_weights=exp_w,
            runs=int(spec.get("runs", 1)),
            subroutine_kind=spec.get("subroutine_kind"),
            subroutine_steps=int(spec["subroutine_steps"]) if spec.get("subroutine_steps") is not None else None,
        )
        return {
            "schema": RECEIPT_SCHEMA,
            "target_kind": "masked_ladder",
            "source_name": source_name,
            "target_sha256": target_hash,
            "verdict": "CERTIFIED" if l_res.valid else "REFUSED",
            "ladder_verification": l_res.to_dict(),
        }

    raise ReceiptRefusalError("UNKNOWN-SPEC-KIND", f"Cannot certify unrecognized JSON spec kind {kind!r}")


def verify_receipt_file(
    receipt_path: str | Path,
    *,
    target_path: str | Path | None = None,
) -> tuple[bool, list[str]]:
    """Standalone zero-float verifier for `ftqc_uncompute-receipt-1` receipts."""
    rpath = Path(receipt_path)
    raw_text = rpath.read_text(encoding="utf-8")
    receipt = parse_strict_json(raw_text)

    errors: list[str] = []
    if receipt.get("schema") != RECEIPT_SCHEMA:
        errors.append(f"WRONG-SCHEMA: expected {RECEIPT_SCHEMA!r}, got {receipt.get('schema')!r}")
        return False, errors

    verdict = receipt.get("verdict")
    if verdict not in ("CERTIFIED", "REFUSED"):
        errors.append(f"INVALID-VERDICT: {verdict!r}")
        return False, errors

    source = receipt.get("source_name")
    if not isinstance(source, str):
        return False, ["MISSING-SOURCE-NAME"]
    candidates = [Path(target_path)] if target_path is not None else [Path(source), rpath.parent / source]
    target = next((candidate for candidate in candidates if candidate.is_file()), None)
    if target is None:
        return False, ["TARGET-REQUIRED: supply the original input with --target"]
    raw = target.read_bytes()
    if receipt.get("target_kind") == "kmx_circuit":
        circuit = KickMixCircuit.from_kmx(raw.decode("utf-8"))
        if circuit.sha256() != receipt.get("target_sha256"):
            return False, ["TARGET-HASH-MISMATCH"]
        legacy = "verification_options" not in receipt
        options = receipt.get("verification_options")
        if legacy:
            trajectory = receipt.get("fiat_shamir_trajectories")
            if not isinstance(trajectory, dict):
                return False, ["MISSING-TRAJECTORIES"]
            options = {"shots": trajectory.get("shots"), "salt": "ftqc_uncompute-fs-v1",
                       "run_kraus": receipt.get("kraus_verification") is not None}
        if (not isinstance(options, dict) or set(options) != {"shots", "salt", "run_kraus"}
                or type(options["shots"]) is not int or options["shots"] < 1
                or not isinstance(options["salt"], str) or type(options["run_kraus"]) is not bool):
            return False, ["INVALID-VERIFICATION-OPTIONS"]
        rebuilt = build_receipt_for_kmx(circuit, source_name=source, **options)
        if legacy:
            del rebuilt["verification_options"]
    else:
        if _sha256_bytes(raw) != receipt.get("target_sha256"):
            return False, ["TARGET-HASH-MISMATCH"]
        rebuilt = build_receipt_for_spec_dict(parse_strict_json(raw.decode("utf-8")),
                                              raw_bytes=raw, source_name=source)

    # Canonical JSON distinguishes booleans from integers and ignores object key order.
    for key in sorted(set(receipt) | set(rebuilt)):
        if (key not in receipt or key not in rebuilt or
                json.dumps(receipt[key], sort_keys=True) != json.dumps(rebuilt[key], sort_keys=True)):
            errors.append(f"EVIDENCE-MISMATCH: {key}")
    return not errors and verdict == "CERTIFIED", errors
