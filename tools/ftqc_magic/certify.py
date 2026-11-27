"""Deterministic receipt builder & independent replay verifier (`ftqc_magic-receipt-1`)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .bridges import (
    audit_ppr_check_frame,
    to_compose_ancilla_contract,
    to_compose_logical_t_summary,
    to_extractors_magic_record,
    to_uncompute_pricing_kwargs,
)
from .catalog import get_catalog_factory
from .ir import CorrelationMode, MagicRefusalError, MagicSourceContract, PhasePolynomialTarget
from .limits import SignedStabilizerFactory, analytical_lower_bound_n
from .matrix import BinaryFactory
from .pipeline import compose_magic_pipeline
from .qldpc_factory import verify_quadcycle_jump_and_singleshot
from .recycle import compress_active_qubits
from .rings import parse_exact_fraction

SCHEMA = "ftqc_magic-receipt-1"


def canonical_json_bytes(obj: Any) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")


def canonical_json_sha256(obj: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(obj)).hexdigest()


def compute_package_source_sha256() -> str:
    pkg_root = Path(__file__).resolve().parent
    h = hashlib.sha256()
    for p in sorted(pkg_root.rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        rel = p.relative_to(pkg_root).as_posix()
        h.update(rel.encode("utf-8") + b"\0")
        h.update(p.read_bytes() + b"\0")
    return h.hexdigest()


def _export_bundle(contract: MagicSourceContract) -> dict[str, Any]:
    exports: dict[str, Any] = {
        "ftqc_compose_summary": to_compose_logical_t_summary(contract),
        "ftqc_extractors_record": to_extractors_magic_record(contract),
        "ftqc_uncompute_pricing": to_uncompute_pricing_kwargs(contract),
    }
    if contract.correlation_mode == CorrelationMode.INDEPENDENT_SINGLE:
        exports["ftqc_compose_ancilla_contract"] = to_compose_ancilla_contract(contract, atom_tag="t1")
    elif contract.correlation_mode == CorrelationMode.BATCH_CORRELATED:
        exports["ftqc_compose_ancilla_contract_interleaved"] = to_compose_ancilla_contract(
            contract, atom_tag="magic_batch1", interleaving="distinct-instances"
        )
    else:
        exports["ftqc_compose_terminal_only_contract"] = to_compose_ancilla_contract(
            contract, atom_tag="t_cat1", terminal_single_shot=True
        )
    return exports


def build_receipt(
    spec: Mapping[str, Any],
    *,
    crosscheck_external: bool = False,
) -> dict[str, Any]:
    """Build a source-pinned `ftqc_magic-receipt-1` certificate from `spec`."""
    kind = str(spec.get("kind", "pipeline")).strip().lower()
    p_in = parse_exact_fraction(spec.get("p_phys", spec.get("p_input", "1/1000")))

    contract: MagicSourceContract
    details: dict[str, Any] = {}
    optimality_gap: dict[str, Any] | None = None

    if kind == "factory":
        if "catalog" in spec:
            fac = get_catalog_factory(str(spec["catalog"]))
            if not isinstance(fac, BinaryFactory):
                raise ValueError(f"Catalog entry {spec['catalog']!r} is not a BinaryFactory; use kind='signed_stabilizer'")
        else:
            fac = BinaryFactory.from_row_strings(
                output_rows=spec["output_rows"],
                check_rows=spec["check_rows"],
                gamma=spec.get("gamma"),
                name=str(spec.get("name", "custom_binary_factory")),
            )
        rec_info = None
        active_ovr = None
        if bool(spec.get("recycle", False)):
            rsched = compress_active_qubits(
                fac,
                max_bnb_nodes=int(spec.get("max_bnb_nodes", 80_000)),
                certify_lrat=bool(spec.get("certify_lrat", False)) and crosscheck_external,
            )
            active_ovr = rsched.active_qubits
            rec_info = rsched.to_dict()

        contract = fac.to_source_contract(
            p_in,
            mode=str(spec.get("mode", "postselect")),
            active_qubits_override=active_ovr,
            use_sat_if_available=crosscheck_external,
        )
        if spec.get("claimed_distance") is not None and contract.distance < int(spec["claimed_distance"]):
            raise MagicRefusalError(
                "CLAIMED-DISTANCE-REFUTED",
                f"Factory {fac.name!r} certified distance {contract.distance} < claimed {spec['claimed_distance']}",
            )
        lb_dict = analytical_lower_bound_n(contract.target, contract.distance)
        optimality_gap = {
            **lb_dict,
            "actual_n": fac.n,
            "gap_inputs": fac.n - lb_dict["lower_bound_n"],
        }
        details = {
            "factory": contract.metadata,
            "recycling": rec_info,
            "faults_frame_audit": audit_ppr_check_frame(fac) if crosscheck_external else None,
        }

    elif kind == "signed_stabilizer":
        if "catalog" in spec:
            sfac = get_catalog_factory(str(spec["catalog"]))
            if not isinstance(sfac, SignedStabilizerFactory):
                raise ValueError(f"Catalog entry {spec['catalog']!r} is not a SignedStabilizerFactory")
        else:
            sfac = SignedStabilizerFactory(
                n=int(spec["n"]),
                k=int(spec["k"]),
                stabilizers=tuple(spec["stabilizers"]),
                name=str(spec.get("name", "custom_signed_stabilizer")),
            )
        tgt = PhasePolynomialTarget.named(str(spec.get("target", "T")))
        d_claim = int(spec.get("claimed_distance", 2 if sfac.n == 9 else 1))
        contract, col_rep = sfac.certify_collapse(tgt, claimed_distance=d_claim, p_input=p_in)
        lb_dict = analytical_lower_bound_n(tgt, contract.distance)
        optimality_gap = {
            **lb_dict,
            "actual_n": sfac.n,
            "gap_inputs": sfac.n - lb_dict["lower_bound_n"],
        }
        details = {"collapse": col_rep}

    elif kind in ("pipeline", "cultivation", "qldpc"):
        if kind == "cultivation":
            stages = [{"kind": "cultivation", **{k: v for k, v in spec.items() if k != "kind"}}]
        elif kind == "qldpc":
            stages = [{"kind": "qldpc", **{k: v for k, v in spec.items() if k != "kind"}}]
        else:
            stages = list(spec["stages"])
        contract, ledger = compose_magic_pipeline(
            stages, p_phys=p_in, crosscheck_external=crosscheck_external
        )
        details = {"stages": ledger}
        if kind == "qldpc" and bool(spec.get("verify_quadcycle_jump", True)):
            details["quadcycle_jump"] = verify_quadcycle_jump_and_singleshot(
                L=int(spec.get("quadcycle_L", 2)), prior=p_in
            )

    else:
        raise ValueError(f"Unsupported spec kind {kind!r}")

    return {
        "schema": SCHEMA,
        "kind": kind,
        "spec_sha256": canonical_json_sha256(spec),
        "source_sha256": compute_package_source_sha256(),
        "contract": contract.to_dict(),
        "optimality_gap": optimality_gap,
        "exports": _export_bundle(contract),
        "details": details,
    }


def verify_receipt(
    spec: Mapping[str, Any],
    receipt: Mapping[str, Any],
    *,
    strict_source_sha: bool = True,
) -> dict[str, Any]:
    """Independently replay `spec` from scratch and certify exact agreement with `receipt`!"""
    if receipt.get("schema") != SCHEMA:
        raise MagicRefusalError(
            "RECEIPT-SCHEMA-MISMATCH",
            f"Expected schema {SCHEMA!r}, got {receipt.get('schema')!r}",
        )
    expected_spec_sha = canonical_json_sha256(spec)
    if receipt.get("spec_sha256") != expected_spec_sha:
        raise MagicRefusalError(
            "RECEIPT-SPEC-SHA-MISMATCH",
            f"Receipt spec_sha256 {receipt.get('spec_sha256')} != {expected_spec_sha}",
        )
    cur_src_sha = compute_package_source_sha256()
    if strict_source_sha and receipt.get("source_sha256") != cur_src_sha:
        raise MagicRefusalError(
            "RECEIPT-SOURCE-SHA-MISMATCH",
            f"Receipt source_sha256 {receipt.get('source_sha256')} != current package {cur_src_sha}",
        )

    replayed = build_receipt(spec, crosscheck_external=False)
    for field_name in ("kind", "contract", "optimality_gap", "exports"):
        if receipt.get(field_name) != replayed.get(field_name):
            raise MagicRefusalError(
                "RECEIPT-VERIFICATION-MISMATCH",
                f"Replayed field {field_name!r} diverged from supplied receipt!",
                {"field": field_name, "replayed": replayed.get(field_name), "supplied": receipt.get(field_name)},
            )

    return {
        "verified": True,
        "schema": SCHEMA,
        "spec_sha256": expected_spec_sha,
        "source_sha256": cur_src_sha,
        "factory_or_pipeline": replayed["contract"]["name"],
    }
