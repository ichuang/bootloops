"""Addressable high-rate QEC logical-T delivery, occupied-slot swap, and batch correlation checks.

Joins verified `ftqc_magic` single-output `T` receipts (`read_delivery`) with verified
`ftqc_codes` + `ftqc_codes-basis-1` addressable blocks (`LP(6,2)` `c6860dd462fe`),
`ftqc_gadgets.occupied_transfer.occupied_swap`, 4-branch adaptive PBC / Kraus injection,
retry & concurrent wire-tick reconciliation, and side-by-side comparison against the
8-qubit HGP `magic_handoff.consume_t` control.
"""

import copy
from fractions import Fraction
from math import comb

from ftqc_codes import (
    basis as fc_basis,
    certify as fc_certify,
    codes as fc_codes,
    gf2 as fc_gf2,
)
from ftqc_extractors import addressable_surgery
from ftqc_gadgets.occupied_transfer import occupied_swap
from ftqc_gadgets.surgery import Instrument, receipt as surgery_receipt
from ftqc_gadgets.surgery_build import graph as build_gauging_graph
from ftqc_magic.delivery import read_delivery
from .logical_t import instrument as verify_t_kraus_instrument
from .supports import fields


def consume_addressable_lp_t(spec):
    """Verify single-output `T` delivery into a named resource slot of an addressable CSS block,
    followed by 4-branch adaptive `T` consumption onto a surviving named data slot.
    """
    fields(
        spec,
        "kind",
        "code_receipt",
        "basis_receipt",
        "data_slot",
        "resource_slot",
        "factory",
        "factory_receipt",
        "attempts",
        "delivery_request",
        "frame",
        "correction",
        "latency_ticks",
        "transfer_p",
        "online_p",
        "surgery_rounds",
    )
    if spec["kind"] != "addressable-magic-t":
        raise ValueError(f"expected kind 'addressable-magic-t', got {spec['kind']!r}")

    # 1. Verify code + basis receipts and distinct named slots
    code, digests = addressable_surgery.verify_bound_code_and_basis(
        spec["code_receipt"], spec["basis_receipt"]
    )
    basis_rec = spec["basis_receipt"]
    data_slot = spec["data_slot"]
    res_slot = spec["resource_slot"]
    if data_slot == res_slot:
        raise ValueError("data_slot and resource_slot must be distinct logical slots")
    data_idx = fc_basis.slot_index(basis_rec, data_slot)
    res_idx = fc_basis.slot_index(basis_rec, res_slot)

    # 2. Verify ftqc_magic source receipt via `read_delivery` (rejects CCZ-for-T or scalar-only summary!)
    source = read_delivery(spec["factory"], spec["factory_receipt"])
    if source["state"] != "T" or source["groups"] != [[0]] or source["law"] != "exact-accepted-Z":
        raise ValueError("addressable T consumption requires one verified T output with exact-accepted-Z law")
    if source["correlation"] != "INDEPENDENT_SINGLE":
        raise ValueError(
            f"single-slot T handoff requires INDEPENDENT_SINGLE correlation, got {source['correlation']!r}"
        )
    if source["active_qubits"] is None or source["expected_wire_ticks_per_output_upper"] is None:
        raise ValueError("addressable T handoff requires priced physical factory qubits and wire-ticks")

    # 3. Verify frame & correction match (rejects unproved reversed frame!)
    frame = spec["frame"]
    correction = spec["correction"]
    if frame not in ("I", "Z") or correction != frame:
        raise ValueError(f"unproved or mismatched resource frame={frame!r}, correction={correction!r}")

    # 4. Verify occupied-slot Clifford SWAP (`occupied_swap`) preserving untouched data slot & Bell reference!
    checks = addressable_surgery.code_checks_pauli(code)
    basis_pairs = [[s["x"], s["z"]] for s in basis_rec["slots"]]
    swap_proof = occupied_swap({"n": code.n, "checks": checks, "basis": basis_pairs}, res_idx)
    if data_idx not in swap_proof["unchanged_logicals"]:
        raise ValueError("occupied_swap failed to leave data_slot untouched")

    # 5. Validate attempt transcript (including rejected retry attempts!) and arrival tick binding
    latency = spec["latency_ticks"]
    if type(latency) is not int or latency < 0:
        raise ValueError("latency_ticks must be a nonnegative int")

    attempts = spec["attempts"]
    if not isinstance(attempts, list) or not attempts:
        raise ValueError("attempts must be a non-empty list")

    attempt_map = {}
    factory_wire_ticks_realized = 0
    earliest_start = None
    for att in attempts:
        fields(att, "id", "start", "end", "wires", "raw", "accepted")
        aid = att["id"]
        if (
            not isinstance(aid, str)
            or not aid
            or aid in attempt_map
            or type(att["start"]) is not int
            or type(att["end"]) is not int
            or att["start"] < 0
            or att["start"] >= att["end"]
            or len(att["wires"]) != source["active_qubits"]
            or len(set(att["wires"])) != len(att["wires"])
            or min(att["wires"]) < code.n
            or att["raw"] not in att["wires"]
            or type(att["accepted"]) is not bool
        ):
            raise ValueError("invalid factory attempt specification")
        if not att["accepted"] and source["acceptance_lower"] == "1":
            raise ValueError("impossible rejected attempt for deterministic correcting factory")
        attempt_map[aid] = att
        factory_wire_ticks_realized += (att["end"] - att["start"]) * len(att["wires"])
        if earliest_start is None or att["start"] < earliest_start:
            earliest_start = att["start"]

    req = spec["delivery_request"]
    fields(req, "attempt", "slot", "arrival_tick", "accepted")
    if req["slot"] != res_slot:
        raise ValueError("delivery_request.slot does not match resource_slot")
    if req["attempt"] not in attempt_map:
        raise ValueError("delivery_request references unknown attempt")
    chosen_att = attempt_map[req["attempt"]]
    if req["accepted"] is not True or chosen_att["accepted"] is not True:
        raise ValueError("cannot consume a rejected magic-state attempt")
    arrival_tick = req["arrival_tick"]
    if type(arrival_tick) is not int or arrival_tick < chosen_att["end"] + latency:
        raise ValueError("arrival_tick precedes factory completion plus classical latency")

    # 6. Build and verify all 4 adaptive feed-forward branches + surviving data-slot readout on `code`!
    # Step A: joint parity `Z_{data} Z_{resource}` (`p in {0, 1}`).
    # Step B: resource erasure `X_{resource}` (`m in {0, 1}`).
    # Step C: adaptive PBC data readout on `data_slot`:
    #   - when `p == 0`, logical `Z^m` frame -> measure `(-1)^m X_{data}`;
    #   - when `p == 1`, logical `S Z^m` Clifford correction conjugates `X_{data}` to `(-1)^m Y_{data}`!
    rounds = int(spec["surgery_rounds"])
    zz_word = fc_basis.embed_named_pauli(basis_rec, [(data_slot, "Z"), (res_slot, "Z")], sign=0)
    xr_word = fc_basis.embed_named_pauli(basis_rec, [(res_slot, "X")], sign=0)

    def make_graph_inst(t_word):
        px, pz, _ = fc_basis.parse_pauli(t_word, code.n)
        sites = fc_gf2.supp(px | pz)
        edges = [(i, i + 1) for i in range(len(sites) - 1)]
        sp = build_gauging_graph(code.n, checks, t_word, edges, method="cone", rounds=rounds)
        rec = surgery_receipt(sp, order=1, cap=200_000)
        if rec["faults"]["status"] != "proved":
            raise ValueError("adaptive T surgery leg failed distance-2 verification")
        return sp, rec

    zz_spec, zz_rec = make_graph_inst(zz_word)
    xr_spec, xr_rec = make_graph_inst(xr_word)

    branch_verifications = []
    for p_bit in (0, 1):
        for m_bit in (0, 1):
            # Adaptive rotated data-slot measurement after `S^p Z^m` feed-forward!
            read_axis = "Y" if p_bit == 1 else "X"
            data_read_word = fc_basis.embed_named_pauli(basis_rec, [(data_slot, read_axis)], sign=m_bit)
            fc_basis.verify_logical_target(
                code, basis_rec, [(data_slot, read_axis)], data_read_word, sign=m_bit
            )
            dr_spec, dr_rec = make_graph_inst(data_read_word)
            branch_verifications.append(
                {
                    "parity_p": p_bit,
                    "erasure_m": m_bit,
                    "clifford_fixup": ("S*Z" if m_bit else "S") if p_bit else ("Z" if m_bit else "I"),
                    "adaptive_data_measurement": f"{'-' if m_bit else '+'}{read_axis}_{data_slot}",
                    "physical_target": data_read_word,
                    "readout_locations": sum(1 for _ in Instrument(dr_spec).c.locations()),
                    "fault_lower": dr_rec["faults"]["lower"],
                }
            )

    # Also run exact 8th-root-of-unity Kraus completeness & Z-error commutation proof (`logical_t.instrument`)!
    from ftqc_extractors.magic import branches as ideal_t_branches
    kraus_check = verify_t_kraus_instrument(ideal_t_branches(resource_z=0))

    # 7. Metric conversion justification + conditioned & unnormalized failure propagation!
    # Because `source['law'] == 'exact-accepted-Z'`, the accepted single-qubit resource state is an exact
    # convex mixture `(1 - delta) |T><T| + delta Z|T><T|Z` (up to global phase), so `D_tr(rho, |T><T|) == delta`.
    # Under the verified Kraus instrument (`resource_Z_action == 'Z T'`), a `Z` error on the resource commutes
    # through `Z_{data} Z_{resource}` and shifts `m -> m ^ 1`, inducing an exact logical `Z` error on `data_slot`
    # with probability `delta`! For a Pauli channel `(1 - delta) I + delta Z`, half-diamond distance equals `delta`!
    delta_cond = Fraction(source["conditional_trace_distance_upper"])
    acc_lo = Fraction(source["acceptance_lower"])
    p_trans = Fraction(spec["transfer_p"])
    p_online = Fraction(spec["online_p"])
    if not (0 <= p_trans <= 1 and 0 <= p_online <= 1):
        raise ValueError("transfer_p and online_p must lie in [0, 1]")

    frame_ops_cnt = 2 if frame == "Z" else 0
    swap_ops_cnt = len(swap_proof["operations"]) + frame_ops_cnt
    transfer_tail = min(Fraction(1), swap_ops_cnt * p_trans)

    locs_zz = sum(1 for _ in Instrument(zz_spec).c.locations())
    locs_xr = sum(1 for _ in Instrument(xr_spec).c.locations())
    locs_dr = max(b["readout_locations"] for b in branch_verifications)
    d_zz = zz_rec["faults"]["lower"]
    d_xr = xr_rec["faults"]["lower"]
    d_dr = min(b["fault_lower"] for b in branch_verifications)

    online_gate_tail = min(
        Fraction(1),
        comb(locs_zz, d_zz) * (p_online ** d_zz) + comb(locs_xr, d_xr) * (p_online ** d_xr),
    )
    online_read_tail = min(Fraction(1), comb(locs_dr, d_dr) * (p_online ** d_dr))

    cond_gate_failure = min(Fraction(1), delta_cond + transfer_tail + online_gate_tail)
    cond_total_failure = min(Fraction(1), cond_gate_failure + online_read_tail)
    unnorm_gate_failure = acc_lo * cond_gate_failure
    unnorm_total_failure = acc_lo * cond_total_failure

    # 8. Concurrent factory + memory + transfer wire-tick reconciliation!
    swap_ticks = swap_ops_cnt
    consume_start = arrival_tick + swap_ticks
    zz_ticks = len(zz_spec["ticks"])
    xr_ticks = len(xr_spec["ticks"])
    dr_ticks = len(dr_spec["ticks"])
    total_elapsed_ticks = consume_start + zz_ticks + latency + xr_ticks + latency + dr_ticks

    buffer_wire_ticks = arrival_tick - chosen_att["end"]
    # Live LP(6,2) memory (`code.n` data wires) is held alive from `earliest_start` through `total_elapsed_ticks`!
    memory_live_wire_ticks = (total_elapsed_ticks - earliest_start) * code.n
    swap_raw_wire_ticks = swap_ticks  # raw port wire occupied during swap before release
    surgery_ancilla_wire_ticks = (
        zz_ticks * (Instrument(zz_spec).c.n - code.n)
        + xr_ticks * (Instrument(xr_spec).c.n - code.n)
        + dr_ticks * (Instrument(dr_spec).c.n - code.n)
    )
    total_realized_wire_ticks = (
        factory_wire_ticks_realized
        + buffer_wire_ticks
        + memory_live_wire_ticks
        + swap_raw_wire_ticks
        + surgery_ancilla_wire_ticks
    )
    expected_factory_wire_ticks_per_accepted = Fraction(source["expected_wire_ticks_per_output_upper"])

    return {
        "status": "proved",
        "schema": "ftqc_compose-addressable-magic-t-1",
        "code": {"name": code.name, "n": code.n, "k": basis_rec["k"], **digests},
        "slots": {"data_slot": data_slot, "resource_slot": res_slot, "unchanged_slots_during_swap": [basis_rec["slots"][i]["name"] for i in swap_proof["unchanged_logicals"]]},
        "source": source,
        "transfer": {
            "schema": swap_proof["schema"],
            "frame": frame,
            "correction": correction,
            "swap_operations": swap_ops_cnt,
            "arrival_tick": arrival_tick,
            "consume_start_tick": consume_start,
        },
        "instrument": {
            "kraus_verification": kraus_check,
            "branches_verified": len(branch_verifications),
            "branches": branch_verifications,
        },
        "metric_conversion": {
            "source_law": source["law"],
            "justification": (
                "Exact accepted Z-error mixture `(1-delta)|T><T| + delta Z|T><T|Z` maps through "
                "`resource_Z_action == Z*T` to Pauli-Z error rate `delta` on `data_slot`, giving "
                "exact trace-distance and half-diamond equality `1/2 ||E_ideal - T||_diamond == delta`."
            ),
            "fault_free_half_diamond_upper": str(delta_cond),
        },
        "error_propagation": {
            "acceptance_lower": str(acc_lo),
            "expected_attempts_upper": str(Fraction(1, 1) / acc_lo),
            "source_conditional_z_upper": str(delta_cond),
            "transfer_fault_upper": str(transfer_tail),
            "online_gate_fault_upper": str(online_gate_tail),
            "online_readout_fault_upper": str(online_read_tail),
            "conditional_gate_failure_upper": str(cond_gate_failure),
            "conditional_total_failure_upper": str(cond_total_failure),
            "unnormalized_attempt_joint_gate_failure_upper": str(unnorm_gate_failure),
            "unnormalized_attempt_joint_total_failure_upper": str(unnorm_total_failure),
        },
        "spacetime_reconciliation": {
            "attempts_in_transcript": len(attempts),
            "rejected_retries_in_transcript": sum(1 for a in attempts if not a["accepted"]),
            "total_elapsed_ticks": total_elapsed_ticks,
            "factory_realized_wire_ticks": factory_wire_ticks_realized,
            "factory_expected_wire_ticks_per_accepted_output": str(expected_factory_wire_ticks_per_accepted),
            "buffer_wire_ticks": buffer_wire_ticks,
            "memory_live_wire_ticks": memory_live_wire_ticks,
            "swap_raw_wire_ticks": swap_raw_wire_ticks,
            "surgery_ancilla_wire_ticks": surgery_ancilla_wire_ticks,
            "total_realized_wire_ticks": total_realized_wire_ticks,
        },
    }


def reject_unjustified_correlated_batch_as_independent(batch_factory_spec, batch_receipt):
    """Enforce that `BATCH_CORRELATED` or `CATALYST_POISONED` outputs cannot be consumed
    as independent single-slot resources without an explicit joint law (`consume_batch`).
    """
    src = read_delivery(batch_factory_spec, batch_receipt)
    if src["correlation"] != "INDEPENDENT_SINGLE" or len(src["groups"]) != 1:
        raise ValueError(
            f"cannot consume {src['correlation']} multi-output batch (groups={src['groups']}) "
            "as independent single-slot resources; use consume_batch with a joint Z-mask law"
        )
    return src


def compare_lp6_2_vs_hgp8_logical_t_control():
    """Compare addressable `LP(6,2)` (`n=6, k=2`, data+resource in one block) against the
    8-qubit HGP `magic_handoff.consume_t` control (`n=8, k=1` fresh resource block) under
    matched `bk_15_1_3_t` source receipts (`correcting` and `postselect`).
    """
    from ftqc_magic.certify import build_receipt as build_magic_receipt
    from .magic_handoff import consume_t, transfer_for
    from .test_logical_t import example as hgp8_example

    lp6 = fc_codes.load("lp6_2")
    lp6_rec = fc_certify.build(lp6)
    lp6_basis = fc_basis.derive_paired_basis(lp6)

    out = {}
    for mode in ("correcting", "postselect"):
        fac_spec = {"kind": "factory", "catalog": "bk_15_1_3_t", "p_phys": "1/1000", "mode": mode}
        fac_rec = build_magic_receipt(fac_spec)

        # 8-qubit HGP control (`test_logical_t` + `magic_handoff.consume_t`)
        hgp_s = hgp8_example()
        hgp_s.pop("handoff")
        hgp_s.update(
            kind="magic-logical-t",
            factory=fac_spec,
            factory_receipt=fac_rec,
            attempt=f"hgp8-{mode}-1",
            port={"block": "HGP8/resource", "pair": [3, 5], "frame": "I", "correction": "I"},
            encoder_p="1/1000000",
            online_p="1/100000",
        )
        hgp_s["transfer_receipt"] = transfer_for(hgp_s["gate"], hgp_s["port"])
        hgp_res = consume_t(hgp_s)

        # Addressable LP(6,2) (`q0` = data with untouched reference, `q1` = resource)
        # Recycled `bk_15_1_3_t` uses `active_qubits == 5` (`wires = [6,7,8,9,10]`, `15` ticks/attempt = `75` wire-ticks)!
        fac_wires = list(range(6, 6 + int(read_delivery(fac_spec, fac_rec)["active_qubits"])))
        attempts = (
            [{"id": "att-1", "start": 0, "end": 15, "wires": fac_wires, "raw": 6, "accepted": True}]
            if mode == "correcting"
            else [
                {"id": "att-0", "start": 0, "end": 15, "wires": fac_wires, "raw": 6, "accepted": False},
                {"id": "att-1", "start": 15, "end": 30, "wires": fac_wires, "raw": 6, "accepted": True},
            ]
        )
        arr_tick = attempts[-1]["end"] + 2
        lp_spec = {
            "kind": "addressable-magic-t",
            "code_receipt": lp6_rec,
            "basis_receipt": lp6_basis,
            "data_slot": "q0",
            "resource_slot": "q1",
            "factory": fac_spec,
            "factory_receipt": fac_rec,
            "attempts": attempts,
            "delivery_request": {"attempt": "att-1", "slot": "q1", "arrival_tick": arr_tick, "accepted": True},
            "frame": "I",
            "correction": "I",
            "latency_ticks": 2,
            "transfer_p": "1/1000000",
            "online_p": "1/100000",
            "surgery_rounds": 2,
        }
        lp_res = consume_addressable_lp_t(lp_spec)

        out[mode] = {
            "shared_source_half_diamond_upper": lp_res["metric_conversion"]["fault_free_half_diamond_upper"],
            "source_matches_hgp8_control": (
                lp_res["metric_conversion"]["fault_free_half_diamond_upper"]
                == hgp_res["fault_free_logical_T_half_diamond_upper"]
            ),
            "hgp8_control": {
                "resource_block_qubits": hgp_res["encoder"]["memory_qubits"],
                "data_plus_resource_qubits": 2 * hgp_res["encoder"]["memory_qubits"],
                "encoder_ops": hgp_res["encoder"]["operation_count"],
                "conditional_gate_failure_upper": hgp_res["cases"]["gate"]["accepted_failure_upper"],
                "unnormalized_gate_failure_upper": hgp_res["cases"]["gate"]["attempt_joint_accepted_failure_upper"],
            },
            "lp6_2_addressable": {
                "shared_block_qubits": lp_res["code"]["n"],
                "logical_slots": lp_res["code"]["k"],
                "swap_ops": lp_res["transfer"]["swap_operations"],
                "conditional_gate_failure_upper": lp_res["error_propagation"]["conditional_gate_failure_upper"],
                "unnormalized_gate_failure_upper": lp_res["error_propagation"]["unnormalized_attempt_joint_gate_failure_upper"],
                "spacetime": lp_res["spacetime_reconciliation"],
            },
        }
    return out
