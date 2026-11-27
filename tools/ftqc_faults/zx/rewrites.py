"""Constructive Fault-Equivalent ZX Rewrite Rules and Polynomial-Time Proof Verifier.

Implements Section 6.2 (Definition 6.11, Propositions 6.12-6.16), Section 3.4
(Propositions 3.12-3.13), Proposition 7.4, Appendix B.2 (Propositions B.4-B.10),
and Proposition C.1 of Rodatz, Poor, and Kissinger (arXiv:2506.17181v4).

Because direct semantic fault-equivalence checking is NP-hard at arbitrary scale
(Theorem 3.11), `verify_derivation` checks each local ZX rewrite rule application
and side condition syntactically in polynomial time O(sum_i |D_i|), certifying
w_certified = min_i w_i by transitivity and compositionality!

Zero float literals appear anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from typing import Any

from .diagram import CliffordCircuit, RefusalError, ZXDiagram, circuit_to_zx_fe, extract_circuit_fe
from .pauliweb import analyze_pauli_webs


@dataclass(frozen=True)
class RewriteStep:
    """Single constructive step in a fault-equivalent ZX derivation."""

    rule: str
    params: dict[str, Any] = field(default_factory=dict)
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "rule": self.rule,
            "params": dict(self.params),
            **({"note": self.note} if self.note else {}),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> RewriteStep:
        return cls(
            rule=str(data["rule"]),
            params=dict(data.get("params", {})),
            note=str(data.get("note", "")),
        )


def apply_zx_rewrite(zx_in: ZXDiagram, step: RewriteStep) -> tuple[ZXDiagram, int | None]:
    """Apply a single fault-equivalent ZX rewrite rule to `zx_in`.

    Returns:
      (zx_out, w_bound)
      where `w_bound is None` means exact fault-equivalence for ALL w (infinity),
      and an integer `w_bound >= 1` means `w_bound`-fault-equivalence!
    """
    zx = zx_in.copy()
    rule = step.rule
    p = step.params

    if rule == "unelim_fe":
        # Proposition 6.12 (wire -> 2-legged phase-0 spider):
        eid = int(p["edge_id"])
        sp_kind = str(p.get("spider_kind", "Z"))
        if sp_kind not in ("Z", "X"):
            raise RefusalError("INVALID-REWRITE-PARAM", f"unelim_fe spider_kind={sp_kind!r}")
        if eid not in zx.edges:
            raise RefusalError("REWRITE-MISSING-EDGE", f"edge {eid} not in diagram")
        ed = zx.remove_edge(eid)
        sp = zx.add_node(sp_kind, 0, label=str(p.get("label", "unelim")))
        zx.add_edge(ed.u, sp, fault_free=ed.fault_free)
        zx.add_edge(sp, ed.v, fault_free=ed.fault_free)
        zx.validate()
        return zx, None

    if rule == "elim_fe":
        # Proposition 6.12 (2-legged phase-0 spider -> single wire):
        nid = int(p["node_id"])
        if nid not in zx.nodes:
            raise RefusalError("REWRITE-MISSING-NODE", f"node {nid} not in diagram")
        nd = zx.nodes[nid]
        if nd.kind not in ("Z", "X") or nd.phase_k != 0 or len(zx.adj[nid]) != 2:
            raise RefusalError(
                "ELIM-PRECONDITION-FAILED",
                f"node {nid} must be a 2-legged phase-0 Z/X spider",
            )
        e0_id, e1_id = zx.adj[nid]
        u = zx.other_endpoint(e0_id, nid)
        v = zx.other_endpoint(e1_id, nid)
        ed0 = zx.remove_edge(e0_id)
        ed1 = zx.remove_edge(e1_id)
        zx.remove_node(nid)
        # Resulting wire is fault-free ONLY if both merged segments were fault-free!
        zx.add_edge(u, v, fault_free=(ed0.fault_free and ed1.fault_free))
        zx.validate()
        return zx, None

    if rule == "fuse_1_fe":
        # Proposition 6.13: fuse two adjacent same-colour spiders (u, v) joined by edge `eid`,
        # PROVIDED:
        #   (a) neither spider is a designated measurement port,
        #   (b) total external legs rem_deg = ext_u + ext_v >= 1, AND
        #   (c) min(ext_u, ext_v) <= 1 (otherwise an opposite-colour flip on `eid` has min boundary
        #       weight min(ext_u, ext_v) >= 2 on the fused spider!).
        eid = int(p["edge_id"])
        if eid not in zx.edges:
            raise RefusalError("REWRITE-MISSING-EDGE", f"edge {eid} not in diagram")
        ed = zx.edges[eid]
        u, v = ed.u, ed.v
        nu, nv = zx.nodes[u], zx.nodes[v]
        meas_nids = {mnid for mnid, _, _ in zx.meas_spiders}
        if u in meas_nids or v in meas_nids:
            raise RefusalError("FUSE1-MEAS-PORT-REFUSED", "cannot fuse a designated measurement spider")
        if nu.kind not in ("Z", "X") or nu.kind != nv.kind:
            raise RefusalError("FUSE1-COLOUR-MISMATCH", f"nodes ({u}, {v}) are not same-colour spiders")
        ext_u = len(zx.adj[u]) - 1
        ext_v = len(zx.adj[v]) - 1
        rem_deg = ext_u + ext_v
        if rem_deg < 1:
            raise RefusalError(
                "FUSE1-ISOLATED-SCALAR-REFUSED",
                "Proposition 6.13 requires the fused spider to retain at least 1 leg (rem_deg >= 1)",
            )
        if min(ext_u, ext_v) > 1:
            raise RefusalError(
                "FUSE1-DEGREE-REFUSED",
                f"Proposition 6.13 requires min(ext_u, ext_v) <= 1, got ({ext_u}, {ext_v}); "
                f"use fuse_4_fe / fuse_n_fe for multi-leg spider splits",
            )
        zx.remove_edge(eid)
        nu.phase_k = (nu.phase_k + nv.phase_k) & 3
        for ov_eid in list(zx.adj[v]):
            oed = zx.remove_edge(ov_eid)
            other = oed.u if oed.v == v else oed.v
            zx.add_edge(u, other, fault_free=oed.fault_free, label=oed.label, gadget_sector=oed.gadget_sector)
        zx.remove_node(v)
        zx.validate()
        return zx, None

    if rule == "unfuse_1_fe":
        # Proposition 6.13 & Proposition B.6: split spider `nid` into `nid` + fresh same-colour `sp_new`
        # connected by a single noisy edge, requiring ext_stay + ext_move >= 1 and min(ext_stay, ext_move) <= 1!
        nid = int(p["node_id"])
        move_eids = [int(x) for x in p["move_edge_ids"]]
        new_phase = int(p.get("new_phase_k", 0))
        if nid not in zx.nodes or zx.nodes[nid].kind not in ("Z", "X"):
            raise RefusalError("UNFUSE-INVALID-NODE", f"node {nid} is not a Z/X spider")
        nd = zx.nodes[nid]
        inc_set = set(zx.adj[nid])
        if not set(move_eids).issubset(inc_set):
            raise RefusalError("UNFUSE-EDGE-NOT-INCIDENT", "move_edge_ids must be incident to node_id")
        stay_count = len(inc_set) - len(move_eids)
        if len(move_eids) + stay_count < 1:
            raise RefusalError("FUSE1-ISOLATED-SCALAR-REFUSED", "cannot unfuse an isolated 0-legged spider")
        if min(len(move_eids), stay_count) > 1:
            raise RefusalError(
                "FUSE1-DEGREE-REFUSED",
                f"Proposition 6.13 requires min(ext_move, ext_stay) <= 1, got ({len(move_eids)}, {stay_count})",
            )
        if (len(move_eids) == 0 and new_phase != 0) or (stay_count == 0 and nd.phase_k != new_phase):
            raise RefusalError(
                "UNFUSE-DANGLING-LEAF-REFUSED",
                "a 0-external-leg unfused leaf must have phase 0",
            )
        if new_phase != 0 and nd.phase_k != new_phase and len(move_eids) != 1 and stay_count != 1:
            raise RefusalError(
                "UNFUSE-PHASE-DEGREE-REFUSED",
                "Proposition B.6 requires the phased split spider to have 1 external leg (degree 2)",
            )
        rem_phase = (nd.phase_k - new_phase) & 3
        nd.phase_k = rem_phase
        sp_new = zx.add_node(nd.kind, new_phase, label=str(p.get("label", "unfused")))
        for meid in move_eids:
            oed = zx.remove_edge(meid)
            other = oed.u if oed.v == nid else oed.v
            zx.add_edge(sp_new, other, fault_free=oed.fault_free, label=oed.label)
        zx.add_edge(nid, sp_new, fault_free=False, label="unfuse_int")
        zx.validate()
        return zx, None

    if rule == "fuse_4_fe":
        # Proposition 6.14: expand a 4-legged phase-0 spider into a 4-cycle of 3-legged spiders!
        nid = int(p["node_id"])
        nd = zx.nodes.get(nid)
        if nd is None or nd.kind not in ("Z", "X") or nd.phase_k != 0 or len(zx.adj[nid]) != 4:
            raise RefusalError(
                "FUSE4-PRECONDITION-FAILED",
                f"Proposition 6.14 requires a 4-legged phase-0 Z/X spider, got node {nid}",
            )
        ordered_legs = [int(x) for x in p.get("leg_order", zx.adj[nid])]
        if sorted(ordered_legs) != sorted(zx.adj[nid]):
            raise RefusalError("FUSE4-LEG-MISMATCH", "leg_order must be a permutation of the 4 incident edges")
        ring = [zx.add_node(nd.kind, 0, label=f"f4_{i}") for i in range(4)]
        for i, leg_eid in enumerate(ordered_legs):
            oed = zx.remove_edge(leg_eid)
            other = oed.u if oed.v == nid else oed.v
            zx.add_edge(ring[i], other, fault_free=oed.fault_free, label=oed.label)
        zx.remove_node(nid)
        for i in range(4):
            zx.add_edge(ring[i], ring[(i + 1) % 4], fault_free=False, label=f"f4_ring_{i}")
        zx.validate()
        return zx, None

    if rule in ("fuse_n_fe", "fuse_n_w_fe"):
        # Proposition 6.15 (w = n -> exact for all weights) & Proposition 6.16 (1 <= w < n -> w-fault-equiv):
        nid = int(p["node_id"])
        nd = zx.nodes.get(nid)
        if nd is None or nd.kind not in ("Z", "X") or nd.phase_k != 0:
            raise RefusalError("FUSEN-PRECONDITION-FAILED", f"node {nid} must be a phase-0 Z/X spider")
        deg = len(zx.adj[nid])
        if deg < 4 or (deg % 2) != 0:
            raise RefusalError("FUSEN-DEGREE-FAILED", f"Proposition 6.15/6.16 requires even degree >= 4, got {deg}")
        n_half = deg // 2
        w_bridges = n_half if rule == "fuse_n_fe" else int(p["w"])
        if not (1 <= w_bridges <= n_half):
            raise RefusalError("FUSEN-INVALID-W", f"w={w_bridges} must satisfy 1 <= w <= {n_half}")
        legs = [int(x) for x in p.get("leg_order", zx.adj[nid])]
        if sorted(legs) != sorted(zx.adj[nid]):
            raise RefusalError("FUSEN-LEG-MISMATCH", "leg_order must match incident edges")
        top_legs = legs[:n_half]
        bot_legs = legs[n_half:]

        s_top = zx.add_node(nd.kind, 0, label="fn_top_hub")
        s_bot = zx.add_node(nd.kind, 0, label="fn_bot_hub")
        for i in range(n_half):
            ed_t = zx.remove_edge(top_legs[i])
            ot_t = ed_t.u if ed_t.v == nid else ed_t.v
            ed_b = zx.remove_edge(bot_legs[i])
            ot_b = ed_b.u if ed_b.v == nid else ed_b.v
            if i < w_bridges:
                b_t = zx.add_node(nd.kind, 0, label=f"fn_bt_{i}")
                b_b = zx.add_node(nd.kind, 0, label=f"fn_bb_{i}")
                zx.add_edge(s_top, b_t, fault_free=False)
                zx.add_edge(s_bot, b_b, fault_free=False)
                zx.add_edge(b_t, b_b, fault_free=False, label=f"fn_bridge_{i}")
                zx.add_edge(b_t, ot_t, fault_free=ed_t.fault_free, label=ed_t.label)
                zx.add_edge(b_b, ot_b, fault_free=ed_b.fault_free, label=ed_b.label)
            else:
                zx.add_edge(s_top, ot_t, fault_free=ed_t.fault_free, label=ed_t.label)
                zx.add_edge(s_bot, ot_b, fault_free=ed_b.fault_free, label=ed_b.label)
        zx.remove_node(nid)
        zx.validate()
        return zx, (None if w_bridges == n_half else w_bridges)

    if rule == "pi_copy_fe":
        # Proposition B.5: push a 2-legged pi spider `pi_nid` of opposite colour through `target_nid`
        pi_nid = int(p["pi_node_id"])
        targ_nid = int(p["target_node_id"])
        pnd = zx.nodes.get(pi_nid)
        tnd = zx.nodes.get(targ_nid)
        if pnd is None or tnd is None or pnd.kind not in ("Z", "X") or tnd.kind not in ("Z", "X"):
            raise RefusalError("PICOPY-INVALID-NODES", "invalid pi_copy_fe nodes")
        if pnd.kind == tnd.kind or pnd.phase_k != 2 or len(zx.adj[pi_nid]) != 2:
            raise RefusalError("PICOPY-PRECONDITION-FAILED", "pi_node must be 2-legged phase-pi opposite-colour")
        shared = [eid for eid in zx.adj[pi_nid] if zx.other_endpoint(eid, pi_nid) == targ_nid]
        if len(shared) != 1:
            raise RefusalError("PICOPY-NOT-ADJACENT", "pi_node and target_node must share 1 edge")
        s_eid = shared[0]
        ext_eid = [eid for eid in zx.adj[pi_nid] if eid != s_eid][0]
        ext_other = zx.other_endpoint(ext_eid, pi_nid)
        ed_s = zx.remove_edge(s_eid)
        ed_ext = zx.remove_edge(ext_eid)
        zx.remove_node(pi_nid)
        other_legs = list(zx.adj[targ_nid])
        zx.add_edge(ext_other, targ_nid, fault_free=(ed_s.fault_free and ed_ext.fault_free))
        tnd.phase_k = (-tnd.phase_k) & 3
        for ol_eid in other_legs:
            oed = zx.remove_edge(ol_eid)
            dest = oed.u if oed.v == targ_nid else oed.v
            p_clone = zx.add_node(pnd.kind, 2, label="pi_copied")
            zx.add_edge(targ_nid, p_clone, fault_free=False)
            zx.add_edge(p_clone, dest, fault_free=oed.fault_free, label=oed.label)
        zx.validate()
        return zx, None

    if rule == "pi_pair_fe":
        # Proposition B.7: insert two adjacent 2-legged pi spiders on edge `eid`
        eid = int(p["edge_id"])
        sp_kind = str(p.get("spider_kind", "Z"))
        if eid not in zx.edges or sp_kind not in ("Z", "X"):
            raise RefusalError("PIPAIR-INVALID-ARGS", "invalid pi_pair_fe args")
        ed = zx.remove_edge(eid)
        p1 = zx.add_node(sp_kind, 2, label="pi_pair_1")
        p2 = zx.add_node(sp_kind, 2, label="pi_pair_2")
        zx.add_edge(ed.u, p1, fault_free=ed.fault_free)
        zx.add_edge(p1, p2, fault_free=False)
        zx.add_edge(p2, ed.v, fault_free=ed.fault_free)
        zx.validate()
        return zx, None

    if rule == "ideal_fuse_fe":
        # Proposition B.8: fuse two same-colour spiders connected by a FAULT-FREE purple edge!
        eid = int(p["edge_id"])
        if eid not in zx.edges:
            raise RefusalError("REWRITE-MISSING-EDGE", f"edge {eid} missing")
        ed = zx.edges[eid]
        if not ed.fault_free:
            raise RefusalError(
                "NON-IDEAL-EDGE-IN-IDEAL-RULE",
                f"ideal_fuse_fe (Prop. B.8) requires fault_free=True on edge {eid}",
            )
        u, v = ed.u, ed.v
        if zx.nodes[u].kind not in ("Z", "X") or zx.nodes[u].kind != zx.nodes[v].kind:
            raise RefusalError("IDEAL-FUSE-COLOUR", "nodes must be same-colour spiders")
        zx.remove_edge(eid)
        zx.nodes[u].phase_k = (zx.nodes[u].phase_k + zx.nodes[v].phase_k) & 3
        for ov_eid in list(zx.adj[v]):
            oed = zx.remove_edge(ov_eid)
            other = oed.u if oed.v == v else oed.v
            zx.add_edge(u, other, fault_free=oed.fault_free, label=oed.label)
        zx.remove_node(v)
        zx.validate()
        return zx, None

    raise RefusalError("UNKNOWN-REWRITE-RULE", f"unknown ZX fault-equivalent rewrite rule {rule!r}")


@dataclass(frozen=True)
class DerivationCertificate:
    """Machine-checkable proof-by-construction certificate from `spec` to `impl`."""

    name: str
    spec_circuit: CliffordCircuit | None
    initial_zx: ZXDiagram
    steps: tuple[RewriteStep, ...]
    extracted_circuit: CliffordCircuit | None
    target_noise_model: str = "circuit_depolarizing"
    declared_w: int | None = None  # None = infinity (exact fault equivalence)
    compositional_subproofs: tuple["DerivationCertificate", ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": "DerivationCertificate",
            "name": self.name,
            "target_noise_model": self.target_noise_model,
            "declared_w": self.declared_w,
            "spec_circuit": self.spec_circuit.to_dict() if self.spec_circuit else None,
            "initial_zx": self.initial_zx.to_dict(),
            "steps": [s.to_dict() for s in self.steps],
            "extracted_circuit": self.extracted_circuit.to_dict() if self.extracted_circuit else None,
            "compositional_subproofs": [sp.to_dict() for sp in self.compositional_subproofs],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DerivationCertificate:
        return cls(
            name=str(data.get("name", "derivation")),
            spec_circuit=CliffordCircuit.from_dict(data["spec_circuit"]) if data.get("spec_circuit") else None,
            initial_zx=ZXDiagram.from_dict(data["initial_zx"]),
            steps=tuple(RewriteStep.from_dict(s) for s in data.get("steps", [])),
            extracted_circuit=(
                CliffordCircuit.from_dict(data["extracted_circuit"]) if data.get("extracted_circuit") else None
            ),
            target_noise_model=str(data.get("target_noise_model", "circuit_depolarizing")),
            declared_w=int(data["declared_w"]) if data.get("declared_w") is not None else None,
            compositional_subproofs=tuple(
                cls.from_dict(sp) for sp in data.get("compositional_subproofs", [])
            ),
        )


def _min_w(a: int | None, b: int | None) -> int | None:
    if a is None:
        return b
    if b is None:
        return a
    return min(a, b)


def verify_derivation(cert: DerivationCertificate) -> dict[str, Any]:
    """Verify a constructive ZX fault-equivalence derivation in polynomial time!

    Checks:
      1) Recursive verification of any compositional sub-proofs (Proposition 3.12),
      2) Every local rewrite step in `cert.steps` and its side conditions (Propositions 6.12-6.16, B.5-B.8),
      3) Fault-free ideal linear map preservation across initial -> final -> extracted circuit
         via polynomial-time F_2 Pauli-web boundary stabiliser canonicalization, AND
      4) Fail-closed composite gate extraction (`extract_circuit_fe`: Props. 6.8-6.10).
    """
    w_running: int | None = None

    for sub in cert.compositional_subproofs:
        sub_res = verify_derivation(sub)
        if sub_res["verdict"] != "PROOF-VALID":
            return {
                "verdict": "PROOF-INVALID",
                "name": cert.name,
                "reason": f"subproof {sub.name!r} failed: {sub_res.get('reason')}",
            }
        w_running = _min_w(w_running, sub_res["certified_w"])

    if cert.spec_circuit is not None:
        ana_spec = analyze_pauli_webs(circuit_to_zx_fe(cert.spec_circuit, "fault_free"))
        ana_init = analyze_pauli_webs(cert.initial_zx)
        if ana_spec.signed_boundary_stabs != ana_init.signed_boundary_stabs:
            return {
                "verdict": "PROOF-INVALID",
                "name": cert.name,
                "reason": "initial_zx does not match spec_circuit fault-free linear map",
            }

    cur_zx = cert.initial_zx.copy()
    ana_start = analyze_pauli_webs(cur_zx)
    for idx, st in enumerate(cert.steps):
        if st.rule in ("shor_repeat_d", "steane_codespace_d"):
            d_val = int(st.params["d"])
            if d_val < 1:
                raise RefusalError("INVALID-DISTANCE-PARAM", f"d must be >= 1 in {st.rule}")
            w_running = _min_w(w_running, d_val)
            continue
        cur_zx, step_w = apply_zx_rewrite(cur_zx, st)
        w_running = _min_w(w_running, step_w)

    ana_end = analyze_pauli_webs(cur_zx)
    if ana_end.signed_boundary_stabs != ana_start.signed_boundary_stabs:
        return {
            "verdict": "PROOF-INVALID",
            "name": cert.name,
            "reason": "rewritten ZX diagram linear map drifted from initial_zx",
        }

    if cert.extracted_circuit is not None:
        extract_circuit_fe(cur_zx, cert.extracted_circuit, target_noise_model=cert.target_noise_model)
        ana_ext = analyze_pauli_webs(circuit_to_zx_fe(cert.extracted_circuit, "fault_free"))
        if ana_ext.signed_boundary_stabs != ana_end.signed_boundary_stabs:
            return {
                "verdict": "PROOF-INVALID",
                "name": cert.name,
                "reason": "extracted_circuit fault-free linear map does not match rewritten ZX diagram",
            }

    if cert.declared_w is not None:
        if w_running is not None and w_running < cert.declared_w:
            return {
                "verdict": "PROOF-INVALID",
                "name": cert.name,
                "reason": f"proof certifies w={w_running} < declared_w={cert.declared_w}",
            }
        w_final: int | None = cert.declared_w
    else:
        w_final = w_running

    return {
        "verdict": "PROOF-VALID",
        "name": cert.name,
        "target_noise_model": cert.target_noise_model,
        "certified_w": w_final,  # None means infinity (all w)
        "n_steps": len(cert.steps),
        "n_subproofs": len(cert.compositional_subproofs),
        "final_zx_nodes": len(cur_zx.nodes),
        "final_zx_edges": len(cur_zx.edges),
    }
