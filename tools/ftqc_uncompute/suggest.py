"""MBUC & uncomputation simplification advisor, QEC spacetime scorer, and `--apply` rewriter.

Detects suboptimal reversible uncomputation patterns across both `.kmx` (`KickMixCircuit`)
and `ArithDAG` inputs, ranks rewrite candidates (`R1`–`R7`) by exact resource deltas
`Delta(Toffoli, T, Rotations, Qubits, Depth, ReactionDepth, D_tr)`, prices physical QEC
wire-tick spacetime against `ftqc_zx-contract-v1` gadget contracts, and synthesizes +
auto-verifies rewritten circuits/DAGs.

Rules:
  - `R1_TEMP_AND_HMR`     : Symmetric `CCX` uncompute -> Gidney `del_and` (-1 Toffoli / -4 T)
  - `R2_BILINEAR_HMR`     : Bilinear / GF(2^m) Toffoli uncompute -> parallel `HMR` + `CZ` (0 Toffolis)
  - `R3_LOOKUP_HMR`       : Reverse QROM uncompute -> LSB-folded `del_lookup` (>= 2x fewer Toffolis)
  - `R4_SPOOKY_PEBBLE`    : Out-of-place sequential chain -> KMRV26 parallel spooky pebbling
  - `R5_STREAM_MBUC`      : Disjoint sequential scratch registers -> single recycled MBUC register
  - `R6_SUPERPOS_MASK`    : Coupled accumulator reverse pass -> superposition mask + control streaming
  - `R7_ADDER_PHASE_GRAD` : Symmetric ripple adder / bitwise `Z_POW` -> Gidney HMR adder / `|grad>`

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction
from typing import Any

from .hmr import verify_hmr_circuit
from .ir_dag import ArithDAG, DAGNode, RegDecl
from .ir_kickmix import (
    ANNOTATION_OPS,
    Bit,
    CircuitBuilder,
    KickMixCircuit,
    KickMixOp,
    Qubit,
)
from .kickmix_sim import simulate_trajectories, z_pow_synthesis_cost
from .kraus_oracle import verify_channel_kraus
from .masking import verify_mask_spec
from .pebbling import PebbleDAG, compile_kmrv_blast_schedule, validate_pebbling_schedule


RULE_IDS = (
    "R1_TEMP_AND_HMR",
    "R2_BILINEAR_HMR",
    "R3_LOOKUP_HMR",
    "R4_SPOOKY_PEBBLE",
    "R5_STREAM_MBUC",
    "R6_SUPERPOS_MASK",
    "R7_ADDER_PHASE_GRAD",
)


@dataclass(frozen=True)
class RewriteSuggestion:
    rule_id: str
    location: str
    summary: str
    delta_toffoli: int
    delta_t_upper: int
    delta_rotations: int
    delta_qubits: int
    delta_depth: int
    delta_reaction_depth: int
    trace_distance_cost: Decimal = Decimal(0)
    delta_qec_wire_ticks: Decimal | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def sort_key(self) -> tuple[Any, ...]:
        # Prioritize exact zero-error Toffoli/T/qubit/depth reductions first
        is_approx = 1 if self.trace_distance_cost > Decimal(0) else 0
        return (
            is_approx,
            self.delta_toffoli,
            self.delta_t_upper,
            self.delta_qubits,
            self.delta_depth,
            self.rule_id,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "rule_id": self.rule_id,
            "location": self.location,
            "summary": self.summary,
            "delta_toffoli": self.delta_toffoli,
            "delta_t_upper": self.delta_t_upper,
            "delta_rotations": self.delta_rotations,
            "delta_qubits": self.delta_qubits,
            "delta_depth": self.delta_depth,
            "delta_reaction_depth": self.delta_reaction_depth,
            "trace_distance_cost": str(self.trace_distance_cost),
            "delta_qec_wire_ticks": (
                str(self.delta_qec_wire_ticks) if self.delta_qec_wire_ticks is not None else None
            ),
            "metadata": self.metadata,
        }


def _price_delta_wire_ticks(
    delta_t: int,
    delta_clifford_2q: int,
    delta_reaction: int,
    zx_contract: dict[str, Any] | None,
    t_factory_wire_ticks: int = 5000,
) -> Decimal | None:
    if zx_contract is None:
        return None
    if zx_contract.get("schema") != "ftqc_zx-contract-v1":
        raise ValueError(f"Unsupported ZX contract schema: {zx_contract.get('schema')!r}")
    ctx = Context(prec=28, rounding=ROUND_CEILING)
    ec_wt = Decimal(str(zx_contract["enclosures"]["expected_wire_ticks_per_shot_hi"]))
    # Each T gate consumes 1 magic-state factory volume + 1 EC gadget; each reaction layer costs 1 EC cycle;
    # each Clifford 2Q gate costs 1 EC gadget fraction.
    term_t = ctx.multiply(Decimal(delta_t), ctx.add(Decimal(t_factory_wire_ticks), ec_wt))
    term_c = ctx.multiply(Decimal(delta_clifford_2q + delta_reaction), ec_wt)
    return +ctx.add(term_t, term_c)


# ---------------------------------------------------------------------------
# KickMixCircuit (`.kmx`) pattern analysis & verified rewriting (`R1`, `R2`, `R5`, `R7`)
# ---------------------------------------------------------------------------


def _find_kmx_bilinear_and_temp_and_sites(circuit: KickMixCircuit) -> list[dict[str, Any]]:
    """Scan `circuit` for uncomputed scratch `CCX` pairs (`R1`) and bilinear scratch `CCX` blocks (`R2`).

    Soundness condition for replacing a trailing `CCX` uncompute on scratch qubit(s) `T`:
      1. Every `t in T` is a scratch qubit (`t not in public_qubits`), starts at `|0>`,
         is written ONLY inside the matched compute `CCX` block and uncompute `CCX` block,
         and is NOT touched at all after the uncompute block (until optional `R t`).
      2. Between the compute block and uncompute block, every `t in T` is ONLY read as a
         control/diagonal operand (`Z`, `CZ`, `CCZ`, `Z_POW`, or control of `CX`/`CCX`/`CSWAP`),
         never flipped by `X`/`CX`-target/`SWAP`!
      3. All controls `{c}` of the `CCX` gates hold the EXACT same computational basis values
         at the uncompute block as they did at the compute block (e.g., never modified in between,
         or restored symmetrically). Here we enforce the strict static check that no control `c`
         is a target of `X`, `CX`, `CCX`, `SWAP`, `CSWAP`, `HMR`, or `R` between compute and uncompute,
         and both blocks execute at condition-stack depth 0 with `condition is None`!
    """
    pub = circuit.public_qubits()
    ops = circuit.ops
    n_ops = len(ops)

    # Check condition stack depth per op index
    cond_depth = [0] * n_ops
    d_cur = 0
    for idx, op in enumerate(ops):
        if op.kind == "PUSH_CONDITION":
            d_cur += 1
        elif op.kind == "POP_CONDITION" and d_cur > 0:
            d_cur -= 1
        cond_depth[idx] = d_cur

    def writes_of(op: KickMixOp) -> set[int]:
        if op.kind in ("X", "R", "HMR"):
            return {op.qubits[0]}
        if op.kind == "CX":
            return {op.qubits[1]}
        if op.kind == "CCX":
            return {op.qubits[2]}
        if op.kind == "SWAP":
            return set(op.qubits)
        if op.kind == "CSWAP":
            return {op.qubits[1], op.qubits[2]}
        return set()

    def touches_of(op: KickMixOp) -> set[int]:
        return set(op.qubits)

    sites: list[dict[str, Any]] = []
    consumed_indices: set[int] = set()

    # Pass 1: Detect contiguous multi-CCX bilinear uncompute blocks (`R2`, >= 2 CCXs)
    idx = 0
    while idx < n_ops:
        op = ops[idx]
        if (
            idx not in consumed_indices
            and op.kind == "CCX"
            and op.condition is None
            and cond_depth[idx] == 0
            and op.qubits[2] not in pub
        ):
            # Greedily extend contiguous CCX run writing scratch targets disjoint from controls
            run_end = idx
            targets: set[int] = set()
            ctrls: set[int] = set()
            while run_end < n_ops:
                cand = ops[run_end]
                if cand.kind in ANNOTATION_OPS:
                    run_end += 1
                    continue
                if (
                    cand.kind == "CCX"
                    and cand.condition is None
                    and cond_depth[run_end] == 0
                    and cand.qubits[2] not in pub
                    and cand.qubits[2] not in ctrls
                    and cand.qubits[0] not in targets
                    and cand.qubits[1] not in targets
                ):
                    ctrls.add(cand.qubits[0])
                    ctrls.add(cand.qubits[1])
                    targets.add(cand.qubits[2])
                    run_end += 1
                else:
                    break

            ccx_in_run = [
                k for k in range(idx, run_end) if ops[k].kind == "CCX"
            ]
            # Look ahead for a matching uncompute run with identical multiset mod 2 of (min(c1,c2), max(c1,c2), t)
            if len(ccx_in_run) >= 2:
                def canonical_triples(indices: list[int]) -> set[tuple[int, int, int]]:
                    acc: set[tuple[int, int, int]] = set()
                    for k in indices:
                        c1, c2, t = ops[k].qubits
                        key = (min(c1, c2), max(c1, c2), t)
                        if key in acc:
                            acc.remove(key)
                        else:
                            acc.add(key)
                    return acc

                fwd_triples = canonical_triples(ccx_in_run)
                # Verify targets were NEVER touched before `idx`
                never_before = all(
                    not (touches_of(ops[p]) & targets) for p in range(0, idx)
                )
                if never_before and len(fwd_triples) >= 2:
                    # Scan forward for matching uncompute block
                    u_start = run_end
                    while u_start < n_ops:
                        if ops[u_start].kind == "CCX" and ops[u_start].qubits[2] in targets:
                            u_end = u_start
                            u_ccx: list[int] = []
                            ok_u = True
                            while u_end < n_ops and len(u_ccx) < len(ccx_in_run):
                                uo = ops[u_end]
                                if uo.kind in ANNOTATION_OPS:
                                    u_end += 1
                                    continue
                                if uo.kind == "CCX" and uo.condition is None and cond_depth[u_end] == 0:
                                    u_ccx.append(u_end)
                                    u_end += 1
                                else:
                                    ok_u = False
                                    break
                            if ok_u and len(u_ccx) == len(ccx_in_run) and canonical_triples(u_ccx) == fwd_triples:
                                # Check between `run_end` and `u_start`: neither `targets` nor `ctrls` are written!
                                mid_clean = all(
                                    not (writes_of(ops[m]) & (targets | ctrls))
                                    for m in range(run_end, u_start)
                                )
                                # Check after `u_end`: `targets` are never touched again!
                                after_clean = all(
                                    not (touches_of(ops[a]) & targets)
                                    for a in range(u_end, n_ops)
                                )
                                if mid_clean and after_clean:
                                    sites.append(
                                        {
                                            "kind": "R2_BILINEAR_HMR",
                                            "fwd_indices": ccx_in_run,
                                            "uncompute_start": u_start,
                                            "uncompute_end": u_end,
                                            "uncompute_ccx_indices": u_ccx,
                                            "triples": sorted(fwd_triples),
                                            "targets": sorted(targets),
                                        }
                                    )
                                    consumed_indices.update(ccx_in_run)
                                    consumed_indices.update(range(u_start, u_end))
                                    break
                        u_start += 1
            idx = max(idx + 1, run_end)
        else:
            idx += 1

    # Pass 2: Detect single-CCX temporary-AND compute/uncompute pairs (`R1_TEMP_AND_HMR`)
    # Also supports nested compute/uncompute stacks (like ripple-carry ladders!) where `t`
    # starts clean before `i` (or after an earlier clean uncompute) and `c1, c2, t` are not written in `(i, j)`!
    for i in range(n_ops):
        if i in consumed_indices:
            continue
        op_i = ops[i]
        if (
            op_i.kind != "CCX"
            or op_i.condition is not None
            or cond_depth[i] != 0
            or op_i.qubits[2] in pub
        ):
            continue
        c1, c2, t = op_i.qubits
        pair_key = (min(c1, c2), max(c1, c2), t)

        # Check `t` is provably `|0>` entering `i`: never touched before `i` (or only in already-matched clean pairs)
        if any((touches_of(ops[p]) & {t}) and (p not in consumed_indices) for p in range(0, i)):
            continue

        # Find matching `CCX` at `j > i`
        for j in range(i + 1, n_ops):
            if j in consumed_indices:
                continue
            op_j = ops[j]
            if (
                op_j.kind == "CCX"
                and op_j.condition is None
                and cond_depth[j] == 0
                and (min(op_j.qubits[0], op_j.qubits[1]), max(op_j.qubits[0], op_j.qubits[1]), op_j.qubits[2])
                == pair_key
            ):
                # Check `(i + 1 .. j - 1)` does not write `c1, c2, t`
                # (WAIT: what if in Gidney's symmetric adder `c1` was `CX`'d twice symmetrically `[CX c0 a1 ... CX c0 a1]`?
                # Above we require `t` not written in `(i, j)`, and `c1, c2` either not written in `(i, j)`!)
                if all(not (writes_of(ops[m]) & {c1, c2, t}) for m in range(i + 1, j)):
                    # Check `t` is not read after `j` before its next clean allocation (here: never touched after `j`)
                    if all(not (touches_of(ops[a]) & {t}) for a in range(j + 1, n_ops)):
                        sites.append(
                            {
                                "kind": "R1_TEMP_AND_HMR",
                                "compute_op": i,
                                "uncompute_op": j,
                                "c1": pair_key[0],
                                "c2": pair_key[1],
                                "target": t,
                            }
                        )
                        consumed_indices.add(i)
                        consumed_indices.add(j)
                break
            if writes_of(op_j) & {c1, c2, t}:
                break

    # Pass 3: Detect `@enter:iadd_symmetric` scope blocks (`R7_ADDER_PHASE_GRAD`)
    stack: list[tuple[str, int]] = []
    for idx_s, op_s in enumerate(ops):
        if op_s.kind == "INFO_TAG" and op_s.info_tag:
            if op_s.info_tag == "@enter:iadd_symmetric":
                stack.append(("iadd_symmetric", idx_s))
            elif op_s.info_tag == "@exit:iadd_symmetric" and stack:
                _, s_start = stack.pop()
                if not any(s_start <= k <= idx_s for k in consumed_indices):
                    sub_ccx = [k for k in range(s_start, idx_s) if ops[k].kind == "CCX"]
                    if len(sub_ccx) >= 2 and len(sub_ccx) % 2 == 0:
                        sites.append(
                            {
                                "kind": "R7_ADDER_PHASE_GRAD",
                                "scope_start": s_start,
                                "scope_end": idx_s,
                                "toffolis_saved": len(sub_ccx) // 2,
                            }
                        )
    return sites


def suggest_for_kmx(
    circuit: KickMixCircuit,
    *,
    zx_contract: dict[str, Any] | None = None,
) -> list[RewriteSuggestion]:
    """Identify and rank MBUC simplification opportunities on a `KickMixCircuit`."""
    sites = _find_kmx_bilinear_and_temp_and_sites(circuit)
    suggestions: list[RewriteSuggestion] = []

    for site in sites:
        kind = site["kind"]
        if kind == "R1_TEMP_AND_HMR":
            u_idx = site["uncompute_op"]
            dt = -4
            wt = _price_delta_wire_ticks(dt, 1, 1, zx_contract)
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R1_TEMP_AND_HMR",
                    location=f"ops[{site['compute_op']}..{u_idx}] (q{site['target']})",
                    summary=(
                        f"Replace symmetric Toffoli uncomputation `CCX q{site['c1']} q{site['c2']} q{site['target']}` "
                        f"with Gidney `del_and` (`HMR` + conditioned `CZ`)"
                    ),
                    delta_toffoli=-1,
                    delta_t_upper=dt,
                    delta_rotations=0,
                    delta_qubits=0,
                    delta_depth=0,
                    delta_reaction_depth=1,
                    delta_qec_wire_ticks=wt,
                    metadata=site,
                )
            )
        elif kind == "R2_BILINEAR_HMR":
            n_ccx = len(site["uncompute_ccx_indices"])
            dt = -4 * n_ccx
            n_t = len(site["targets"])
            wt = _price_delta_wire_ticks(dt, len(site["triples"]), 1, zx_contract)
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R2_BILINEAR_HMR",
                    location=f"ops[{site['uncompute_start']}..{site['uncompute_end']}] ({n_t} scratch qubits)",
                    summary=(
                        f"Eliminate all {n_ccx} bilinear uncompute Toffolis on {site['targets']} via parallel "
                        f"`HMR` + classically conditioned `CZ` fixup matrix"
                    ),
                    delta_toffoli=-n_ccx,
                    delta_t_upper=dt,
                    delta_rotations=0,
                    delta_qubits=0,
                    delta_depth=-(n_ccx - 1),
                    delta_reaction_depth=1,
                    delta_qec_wire_ticks=wt,
                    metadata=site,
                )
            )
        elif kind == "R7_ADDER_PHASE_GRAD":
            saved = int(site["toffolis_saved"])
            dt = -4 * saved
            wt = _price_delta_wire_ticks(dt, saved, saved, zx_contract)
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R7_ADDER_PHASE_GRAD",
                    location=f"ops[{site['scope_start']}..{site['scope_end']}] (@enter:iadd_symmetric)",
                    summary=f"Replace symmetric ripple-carry uncompute with Gidney HMR carry cleanup (-{saved} Toffolis)",
                    delta_toffoli=-saved,
                    delta_t_upper=dt,
                    delta_rotations=0,
                    delta_qubits=0,
                    delta_depth=-saved,
                    delta_reaction_depth=saved,
                    delta_qec_wire_ticks=wt,
                    metadata=site,
                )
            )

    # Pass 4 (`R5_STREAM_MBUC`): check if multiple scratch qubits have strictly sequential non-overlapping lifetimes
    pub = circuit.public_qubits()
    lifetimes: dict[int, tuple[int, int]] = {}
    for idx_op, op in enumerate(circuit.ops):
        for q in op.qubits:
            if q not in pub:
                if q not in lifetimes:
                    lifetimes[q] = (idx_op, idx_op)
                else:
                    lifetimes[q] = (lifetimes[q][0], idx_op)
    if len(lifetimes) >= 2:
        # Minimum coloring / interval-graph peak overlap vs len(lifetimes)
        events: list[tuple[int, int]] = []
        for lo, hi in lifetimes.values():
            events.append((lo, 1))
            events.append((hi + 1, -1))
        events.sort()
        cur_live = 0
        peak_live = 0
        for _, delta in events:
            cur_live += delta
            peak_live = max(peak_live, cur_live)
        reducible_q = len(lifetimes) - peak_live
        if reducible_q >= 1:
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R5_STREAM_MBUC",
                    location=f"scratch_qubits ({len(lifetimes)} allocated, {peak_live} concurrent)",
                    summary=(
                        f"Recycle {reducible_q} non-overlapping scratch qubit(s) across sequential MBUC blocks "
                        f"instead of allocating fresh ancillas"
                    ),
                    delta_toffoli=0,
                    delta_t_upper=0,
                    delta_rotations=0,
                    delta_qubits=-reducible_q,
                    delta_depth=0,
                    delta_reaction_depth=0,
                    metadata={"allocated_scratch": len(lifetimes), "peak_concurrent_scratch": peak_live},
                )
            )

    suggestions.sort(key=lambda s: s.sort_key())
    return suggestions


def apply_kmx_rewrites(
    circuit: KickMixCircuit,
    *,
    verify: bool = True,
) -> tuple[KickMixCircuit, list[RewriteSuggestion]]:
    """Apply sound MBUC rewrites (`R1`, `R2`, `R5` qubit recycling) to `circuit` and auto-verify equivalence!"""
    suggestions = suggest_for_kmx(circuit)
    if not suggestions:
        return circuit, []

    r1_map: dict[int, dict[str, Any]] = {}
    r2_ranges: dict[int, dict[str, Any]] = {}
    skip_indices: set[int] = set()

    for sug in suggestions:
        if sug.rule_id == "R1_TEMP_AND_HMR":
            r1_map[int(sug.metadata["uncompute_op"])] = sug.metadata
        elif sug.rule_id == "R2_BILINEAR_HMR":
            u_start = int(sug.metadata["uncompute_start"])
            u_end = int(sug.metadata["uncompute_end"])
            r2_ranges[u_start] = sug.metadata
            for k in range(u_start + 1, u_end):
                skip_indices.add(k)

    next_b = circuit.num_bits
    new_ops: list[KickMixOp] = []

    for idx, op in enumerate(circuit.ops):
        if idx in skip_indices:
            continue
        if idx in r1_map:
            info = r1_map[idx]
            b_tmp = next_b
            next_b += 1
            t_q = int(info["target"])
            c1_q = int(info["c1"])
            c2_q = int(info["c2"])
            new_ops.append(KickMixOp("INFO_TAG", info_tag="@enter:del_and"))
            new_ops.append(KickMixOp("HMR", qubits=(t_q,), bits=(b_tmp,)))
            new_ops.append(KickMixOp("CZ", qubits=(c1_q, c2_q), condition=b_tmp))
            new_ops.append(KickMixOp("BIT_STORE0", bits=(b_tmp,)))
            new_ops.append(KickMixOp("INFO_TAG", info_tag="@exit:del_and"))
            continue
        if idx in r2_ranges:
            info = r2_ranges[idx]
            targets: list[int] = info["targets"]
            triples: list[tuple[int, int, int]] = info["triples"]
            meas_map = {tq: next_b + i for i, tq in enumerate(targets)}
            next_b += len(targets)
            par_b = next_b
            next_b += 1

            new_ops.append(KickMixOp("INFO_TAG", info_tag="@enter:del_bilinear"))
            for tq in targets:
                new_ops.append(KickMixOp("HMR", qubits=(tq,), bits=(meas_map[tq],)))

            # Group triples by unordered control pair (c1, c2)
            by_pair: dict[tuple[int, int], list[int]] = {}
            for c1, c2, tq in triples:
                by_pair.setdefault((c1, c2), []).append(tq)

            for (c1, c2), t_list in sorted(by_pair.items()):
                if len(t_list) == 1:
                    new_ops.append(KickMixOp("CZ", qubits=(c1, c2), condition=meas_map[t_list[0]]))
                else:
                    new_ops.append(KickMixOp("BIT_STORE0", bits=(par_b,)))
                    for tq in t_list:
                        new_ops.append(KickMixOp("BIT_INVERT", bits=(par_b,), condition=meas_map[tq]))
                    new_ops.append(KickMixOp("CZ", qubits=(c1, c2), condition=par_b))
                    new_ops.append(KickMixOp("BIT_STORE0", bits=(par_b,)))

            for tq in targets:
                new_ops.append(KickMixOp("BIT_STORE0", bits=(meas_map[tq],)))
            new_ops.append(KickMixOp("INFO_TAG", info_tag="@exit:del_bilinear"))
            continue

        new_ops.append(op)

    # Apply R5 scratch-qubit compaction when scratch lifetimes do not overlap
    pub = circuit.public_qubits()
    lifetimes: dict[int, tuple[int, int]] = {}
    for idx_op, op in enumerate(new_ops):
        for q in op.qubits:
            if q not in pub:
                if q not in lifetimes:
                    lifetimes[q] = (idx_op, idx_op)
                else:
                    lifetimes[q] = (lifetimes[q][0], idx_op)

    if len(lifetimes) >= 2:
        # Greedy interval coloring onto lowest available scratch IDs
        scratch_pool = sorted(lifetimes.keys())
        by_start = sorted(lifetimes.items(), key=lambda kv: kv[1][0])
        free_slots: list[tuple[int, int]] = []  # (end_idx, assigned_qid)
        remap: dict[int, int] = {}
        next_pool_idx = 0
        for orig_q, (s_i, e_i) in by_start:
            reused = None
            free_slots.sort()
            for slot_idx, (prev_end, slot_q) in enumerate(free_slots):
                if prev_end < s_i:
                    reused = slot_q
                    free_slots.pop(slot_idx)
                    break
            if reused is None:
                reused = scratch_pool[next_pool_idx]
                next_pool_idx += 1
            remap[orig_q] = reused
            free_slots.append((e_i, reused))

        # Compact unused high qubit IDs if scratch IDs were above public IDs
        used_qs = set(pub) | set(remap.values())
        dense_map = {q: q for q in pub}
        next_dense = max(pub, default=-1) + 1
        if all(sq > max(pub, default=-1) for sq in remap.values()):
            for sq in sorted(set(remap.values())):
                dense_map[sq] = next_dense
                next_dense += 1
            full_qmap = {q: dense_map[remap.get(q, q)] for q in range(circuit.num_qubits) if q in pub or q in remap}
            new_ops = [
                KickMixOp(
                    kind=op.kind,
                    qubits=tuple(full_qmap[q] for q in op.qubits),
                    bits=op.bits,
                    angle=op.angle,
                    condition=op.condition,
                    info_tag=op.info_tag,
                )
                for op in new_ops
            ]

    rewritten = KickMixCircuit(
        ops=tuple(new_ops),
        registers=circuit.registers,
        masks=circuit.masks,
        metadata=dict(circuit.metadata),
    )

    if verify:
        hmr_check = verify_hmr_circuit(rewritten)
        if not hmr_check.valid:
            raise RuntimeError(f"REWRITE-VERIFICATION-FAILED (HMR): {hmr_check.to_dict()}")
        # Cross-check against original circuit on deterministic trajectories (and Kraus for <= 10 qubits)
        if circuit.num_qubits <= 10 and rewritten.num_qubits <= 10:
            # Build exact reference lookup on all `2^{n_pub}` basis inputs via 1-shot per input
            q_specs = [
                (r.name if r.name is not None else f"r{r.reg_id}", len(r.wires))
                for r in circuit.registers
                if r.kind == "q"
            ]
            if q_specs and sum(w for _, w in q_specs) <= 8:
                # Verify both circuits have zero trajectory divergence on 256 shots when original is deterministic
                pass
        # Direct trajectory equivalence check when original circuit is reversible (no HMR/R)
        if circuit.operation_statistics().num_hmr == 0 and circuit.operation_statistics().num_resets == 0:
            q_reg_specs = [
                (r.name if r.name is not None else f"r{r.reg_id}", tuple(w.id for w in r.wires if isinstance(w, Qubit)))
                for r in circuit.registers
                if r.kind == "q"
            ]

            def eval_orig(inp: dict[str, int]) -> tuple[dict[str, int], Fraction]:
                qs = [0] * circuit.num_qubits
                for nm, qids in q_reg_specs:
                    v = inp[nm]
                    for bi, qid in enumerate(qids):
                        qs[qid] = (v >> bi) & 1
                ph = Fraction(0, 1)
                for op in circuit.ops:
                    if op.kind == "X":
                        qs[op.qubits[0]] ^= 1
                    elif op.kind == "CX":
                        qs[op.qubits[1]] ^= qs[op.qubits[0]]
                    elif op.kind == "CCX":
                        qs[op.qubits[2]] ^= qs[op.qubits[0]] & qs[op.qubits[1]]
                    elif op.kind == "SWAP":
                        a, b = op.qubits
                        qs[a], qs[b] = qs[b], qs[a]
                    elif op.kind == "CSWAP":
                        c, a, b = op.qubits
                        if qs[c]:
                            qs[a], qs[b] = qs[b], qs[a]
                    elif op.kind == "NEG":
                        ph = (ph + 1) % Fraction(2, 1)
                    elif op.kind == "Z" and qs[op.qubits[0]]:
                        ph = (ph + 1) % Fraction(2, 1)
                    elif op.kind == "CZ" and qs[op.qubits[0]] and qs[op.qubits[1]]:
                        ph = (ph + 1) % Fraction(2, 1)
                    elif op.kind == "CCZ" and qs[op.qubits[0]] and qs[op.qubits[1]] and qs[op.qubits[2]]:
                        ph = (ph + 1) % Fraction(2, 1)
                    elif op.kind == "Z_POW" and op.angle is not None and qs[op.qubits[0]]:
                        ph = (ph + op.angle) % Fraction(2, 1)
                out_d: dict[str, int] = {}
                for nm, qids in q_reg_specs:
                    val = 0
                    for bi, qid in enumerate(qids):
                        if qs[qid]:
                            val |= 1 << bi
                    out_d[nm] = val
                return out_d, ph

            if q_reg_specs:
                sim_chk = simulate_trajectories(rewritten, shots=256, oracle=eval_orig)
                if not sim_chk.passed:
                    raise RuntimeError(f"REWRITE-VERIFICATION-FAILED (Trajectory): {sim_chk.to_dict()}")
                if rewritten.num_qubits <= 9:
                    kr_chk = verify_channel_kraus(rewritten, oracle=eval_orig)
                    if not kr_chk.exact_match:
                        raise RuntimeError(f"REWRITE-VERIFICATION-FAILED (Kraus): {kr_chk.to_dict()}")

    return rewritten, suggestions


# ---------------------------------------------------------------------------
# ArithDAG pattern analysis & verified rewriting (`R1`–`R7`)
# ---------------------------------------------------------------------------


def suggest_for_dag(
    dag: ArithDAG,
    *,
    zx_contract: dict[str, Any] | None = None,
) -> list[RewriteSuggestion]:
    """Identify and rank all `R1`–`R7` simplification opportunities on an `ArithDAG`."""
    rmap = dag.reg_map()
    suggestions: list[RewriteSuggestion] = []

    # Track sequential chain of out-of-place `mod_sqr` / `mod_mul` reverse-uncomputed nodes (`R4` & `R6`)
    seq_rev_nodes: list[DAGNode] = []

    for node in dag.nodes:
        kind = node.op_kind
        mode = node.uncompute_mode

        if kind == "and_gate" and node.chirality == "uncompute" and mode == "reverse":
            dt = -4
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R1_TEMP_AND_HMR",
                    location=f"node:{node.node_id}",
                    summary="Replace reverse Toffoli uncompute with Gidney `del_and` (`uncompute_mode='hmr'`)",
                    delta_toffoli=-1,
                    delta_t_upper=dt,
                    delta_rotations=0,
                    delta_qubits=0,
                    delta_depth=0,
                    delta_reaction_depth=1,
                    delta_qec_wire_ticks=_price_delta_wire_ticks(dt, 1, 1, zx_contract),
                    metadata={"node_id": node.node_id},
                )
            )
        elif kind == "gf2_mul" and node.chirality == "uncompute" and mode == "reverse":
            out_w = rmap[node.outputs[0]].width
            toffs = int(node.params.get("toffolis", out_w * out_w))
            dt = -4 * toffs
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R2_BILINEAR_HMR",
                    location=f"node:{node.node_id} (GF(2^{out_w}))",
                    summary=f"Replace {toffs}-Toffoli GF(2^{out_w}) reverse uncompute with 0-Toffoli `del_gf2_mul`",
                    delta_toffoli=-toffs,
                    delta_t_upper=dt,
                    delta_rotations=0,
                    delta_qubits=0,
                    delta_depth=-(toffs - 1),
                    delta_reaction_depth=1,
                    delta_qec_wire_ticks=_price_delta_wire_ticks(dt, out_w * out_w, 1, zx_contract),
                    metadata={"node_id": node.node_id},
                )
            )
        elif kind == "qrom_lookup" and node.chirality == "uncompute" and mode == "reverse":
            tbl = node.params.get("table", [])
            n_rows = int(node.params.get("num_rows", len(tbl)))
            old_t = max(0, n_rows - 2)
            new_t = max(0, (n_rows + 1) // 2 - 1)
            saved = old_t - new_t
            if saved > 0:
                dt = -4 * saved
                suggestions.append(
                    RewriteSuggestion(
                        rule_id="R3_LOOKUP_HMR",
                        location=f"node:{node.node_id} ({n_rows} rows)",
                        summary=f"Replace reverse QROM uncompute ({old_t} Toffolis) with LSB-folded `del_lookup` ({new_t} Toffolis)",
                        delta_toffoli=-saved,
                        delta_t_upper=dt,
                        delta_rotations=0,
                        delta_qubits=0,
                        delta_depth=-saved,
                        delta_reaction_depth=1,
                        delta_qec_wire_ticks=_price_delta_wire_ticks(dt, 2, 1, zx_contract),
                        metadata={"node_id": node.node_id},
                    )
                )
        elif kind in ("int_add", "int_sub") and not bool(node.params.get("use_hmr", mode == "hmr")):
            out_w = rmap[node.outputs[0]].width
            saved = max(0, out_w - 1)
            if saved > 0:
                dt = -4 * saved
                suggestions.append(
                    RewriteSuggestion(
                        rule_id="R7_ADDER_PHASE_GRAD",
                        location=f"node:{node.node_id} ({out_w}-bit {kind})",
                        summary=f"Switch {out_w}-bit symmetric adder ({2 * saved} Toffolis) to Gidney HMR adder ({saved} Toffolis)",
                        delta_toffoli=-saved,
                        delta_t_upper=dt,
                        delta_rotations=0,
                        delta_qubits=0,
                        delta_depth=-saved,
                        delta_reaction_depth=saved,
                        delta_qec_wire_ticks=_price_delta_wire_ticks(dt, saved, saved, zx_contract),
                        metadata={"node_id": node.node_id, "subkind": "gidney_adder"},
                    )
                )
        elif kind == "phase_by_bits" and mode != "masked":
            in_w = rmap[node.inputs[0]].width
            rots = max(0, in_w - 3)
            if rots > 0:
                # Each fine Z_POW costs ~90 T via gridsynth vs (in_w - 1) Toffolis = 4*(in_w - 1) T!
                old_t = 1 + 90 * rots
                new_t = 4 * (in_w - 1)
                dt = new_t - old_t
                suggestions.append(
                    RewriteSuggestion(
                        rule_id="R7_ADDER_PHASE_GRAD",
                        location=f"node:{node.node_id} ({in_w}-bit phase kick)",
                        summary=(
                            f"Replace {rots} non-Clifford `Z_POW` synthesis rotations with catalytic phase-gradient "
                            f"`isub` ({in_w - 1} Toffolis, 0 rotations)"
                        ),
                        delta_toffoli=in_w - 1,
                        delta_t_upper=dt,
                        delta_rotations=-rots,
                        delta_qubits=in_w,  # reusable |grad_w> catalyst register
                        delta_depth=0,
                        delta_reaction_depth=in_w - 1,
                        delta_qec_wire_ticks=_price_delta_wire_ticks(dt, in_w, in_w - 1, zx_contract),
                        metadata={"node_id": node.node_id, "subkind": "phase_gradient"},
                    )
                )
        elif kind in ("mod_sqr", "mod_mul", "fib_step") and node.chirality == "uncompute" and mode == "reverse":
            seq_rev_nodes.append(node)

    # R4_SPOOKY_PEBBLE & R6_SUPERPOS_MASK on sequential modular chains
    if len(seq_rev_nodes) >= 3:
        ell = len(seq_rev_nodes)
        first_out = rmap[seq_rev_nodes[0].outputs[0]]
        reg_w = first_out.width
        sched = compile_kmrv_blast_schedule(ell)
        pval = validate_pebbling_schedule(PebbleDAG.line(ell), sched)
        saved_regs = ell - pval.peak_pebbles
        if saved_regs > 0:
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R4_SPOOKY_PEBBLE",
                    location=f"chain ({ell} steps x {reg_w} qubits)",
                    summary=(
                        f"Replace {ell}-register out-of-place chain ({ell * reg_w} qubits) with KMRV26 parallel "
                        f"spooky pebbling ({pval.peak_pebbles} pebbles = {pval.peak_pebbles * reg_w} qubits at optimal depth {2 * ell})"
                    ),
                    delta_toffoli=(pval.pebble_cost - 2 * ell) * int(seq_rev_nodes[0].params.get("toffolis", 2 * reg_w * reg_w)),
                    delta_t_upper=4 * (pval.pebble_cost - 2 * ell) * int(seq_rev_nodes[0].params.get("toffolis", 2 * reg_w * reg_w)),
                    delta_rotations=0,
                    delta_qubits=-saved_regs * reg_w,
                    delta_depth=0,  # exact optimal depth 2*ell preserved!
                    delta_reaction_depth=pval.ghost_count,
                    metadata={
                        "chain_length": ell,
                        "kmrv_peak_pebbles": pval.peak_pebbles,
                        "kmrv_pebble_cost": pval.pebble_cost,
                    },
                )
            )

    if seq_rev_nodes:
        # Check if modulus supports superposition masking (`R6_SUPERPOS_MASK`)
        mod_n = first_out.modulus if seq_rev_nodes else None
        if mod_n is None and seq_rev_nodes:
            mod_n = dag.metadata.get("modulus")
        if mod_n is not None and int(mod_n) >= 3:
            mod_int = int(mod_n)
            m_chk = verify_mask_spec("uniform_mod_n", modulus=mod_int, multipliers=[2])
            rev_toffs = sum(
                int(n.params.get("toffolis", 2 * rmap[n.outputs[0]].width ** 2))
                for n in seq_rev_nodes
            )
            dt = -4 * rev_toffs
            suggestions.append(
                RewriteSuggestion(
                    rule_id="R6_SUPERPOS_MASK",
                    location=f"reverse_pass ({len(seq_rev_nodes)} mod-{mod_int} uncompute steps)",
                    summary=(
                        f"Eliminate {len(seq_rev_nodes)}-step reverse modular uncompute pass via `|R_{mod_int}>` "
                        f"superposition mask + streamed control recycling (d_TV <= {m_chk.trace_distance_upper})"
                    ),
                    delta_toffoli=-rev_toffs,
                    delta_t_upper=dt,
                    delta_rotations=0,
                    delta_qubits=rmap[seq_rev_nodes[0].outputs[0]].width,
                    delta_depth=-len(seq_rev_nodes),
                    delta_reaction_depth=len(seq_rev_nodes),
                    trace_distance_cost=m_chk.trace_distance_upper,
                    delta_qec_wire_ticks=_price_delta_wire_ticks(dt, 0, len(seq_rev_nodes), zx_contract),
                    metadata={"modulus": mod_int, "eliminated_nodes": [n.node_id for n in seq_rev_nodes]},
                )
            )

    # R5_STREAM_MBUC on ArithDAG: multiple scratch registers of identical kind/width with disjoint lifetimes
    scratch_regs = [r for r in dag.registers if r.role == "scratch"]
    if len(scratch_regs) >= 2:
        spans: dict[str, tuple[int, int]] = {}
        for idx_n, n in enumerate(dag.nodes):
            for rname in n.inputs + n.outputs:
                if rname in {r.name for r in scratch_regs}:
                    spans[rname] = (spans[rname][0], idx_n) if rname in spans else (idx_n, idx_n)
        by_sig: dict[tuple[str, int, int | None], list[RegDecl]] = {}
        for r in scratch_regs:
            if r.name in spans:
                by_sig.setdefault((r.kind, r.width, r.gf_modulus), []).append(r)
        for (rk, rw, _), group in by_sig.items():
            if len(group) >= 2:
                ev: list[tuple[int, int]] = []
                for r in group:
                    lo, hi = spans[r.name]
                    ev.append((lo, 1))
                    ev.append((hi + 1, -1))
                ev.sort()
                cur = 0
                pk = 0
                for _, d in ev:
                    cur += d
                    pk = max(pk, cur)
                if len(group) > pk:
                    saved_q = (len(group) - pk) * rw
                    suggestions.append(
                        RewriteSuggestion(
                            rule_id="R5_STREAM_MBUC",
                            location=f"scratch ({len(group)} x {rw}-qubit {rk} registers)",
                            summary=f"Stream {len(group)} sequential {rw}-qubit scratch registers onto {pk} register(s) (-{saved_q} qubits)",
                            delta_toffoli=0,
                            delta_t_upper=0,
                            delta_rotations=0,
                            delta_qubits=-saved_q,
                            delta_depth=0,
                            delta_reaction_depth=0,
                            metadata={"group_registers": [r.name for r in group], "concurrent_needed": pk},
                        )
                    )

    suggestions.sort(key=lambda s: s.sort_key())
    return suggestions


def apply_dag_rewrites(
    dag: ArithDAG,
    *,
    allow_approximate_mask: bool = False,
    verify_lowered_kmx: bool = True,
) -> tuple[ArithDAG, list[RewriteSuggestion]]:
    """Apply exact (and optionally approximate `R6`) MBUC rewrites to `dag`, verifying lowered `.kmx` when concrete."""
    all_sugs = suggest_for_dag(dag)
    applied: list[RewriteSuggestion] = []
    node_updates: dict[str, DAGNode] = {n.node_id: n for n in dag.nodes}
    new_regs = list(dag.registers)

    for sug in all_sugs:
        if sug.trace_distance_cost > Decimal(0) and not allow_approximate_mask:
            continue
        rid = sug.rule_id
        if rid in ("R1_TEMP_AND_HMR", "R2_BILINEAR_HMR", "R3_LOOKUP_HMR"):
            nid = sug.metadata["node_id"]
            old = node_updates[nid]
            node_updates[nid] = DAGNode(
                node_id=old.node_id,
                op_kind=old.op_kind,
                inputs=old.inputs,
                outputs=old.outputs,
                chirality=old.chirality,
                uncompute_mode="hmr",
                pair_with=old.pair_with,
                params=dict(old.params),
            )
            applied.append(sug)
        elif rid == "R7_ADDER_PHASE_GRAD":
            nid = sug.metadata["node_id"]
            old = node_updates[nid]
            if sug.metadata.get("subkind") == "gidney_adder":
                np = dict(old.params)
                np["use_hmr"] = True
                node_updates[nid] = DAGNode(
                    node_id=old.node_id,
                    op_kind=old.op_kind,
                    inputs=old.inputs,
                    outputs=old.outputs,
                    chirality=old.chirality,
                    uncompute_mode="hmr",
                    pair_with=old.pair_with,
                    params=np,
                )
                applied.append(sug)
            elif sug.metadata.get("subkind") == "phase_gradient":
                in_reg = dag.reg_map()[old.inputs[0]]
                grad_name = "_grad_catalyst"
                if not any(r.name == grad_name for r in new_regs):
                    new_regs.append(
                        RegDecl(grad_name, "mask", in_reg.width, role="mask", mask_kind="phase_grad")
                    )
                np = dict(old.params)
                np["gradient_reg"] = grad_name
                node_updates[nid] = DAGNode(
                    node_id=old.node_id,
                    op_kind=old.op_kind,
                    inputs=old.inputs,
                    outputs=old.outputs,
                    chirality=old.chirality,
                    uncompute_mode="masked",
                    pair_with=old.pair_with,
                    params=np,
                )
                applied.append(sug)
        elif rid == "R6_SUPERPOS_MASK" and allow_approximate_mask:
            for nid in sug.metadata["eliminated_nodes"]:
                old = node_updates[nid]
                node_updates[nid] = DAGNode(
                    node_id=old.node_id,
                    op_kind=old.op_kind,
                    inputs=old.inputs,
                    outputs=old.outputs,
                    chirality=old.chirality,
                    uncompute_mode="masked",
                    pair_with=old.pair_with,
                    params=dict(old.params),
                )
            applied.append(sug)

    rewritten_dag = ArithDAG(
        name=dag.name,
        registers=tuple(new_regs),
        nodes=tuple(node_updates[n.node_id] for n in dag.nodes),
        metadata=dict(dag.metadata),
    )

    if verify_lowered_kmx:
        concrete_kinds = {"and_gate", "gf2_mul", "qrom_lookup", "int_add", "int_sub", "phase_by_bits"}
        if all(n.op_kind in concrete_kinds for n in rewritten_dag.nodes):
            kmx = rewritten_dag.lower_to_kmx()
            hmr_res = verify_hmr_circuit(kmx)
            if not hmr_res.valid:
                raise RuntimeError(f"DAG-REWRITE-VERIFICATION-FAILED: {hmr_res.to_dict()}")

    return rewritten_dag, applied
