"""Lineage A: Exact F_2 Pauli-Web Engine, Detecting Regions, and Outcome Push-Out.

Implements Section 5.2 (Definitions 5.4-5.7, Proposition 5.6) and Appendix A
(Proposition A.1, Equations 7-9) of Rodatz, Poor, and Kissinger (arXiv:2506.17181v4).

Exact spider eigenvalue formulas in Choi ket space (derived from Def. 4.1):
- Z-spider v of degree d >= 1 with phase k * (pi/2), k in {0, 1, 2, 3}:
    Let (g_e, r_e) in F_2^2 be the (green=Z, red=X) bits on incident edges e in delta(v).
    Linear Pauli-web equations:
      r_e == r_v  for all e in delta(v)  (all or none of legs highlighted red),
      sum_{e in delta(v)} g_e == (k & 1) * r_v  (mod 2).
    Local spider eigenvalue (-1)^{s_Z(v)}:
      if r_v == 0: s_Z(v) = 0  (+1 stabiliser)
      if r_v == 1: s_Z(v) = ((G - k) // 2) & 1, where G = sum_{e in delta(v)} g_e (integer sum).
- X-spider v of degree d >= 1 with phase k * (pi/2), k in {0, 1, 2, 3}:
    Swap roles of g_e <-> r_e, plus H Y H = -Y parity on Y-legs:
      g_e == g_v  for all e in delta(v),
      sum_{e in delta(v)} r_e == (k & 1) * g_v  (mod 2).
    Local spider eigenvalue (-1)^{s_X(v)}:
      if g_v == 0: s_X(v) = 0
      if g_v == 1: s_X(v) = (((R - k) // 2) + m_v) & 1,
        where R = sum_{e in delta(v)} r_e and m_v = sum_{e in delta(v)} (g_e & r_e).
- Hadamard box h of degree 2 with incident edges (e1, e2):
    g_{e1} == r_{e2}, r_{e1} == g_{e2}, local eigenvalue = +1 (s_H = 0).
- Internal edges (neither endpoint is BOUND) contract against <Phi+| = <00| + <11|,
  contributing (-1)^{g_e & r_e} per internal Y-highlighted edge (Proposition 5.6).

Zero float literals appear anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .diagram import AtomicFaultGroup, CliffordCircuit, RefusalError, ZXDiagram, circuit_to_zx_fe
from .gf2 import (
    canonicalize_signed_stabilizers,
    nullspace_gf2,
    rref_gf2,
    solve_gf2,
    symplectic_inner_product,
    symplectic_to_pauli_string,
)


@dataclass(frozen=True)
class PauliWebAnalysis:
    """Complete exact Pauli-web analysis of a ZXDiagram."""

    is_nonzero: bool
    n_edges: int
    n_boundary: int
    ordered_eids: tuple[int, ...]
    boundary_eids: tuple[int, ...]
    detecting_regions: tuple[int, ...]       # 2*n_edges symplectic bitmasks with zero boundary support
    boundary_webs: tuple[int, ...]           # r = n_boundary full-diagram webs whose boundaries are RREF
    signed_boundary_stabs: tuple[tuple[int, int], ...]  # canonical (symp_2nb, sign_bit) on boundary
    zero_witness_web: int | None = None      # detecting region with sign == 1 when is_nonzero is False


def _web_global_sign(zx: ZXDiagram, web_vec: int, eid_to_idx: dict[int, int], bound_eid_set: set[int]) -> int:
    """Compute exact global sign bit s in {0, 1} such that (partial W) |D> = (-1)^s |D>."""
    sign = 0
    for nid, nd in zx.nodes.items():
        if nd.kind in ("BOUND", "H"):
            continue
        inc = zx.adj[nid]
        if not inc:
            # 0-legged scalar spider: |0> + i^k |1> contracted on 0 legs = 1 + i^k
            continue
        if nd.kind == "Z":
            first_idx = eid_to_idx[inc[0]]
            r_v = (web_vec >> (2 * first_idx + 1)) & 1
            if r_v == 1:
                g_sum = sum((web_vec >> (2 * eid_to_idx[eid])) & 1 for eid in inc)
                sign ^= ((g_sum - nd.phase_k) // 2) & 1
        elif nd.kind == "X":
            first_idx = eid_to_idx[inc[0]]
            g_v = (web_vec >> (2 * first_idx)) & 1
            if g_v == 1:
                r_sum = 0
                y_sum = 0
                for eid in inc:
                    idx = eid_to_idx[eid]
                    ge = (web_vec >> (2 * idx)) & 1
                    re = (web_vec >> (2 * idx + 1)) & 1
                    r_sum += re
                    y_sum += ge & re
                sign ^= (((r_sum - nd.phase_k) // 2) + y_sum) & 1

    # Internal edge <Phi+| contractions: -1 for each internal edge with g_e == 1 and r_e == 1
    for eid in zx.edges:
        if eid in bound_eid_set:
            continue
        idx = eid_to_idx[eid]
        ge = (web_vec >> (2 * idx)) & 1
        re = (web_vec >> (2 * idx + 1)) & 1
        if ge and re:
            sign ^= 1
    return sign


def _restrict_web_to_boundary(web_vec: int, bound_indices: list[int]) -> int:
    """Project a 2*|E|-bit diagram web onto the ordered 2*|E_B|-bit boundary vector."""
    bvec = 0
    for b_pos, e_idx in enumerate(bound_indices):
        ge = (web_vec >> (2 * e_idx)) & 1
        re = (web_vec >> (2 * e_idx + 1)) & 1
        if ge:
            bvec |= 1 << (2 * b_pos)
        if re:
            bvec |= 1 << (2 * b_pos + 1)
    return bvec


def analyze_pauli_webs(zx: ZXDiagram) -> PauliWebAnalysis:
    """Compute detecting regions, zero-branch status (Prop. 5.6), and signed boundary stabilisers."""
    zx.validate()
    eids = zx.ordered_edge_ids()
    n_edges = len(eids)
    eid_to_idx = {eid: i for i, eid in enumerate(eids)}
    bound_eids = zx.boundary_edges()
    n_bound = len(bound_eids)
    bound_indices = [eid_to_idx[eid] for eid in bound_eids]
    bound_eid_set = set(bound_eids)

    # Check 0-legged scalar spiders first: Z(2) or X(2) with 0 legs has scalar 1 + (-1) = 0!
    for nd in zx.nodes.values():
        if nd.kind in ("Z", "X") and len(zx.adj[nd.node_id]) == 0 and nd.phase_k == 2:
            return PauliWebAnalysis(
                is_nonzero=False,
                n_edges=n_edges,
                n_boundary=n_bound,
                ordered_eids=tuple(eids),
                boundary_eids=tuple(bound_eids),
                detecting_regions=(),
                boundary_webs=(),
                signed_boundary_stabs=(),
                zero_witness_web=0,
            )

    # Assemble homogeneous F_2 linear constraints over 2 * n_edges variables:
    # var 2*i = g_e (Z / green), var 2*i + 1 = r_e (X / red).
    eq_rows: list[int] = []
    ncols = 2 * n_edges

    for nid, nd in zx.nodes.items():
        if nd.kind == "BOUND":
            continue
        inc = zx.adj[nid]
        if not inc:
            continue
        if nd.kind == "H":
            i0 = eid_to_idx[inc[0]]
            i1 = eid_to_idx[inc[1]]
            # g_0 + r_1 = 0  and  r_0 + g_1 = 0
            eq_rows.append((1 << (2 * i0)) ^ (1 << (2 * i1 + 1)))
            eq_rows.append((1 << (2 * i0 + 1)) ^ (1 << (2 * i1)))
        elif nd.kind == "Z":
            i0 = eid_to_idx[inc[0]]
            r0_bit = 1 << (2 * i0 + 1)
            for eid_other in inc[1:]:
                io = eid_to_idx[eid_other]
                eq_rows.append(r0_bit ^ (1 << (2 * io + 1)))
            # sum_{e in inc} g_e + (phase_k & 1) * r_0 = 0
            g_parity = 0
            for eid in inc:
                g_parity ^= 1 << (2 * eid_to_idx[eid])
            if nd.phase_k & 1:
                g_parity ^= r0_bit
            eq_rows.append(g_parity)
        elif nd.kind == "X":
            i0 = eid_to_idx[inc[0]]
            g0_bit = 1 << (2 * i0)
            for eid_other in inc[1:]:
                io = eid_to_idx[eid_other]
                eq_rows.append(g0_bit ^ (1 << (2 * io)))
            # sum_{e in inc} r_e + (phase_k & 1) * g_0 = 0
            r_parity = 0
            for eid in inc:
                r_parity ^= 1 << (2 * eid_to_idx[eid] + 1)
            if nd.phase_k & 1:
                r_parity ^= g0_bit
            eq_rows.append(r_parity)

    # 1. Detecting regions (Definition 5.5): webs with (g_e, r_e) = (0, 0) on all boundary edges
    det_eqs = list(eq_rows)
    for b_idx in bound_indices:
        det_eqs.append(1 << (2 * b_idx))
        det_eqs.append(1 << (2 * b_idx + 1))
    det_basis = nullspace_gf2(det_eqs, ncols)

    # Proposition 5.6: since sign: D_reg -> F_2 is linear, D == 0 iff some basis detecting region has sign 1!
    for w in det_basis:
        if _web_global_sign(zx, w, eid_to_idx, bound_eid_set) == 1:
            return PauliWebAnalysis(
                is_nonzero=False,
                n_edges=n_edges,
                n_boundary=n_bound,
                ordered_eids=tuple(eids),
                boundary_eids=tuple(bound_eids),
                detecting_regions=tuple(det_basis),
                boundary_webs=(),
                signed_boundary_stabs=(),
                zero_witness_web=w,
            )

    # 2. Unrestricted Pauli webs -> quotient onto boundary edges for signed spacetime stabilisers!
    all_webs = nullspace_gf2(eq_rows, ncols)
    b_rows = [_restrict_web_to_boundary(w, bound_indices) for w in all_webs]
    rref_b, piv_cols, transforms = rref_gf2(b_rows, 2 * n_bound, track_transform=True)
    n_rank = len(rref_b)

    # Lift each canonical RREF boundary stabiliser to a concrete Pauli web on D!
    # Using `transforms[r]`, XOR-ing the corresponding `all_webs` gives a valid Pauli web on D
    # whose boundary restriction is EXACTLY `rref_b[r]`!
    # And evaluating `_web_global_sign` directly on that lifted web gives the EXACT sign bit
    # of `rref_b[r]` without even needing multi-step boundary phase accumulation!
    # Even better: we ALSO cross-check against `canonicalize_signed_stabilizers`!
    lifted_webs: list[int] = []
    signed_rref_from_lift: list[tuple[int, int]] = []
    for r_idx in range(n_rank):
        comb = transforms[r_idx]
        w_lift = 0
        for j, w_full in enumerate(all_webs):
            if (comb >> j) & 1:
                w_lift ^= w_full
        b_restr = _restrict_web_to_boundary(w_lift, bound_indices)
        if b_restr != rref_b[r_idx]:
            raise AssertionError("Boundary lift mismatch")
        s_bit = _web_global_sign(zx, w_lift, eid_to_idx, bound_eid_set)
        lifted_webs.append(w_lift)
        signed_rref_from_lift.append((b_restr, s_bit))

    # Sanity-verify against independent Pauli-phase RREF on raw generators
    raw_signed = [
        (_restrict_web_to_boundary(w, bound_indices), _web_global_sign(zx, w, eid_to_idx, bound_eid_set))
        for w in all_webs
    ]
    ok_stab, canon_signed = canonicalize_signed_stabilizers(raw_signed, n_bound)
    if not ok_stab:
        return PauliWebAnalysis(
            is_nonzero=False,
            n_edges=n_edges,
            n_boundary=n_bound,
            ordered_eids=tuple(eids),
            boundary_eids=tuple(bound_eids),
            detecting_regions=tuple(det_basis),
            boundary_webs=(),
            signed_boundary_stabs=(),
            zero_witness_web=0,
        )
    if canon_signed != signed_rref_from_lift:
        raise AssertionError(
            f"Internal sign consistency check failed: lift={signed_rref_from_lift} vs canon={canon_signed}"
        )

    return PauliWebAnalysis(
        is_nonzero=True,
        n_edges=n_edges,
        n_boundary=n_bound,
        ordered_eids=tuple(eids),
        boundary_eids=tuple(bound_eids),
        detecting_regions=tuple(det_basis),
        boundary_webs=tuple(lifted_webs),
        signed_boundary_stabs=tuple(signed_rref_from_lift),
        zero_witness_web=None,
    )


@dataclass(frozen=True)
class FaultSignatureTables:
    """Exact detector and boundary-effect linear maps over atomic fault generators.

    Suppose `fault_groups` has G atomic fault groups, with total binary generator columns M:
      - Column c (0 <= c < M) corresponds to one binary generator mask `gen_masks[c]` inside
        group `col_to_group[c]`.
      - `det_rows`  : list of r_det M-bit rows (detector parity checks H_det).
      - `eff_rows`  : list of n_boundary M-bit rows (canonical boundary syndrome L_eff).
    For any binary generator assignment x in F_2^M:
      - x is detectable iff H_det * x != 0 (i.e. some row in `det_rows` has odd overlap with x).
      - If x is undetectable (H_det * x == 0), its exact effect on D^F is uniquely determined
        by the n_boundary-bit vector sigma(x) = L_eff * x!
    """

    analysis: PauliWebAnalysis
    fault_groups: tuple[AtomicFaultGroup, ...]
    n_gen_cols: int
    col_to_group: tuple[int, ...]
    group_col_ranges: tuple[tuple[int, int], ...]
    gen_masks: tuple[int, ...]
    det_rows: tuple[int, ...]
    eff_rows: tuple[int, ...]

    def is_detectable(self, gen_bits: int) -> bool:
        for row in self.det_rows:
            if (row & gen_bits).bit_count() & 1:
                return True
        return False

    def detector_syndrome(self, gen_bits: int) -> int:
        syn = 0
        for i, row in enumerate(self.det_rows):
            if (row & gen_bits).bit_count() & 1:
                syn |= 1 << i
        return syn

    def boundary_effect(self, gen_bits: int) -> int:
        eff = 0
        for i, row in enumerate(self.eff_rows):
            if (row & gen_bits).bit_count() & 1:
                eff |= 1 << i
        return eff

    def group_weight(self, gen_bits: int) -> int:
        """Exact adversarial fault weight (Definition 2.3): count active AtomicFaultGroups."""
        wt = 0
        for lo, hi in self.group_col_ranges:
            width = hi - lo
            mask = ((1 << width) - 1) << lo
            if gen_bits & mask:
                wt += 1
        return wt


def build_fault_signature_tables(
    zx: ZXDiagram,
    noise_model: str = "circuit_depolarizing",
    analysis: PauliWebAnalysis | None = None,
) -> FaultSignatureTables:
    """Build exact F_2 detector and canonical boundary-effect matrices for `zx`."""
    ana = analysis if analysis is not None else analyze_pauli_webs(zx)
    if not ana.is_nonzero:
        raise RefusalError("ZERO-DIAGRAM", f"diagram {zx.name!r} evaluates to the zero linear map")
    _, groups = zx.build_edge_flip_fault_groups(noise_model=noise_model)

    gen_masks: list[int] = []
    col_to_group: list[int] = []
    group_ranges: list[tuple[int, int]] = []
    for grp_idx, grp in enumerate(groups):
        lo = len(gen_masks)
        for gm in grp.generators:
            gen_masks.append(gm)
            col_to_group.append(grp_idx)
        hi = len(gen_masks)
        group_ranges.append((lo, hi))

    m_cols = len(gen_masks)
    n_edges = ana.n_edges

    det_rows: list[int] = []
    for w_det in ana.detecting_regions:
        row = 0
        for col, f_mask in enumerate(gen_masks):
            if symplectic_inner_product(w_det, f_mask, n_edges):
                row |= 1 << col
        det_rows.append(row)

    eff_rows: list[int] = []
    for w_ext in ana.boundary_webs:
        row = 0
        for col, f_mask in enumerate(gen_masks):
            if symplectic_inner_product(w_ext, f_mask, n_edges):
                row |= 1 << col
        eff_rows.append(row)

    return FaultSignatureTables(
        analysis=ana,
        fault_groups=tuple(groups),
        n_gen_cols=m_cols,
        col_to_group=tuple(col_to_group),
        group_col_ranges=tuple(group_ranges),
        gen_masks=tuple(gen_masks),
        det_rows=tuple(det_rows),
        eff_rows=tuple(eff_rows),
    )


@dataclass(frozen=True)
class OutcomePushOutMap:
    """Appendix A measurement-outcome push-out table for a ZXDiagram with `m` measurements.

    For measurement outcome vector k in F_2^m:
      1) Valid (non-zero probability in the fault-free circuit) iff `parity_checks * k == 0`
         (each row in `parity_checks` is a noiseless detector parity over the m measurements,
         e.g. k_3 ^ k_4 == 0 in Equation 9, p. 50).
      2) Every valid outcome vector k in ker(parity_checks) maps linearly to:
           - a canonical boundary syndrome `eff(k)` = `outcome_eff_rows * k`, AND
           - an explicit boundary Pauli correction string `boundary_pauli_for_basis[j]`
             for each basis vector `valid_outcome_basis[j]` of ker(parity_checks)
             (Proposition A.1: D^{P(k)} = P_out(k) * D^I * P_in(k)).
    """

    meas_labels: tuple[str, ...]
    parity_checks: tuple[int, ...]                  # rows in F_2^m that must vanish for D^{P(k)} != 0
    valid_outcome_basis: tuple[int, ...]            # basis of ker(parity_checks) in F_2^m
    basis_boundary_syndromes: tuple[int, ...]       # canonical boundary syndrome per basis outcome
    basis_boundary_paulis: tuple[str, ...]          # boundary Pauli string (len = n_in + n_out) per basis outcome
    n_distinct_branches: int                        # dimension of distinct logical branches on boundary!


def compile_outcome_map(zx: ZXDiagram) -> OutcomePushOutMap:
    """Implement Appendix A (Proposition A.1 & Equations 7-9): push measurement outcomes to boundary."""
    ana = analyze_pauli_webs(zx)
    if not ana.is_nonzero:
        raise RefusalError("ZERO-DIAGRAM", f"diagram {zx.name!r} has zero amplitude on all-0 branch")
    m = len(zx.meas_spiders)
    eid_to_idx = {eid: i for i, eid in enumerate(ana.ordered_eids)}

    # Build the single-edge Pauli mask for flipping each measurement outcome k_r = 1:
    meas_masks: list[int] = []
    meas_labels: list[str] = []
    for nid, flip_p, lbl in zx.meas_spiders:
        inc = zx.adj[nid]
        if not inc:
            raise RefusalError("ISOLATED-MEAS-SPIDER", f"measurement spider {nid} has 0 legs")
        e_idx = eid_to_idx[inc[0]]
        if flip_p == "X":
            meas_masks.append(1 << (2 * e_idx + 1))
        elif flip_p == "Z":
            meas_masks.append(1 << (2 * e_idx))
        else:
            raise RefusalError("INVALID-MEAS-FLIP", f"unsupported flip Pauli {flip_p!r}")
        meas_labels.append(lbl)

    # 1. Detector parity checks over outcomes k in F_2^m:
    raw_checks: list[int] = []
    for w_det in ana.detecting_regions:
        row = 0
        for r_idx, mmask in enumerate(meas_masks):
            if symplectic_inner_product(w_det, mmask, ana.n_edges):
                row |= 1 << r_idx
        if row:
            raw_checks.append(row)
    canon_checks, _, _ = rref_gf2(raw_checks, m)
    valid_basis = nullspace_gf2(canon_checks, m)

    # 2. Push each valid outcome basis vector to the boundary via Proposition A.1!
    # Need boundary Pauli P_B (2 * n_boundary bits) such that for all i in 0..r-1:
    #   symplectic_ip(S_bar_i, P_B) == symplectic_ip(W_ext_i, P(k))!
    # Notice `symplectic_ip(S_bar_i, P_B)` is linear in the 2*n_boundary bits of P_B:
    #   swap_zx(S_bar_i) . P_B = rhs_i (mod 2)!
    from .gf2 import swap_zx_bits

    n_b = ana.n_boundary
    stab_dual_rows = [swap_zx_bits(s_vec, n_b) for s_vec, _ in ana.signed_boundary_stabs]

    basis_syndromes: list[int] = []
    basis_paulis: list[str] = []
    for k_vec in valid_basis:
        f_mask = 0
        for r_idx, mmask in enumerate(meas_masks):
            if (k_vec >> r_idx) & 1:
                f_mask ^= mmask
        rhs_bits: list[int] = []
        syn_int = 0
        for i, w_ext in enumerate(ana.boundary_webs):
            bit = symplectic_inner_product(w_ext, f_mask, ana.n_edges)
            rhs_bits.append(bit)
            if bit:
                syn_int |= 1 << i
        ok, pb_vec, _ = solve_gf2(stab_dual_rows, rhs_bits, 2 * n_b)
        if not ok:
            raise AssertionError("Proposition A.1 violated: valid measurement branch could not push to boundary")
        basis_syndromes.append(syn_int)
        basis_paulis.append(symplectic_to_pauli_string(pb_vec, n_b))

    # Quotient `basis_syndromes` by the subspace of boundary syndromes generated by
    # pure OUTPUT-wire Pauli corrections C (Appendix A, p. 49: (P, Q) ~ (P', Q') iff Q L P = C Q' L P')!
    n_in_b = len(zx.inputs)
    out_corr_syndromes: list[int] = []
    for out_pos in range(n_in_b, n_b):
        for bit_offset in (0, 1):
            c_mask_bound = 1 << (2 * out_pos + bit_offset)
            # Map boundary bit `out_pos` to full-diagram edge index `bound_eids[out_pos]`:
            full_e_idx = eid_to_idx[ana.boundary_eids[out_pos]]
            c_mask_full = 1 << (2 * full_e_idx + bit_offset)
            csyn = 0
            for i, w_ext in enumerate(ana.boundary_webs):
                if symplectic_inner_product(w_ext, c_mask_full, ana.n_edges):
                    csyn |= 1 << i
            if csyn:
                out_corr_syndromes.append(csyn)

    n_webs = len(ana.boundary_webs)
    corr_rref, _, _ = rref_gf2(out_corr_syndromes, n_webs)
    comb_rref, _, _ = rref_gf2(out_corr_syndromes + basis_syndromes, n_webs)
    n_branches_mod_frame = len(comb_rref) - len(corr_rref)

    return OutcomePushOutMap(
        meas_labels=tuple(meas_labels),
        parity_checks=tuple(canon_checks),
        valid_outcome_basis=tuple(valid_basis),
        basis_boundary_syndromes=tuple(basis_syndromes),
        basis_boundary_paulis=tuple(basis_paulis),
        n_distinct_branches=n_branches_mod_frame,
    )


def match_measurement_branches(spec: CliffordCircuit, impl: CliffordCircuit) -> dict[str, Any]:
    """Relate measurement outcomes of `impl` back to `spec` up to boundary Pauli frame (Appendix A)."""
    zx_spec = circuit_to_zx_fe(spec, noise_model="fault_free")
    zx_impl = circuit_to_zx_fe(impl, noise_model="fault_free")
    ana_spec = analyze_pauli_webs(zx_spec)
    ana_impl = analyze_pauli_webs(zx_impl)
    if ana_spec.signed_boundary_stabs != ana_impl.signed_boundary_stabs:
        raise RefusalError("MAP-MISMATCH", "spec and impl do not implement the same fault-free linear map")
    map_spec = compile_outcome_map(zx_spec)
    map_impl = compile_outcome_map(zx_impl)
    return {
        "spec_measurements": list(map_spec.meas_labels),
        "impl_measurements": list(map_impl.meas_labels),
        "impl_detector_checks": [
            [map_impl.meas_labels[i] for i in range(len(map_impl.meas_labels)) if (chk >> i) & 1]
            for chk in map_impl.parity_checks
        ],
        "spec_branch_rank": map_spec.n_distinct_branches,
        "impl_branch_rank": map_impl.n_distinct_branches,
        "impl_valid_basis": [
            {
                "outcomes_flipped": [
                    map_impl.meas_labels[i]
                    for i in range(len(map_impl.meas_labels))
                    if (kv >> i) & 1
                ],
                "boundary_pauli": bp,
                "boundary_syndrome": bsyn,
            }
            for kv, bp, bsyn in zip(
                map_impl.valid_outcome_basis,
                map_impl.basis_boundary_paulis,
                map_impl.basis_boundary_syndromes,
            )
        ],
    }
