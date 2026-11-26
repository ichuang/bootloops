"""Constructive FTQC Gadget Synthesizers with Machine-Checkable ZX Proof Certificates.

Implements Section 7 of Rodatz, Poor, and Kissinger (arXiv:2506.17181v4):
  - Cat state preparation:
      * Unflagged 4-cat negative control (§3.1, p. 10)
      * Flagged 4-cat (Proposition 7.1, p. 29)
      * Arbitrary-distance recursive 2n-cat (Proposition 7.2, p. 30)
      * Truncated w-fault-tolerant 2n-cat using 2w CNOTs & w checks (Proposition 7.3, p. 30)
  - Shor-style syndrome measurement & alternatives:
      * Standard Shor Pauli measurement (Proposition 7.5, p. 33)
      * Ancilla-reduced optimised Shor via H|cat>-like map (Proposition 7.6, p. 33)
      * Teleportation-routing syndrome measurement (Proposition 7.7, p. 35)
      * d-round Shor syndrome extraction schedule (Proposition 7.4, p. 32)
  - Steane-style syndrome extraction & leg-bending optimisation:
      * Goto 8-qubit 3-fault-tolerant Steane |0_L> preparation (Proposition 7.10, p. 38)
      * Bent-leg 5-qubit |0>-like linear map (Proposition 7.11, p. 39)
      * Standard 8-ancilla (18-CNOT) vs Optimised 5-ancilla (15-CNOT) Steane extraction
        (Propositions 7.8, 7.9, 7.12, pp. 36-40)
      * Naive 1-round sequential cat_4 Steane negative control (hook-error counterexample).

Zero float literals appear anywhere in this module.
"""

from __future__ import annotations

from .diagram import CircuitOp, CliffordCircuit, RefusalError, ZXDiagram, circuit_to_zx_fe
from .equivalence import CSSCodeSpec, make_steane_7_1_3
from .rewrites import DerivationCertificate, RewriteStep


def make_cat_spec(n_qubits: int, basis: str = "Z") -> CliffordCircuit:
    """Idealised fault-free |cat_n> preparation specification (Definition 3.2)."""
    if n_qubits < 1:
        raise RefusalError("INVALID-CAT-SIZE", f"n_qubits must be >= 1, got {n_qubits}")
    if basis not in ("Z", "X"):
        raise RefusalError("INVALID-CAT-BASIS", f"basis must be 'Z' or 'X', got {basis!r}")
    qs = tuple(range(n_qubits))
    return CliffordCircuit(
        n_qubits=n_qubits,
        input_qubits=(),
        output_qubits=qs,
        ops=[CircuitOp("CAT_PREP", qs, mode="ff", params={"basis": basis})],
        name=f"ideal_cat_{basis}_{n_qubits}",
    )


def synth_unflagged_cat_state(n_qubits: int = 4, basis: str = "Z") -> CliffordCircuit:
    """Non-fault-tolerant cat-state circuit (§3.1 p. 10 counterexample for n >= 4)."""
    if n_qubits < 1:
        raise RefusalError("INVALID-CAT-SIZE", f"n_qubits={n_qubits}")
    ops: list[CircuitOp] = [CircuitOp("PREP_PLUS", (0,))]
    for q in range(1, n_qubits):
        ops.append(CircuitOp("PREP_0", (q,)))
    for q in range(n_qubits - 1):
        ops.append(CircuitOp("CX", (q, q + 1)))
    if basis == "X":
        for q in range(n_qubits):
            ops.append(CircuitOp("H", (q,)))
    return CliffordCircuit(
        n_qubits=n_qubits,
        input_qubits=(),
        output_qubits=tuple(range(n_qubits)),
        ops=ops,
        name=f"unflagged_cat_{basis}_{n_qubits}",
    )


def synth_cat_state(
    n_qubits: int,
    w: int | None = None,
    basis: str = "Z",
) -> tuple[CliffordCircuit, DerivationCertificate]:
    """Constructive fault-tolerant cat-state compiler (Propositions 7.1, 7.2, 7.3).

    - For n_qubits <= 3: unflagged tree is already fault-equivalent for all w (any weight-2
      X error on <=3 qubits equals a weight-1 X error modulo X^{otimes n}).
    - For n_qubits == 4 and (w is None or w >= 2): Proposition 7.1 (1 flag ancilla, 5 CNOTs),
      or if w == 1: unflagged tree.
    - For even n_qubits = 2*m >= 6: builds two m-qubit FT cat states and measures
      `m_checks = m if w is None else min(w, m)` transversal ZZ parity checks into |0>
      flag ancillas (Proposition 7.2 when w is None, Proposition 7.3 when w < m!).
    - For odd n_qubits >= 5: unfuses 1 leg from an (n_qubits - 1)-qubit FT cat state via
      `|0>; CX(last, new)` (`unfuse_1_fe` + colour-dual Prop. 6.10!).
    """
    if n_qubits < 1:
        raise RefusalError("INVALID-CAT-SIZE", f"n_qubits must be >= 1, got {n_qubits}")
    if w is not None and w < 1:
        raise RefusalError("INVALID-WEIGHT-BOUND", f"w must be >= 1, got {w}")
    if basis not in ("Z", "X"):
        raise RefusalError("INVALID-CAT-BASIS", f"basis must be 'Z' or 'X', got {basis!r}")

    spec = make_cat_spec(n_qubits, basis=basis)

    def build_z_cat(n: int, w_target: int | None) -> tuple[CliffordCircuit, DerivationCertificate]:
        z_spec = make_cat_spec(n, basis="Z")
        init_zx = circuit_to_zx_fe(z_spec, "circuit_depolarizing")

        if n <= 3 or w_target == 1:
            circ = synth_unflagged_cat_state(n, basis="Z")
            cert = DerivationCertificate(
                name=f"cat_z_{n}_base",
                spec_circuit=z_spec,
                initial_zx=init_zx,
                steps=(),
                extracted_circuit=circ,
                target_noise_model="circuit_depolarizing",
                declared_w=(1 if (w_target == 1 and n >= 4) else None),
            )
            return circ, cert

        if n == 4:
            # Proposition 7.1 (p. 29): data 0..3, flag 4!
            ops = [
                CircuitOp("PREP_PLUS", (0,)),
                CircuitOp("PREP_0", (1,)),
                CircuitOp("PREP_0", (2,)),
                CircuitOp("PREP_0", (3,)),
                CircuitOp("PREP_0", (4,)),
                CircuitOp("CX", (0, 4)),
                CircuitOp("CX", (0, 1)),
                CircuitOp("CX", (1, 2)),
                CircuitOp("CX", (2, 3)),
                CircuitOp("CX", (3, 4)),
                CircuitOp("MEAS_Z", (4,), label="cat4_flag"),
            ]
            circ = CliffordCircuit(5, (), (0, 1, 2, 3), ops, name="cat4_flagged_prop_7_1")
            # Explicit ZX rewrite step: apply Fuse-4_fe on the central 4-legged Z-spider!
            hub_id = [nid for nid, nd in init_zx.nodes.items() if nd.kind == "Z"][0]
            steps = (RewriteStep("fuse_4_fe", {"node_id": hub_id}, note="Prop 6.14 detecting 4-cycle"),)
            cert = DerivationCertificate(
                name="cat4_prop_7_1",
                spec_circuit=z_spec,
                initial_zx=init_zx,
                steps=steps,
                extracted_circuit=circ,
                target_noise_model="circuit_depolarizing",
                declared_w=None,
            )
            return circ, cert

        if n % 2 == 1:
            # Odd n >= 5: prepare an even (n + 1)-cat state and project leg n onto <+| via MEAS_X!
            # By Proposition 6.13 (unfuse_1_fe with a 0-external-leg Z(0) leaf), terminating one
            # leg of |cat_{n+1}> in <+| (MEAS_X) preserves exact w-fault equivalence for ALL w!
            sub_circ, sub_cert = build_z_cat(n + 1, w_target)
            ops = list(sub_circ.ops) + [CircuitOp("MEAS_X", (n,), label=f"cat{n}_trim_x")]
            circ = CliffordCircuit(
                sub_circ.n_qubits,
                (),
                tuple(range(n)),
                ops,
                name=f"cat_{n}_odd",
            )
            hub_id = [nid for nid, nd in init_zx.nodes.items() if nd.kind == "Z"][0]
            steps = (
                RewriteStep(
                    "unfuse_1_fe",
                    {"node_id": hub_id, "move_edge_ids": [], "new_phase_k": 0},
                    note="Prop 6.13 unfuse 1-legged <+| (MEAS_X) leaf from (n+1)-cat",
                ),
            )
            cert = DerivationCertificate(
                name=f"cat_{n}_odd_cert",
                spec_circuit=z_spec,
                initial_zx=init_zx,
                steps=steps,
                extracted_circuit=circ,
                target_noise_model="circuit_depolarizing",
                declared_w=sub_cert.declared_w,
                compositional_subproofs=(sub_cert,),
            )
            return circ, cert

        # Even n = 2 * m >= 6: Propositions 7.2 (w_target is None or >= m) & 7.3 (w_target < m)!
        m = n // 2
        sub_c1, sub_cert1 = build_z_cat(m, w_target)
        sub_c2, sub_cert2 = build_z_cat(m, w_target)
        n_anc1 = sub_c1.n_qubits - m
        n_anc2 = sub_c2.n_qubits - m
        # Layout qubits cleanly so outputs are 0 .. n-1:
        #   top data: 0 .. m-1
        #   bot data: m .. 2m-1
        #   top internal ancillas: 2m .. 2m + n_anc1 - 1
        #   bot internal ancillas: 2m + n_anc1 .. 2m + n_anc1 + n_anc2 - 1
        #   parity-bridge ancillas: start after that!
        map1 = {q: (q if q < m else 2 * m + (q - m)) for q in range(sub_c1.n_qubits)}
        map2 = {q: (m + q if q < m else 2 * m + n_anc1 + (q - m)) for q in range(sub_c2.n_qubits)}
        ops = []
        for op in sub_c1.ops:
            ops.append(CircuitOp(op.kind, tuple(map1[q] for q in op.qubits), op.mode, op.label, dict(op.params)))
        for op in sub_c2.ops:
            ops.append(CircuitOp(op.kind, tuple(map2[q] for q in op.qubits), op.mode, op.label, dict(op.params)))

        m_checks = m if (w_target is None or w_target >= m) else w_target
        bridge_base = 2 * m + n_anc1 + n_anc2
        for j in range(m_checks):
            fq = bridge_base + j
            ops.append(CircuitOp("PREP_0", (fq,)))
            ops.append(CircuitOp("CX", (j, fq)))
            ops.append(CircuitOp("CX", (m + j, fq)))
            ops.append(CircuitOp("MEAS_Z", (fq,), label=f"cat{n}_zz_{j}"))

        total_q = bridge_base + m_checks
        circ = CliffordCircuit(
            total_q,
            (),
            tuple(range(n)),
            ops,
            name=f"cat_{n}_w{w_target if w_target is not None else 'inf'}",
        )
        hub_id = [nid for nid, nd in init_zx.nodes.items() if nd.kind == "Z"][0]
        if m_checks == m:
            step = RewriteStep("fuse_n_fe", {"node_id": hub_id}, note=f"Prop 6.15 ({m} bridges)")
            decl_w = sub_cert1.declared_w
        else:
            step = RewriteStep("fuse_n_w_fe", {"node_id": hub_id, "w": m_checks}, note=f"Prop 6.16 ({m_checks} bridges)")
            decl_w = m_checks if sub_cert1.declared_w is None else min(m_checks, sub_cert1.declared_w)

        cert = DerivationCertificate(
            name=f"cat_{n}_prop_{'7_2' if m_checks == m else '7_3'}",
            spec_circuit=z_spec,
            initial_zx=init_zx,
            steps=(step,),
            extracted_circuit=circ,
            target_noise_model="circuit_depolarizing",
            declared_w=decl_w,
            compositional_subproofs=(sub_cert1, sub_cert2),
        )
        return circ, cert

    z_circ, z_cert = build_z_cat(n_qubits, w)
    if basis == "Z":
        return z_circ, z_cert

    # X-basis cat state (H |cat>): apply transversal Hadamard on all n_qubits outputs (Corollary 6.6)!
    x_ops = list(z_circ.ops) + [CircuitOp("H", (q,)) for q in range(n_qubits)]
    x_circ = CliffordCircuit(
        z_circ.n_qubits,
        (),
        z_circ.output_qubits,
        x_ops,
        name=f"H_{z_circ.name}",
    )
    x_init_zx = circuit_to_zx_fe(spec, "circuit_depolarizing")
    x_cert = DerivationCertificate(
        name=f"H_{z_cert.name}",
        spec_circuit=spec,
        initial_zx=x_init_zx,
        steps=(),
        extracted_circuit=x_circ,
        target_noise_model="circuit_depolarizing",
        declared_w=z_cert.declared_w,
        compositional_subproofs=(z_cert,),
    )
    return x_circ, x_cert


def make_syndrome_meas_spec(pauli_str: str, mode: str = "ft") -> CliffordCircuit:
    """Idealised fault-tolerant ('ft') or fault-free ('ff') Pauli measurement spec (Definition 3.5)."""
    n = len(pauli_str)
    if n < 1 or any(ch not in ("I", "X", "Y", "Z") for ch in pauli_str):
        raise RefusalError("INVALID-PAULI-STRING", f"invalid Pauli word {pauli_str!r}")
    qs = tuple(range(n))
    return CliffordCircuit(
        n_qubits=n,
        input_qubits=qs,
        output_qubits=qs,
        ops=[CircuitOp("MPP", qs, mode=mode, label=f"spec_{pauli_str}", params={"pauli": pauli_str})],
        name=f"spec_meas_{mode}_{pauli_str}",
    )


def _basis_pre_ops(active_qubits: list[int], active_paulis: list[str]) -> tuple[list[CircuitOp], list[CircuitOp]]:
    """Single-qubit Cliffords mapping active_paulis -> Z and back (Corollary 6.6)."""
    pre: list[CircuitOp] = []
    post: list[CircuitOp] = []
    for q, pch in zip(active_qubits, active_paulis):
        if pch == "Z":
            continue
        if pch == "X":
            pre.append(CircuitOp("H", (q,)))
            post.append(CircuitOp("H", (q,)))
        elif pch == "Y":
            pre.append(CircuitOp("SDG", (q,)))
            pre.append(CircuitOp("H", (q,)))
            post.append(CircuitOp("H", (q,)))
            post.append(CircuitOp("S", (q,)))
    return pre, post


def synth_shor_syndrome(
    pauli_str: str,
    variant: str = "optimised",
) -> tuple[CliffordCircuit, DerivationCertificate]:
    """Synthesize a fault-tolerant Pauli syndrome measurement circuit (Props. 7.5, 7.6, 7.7).

    Variants:
      - 'standard' : Proposition 7.5 — fault-tolerant H|cat_r> prep + r CNOTs + r MEAS_Z.
      - 'optimised': Proposition 7.6 — ancilla-halved H|cat>-like linear map using 2 |+> -> CX -> MEAS_X
                     bridges per pair (on weight-4 ZZZZ: 3 concurrent ancillas instead of 5!).
      - 'teleport' : Proposition 7.7 — dynamic teleportation-routing syndrome measurement
                     (supported on weight-4 Pauli checks).
    """
    n_data = len(pauli_str)
    active = [(q, ch) for q, ch in enumerate(pauli_str) if ch != "I"]
    r = len(active)
    if r < 2:
        raise RefusalError("TRIVIAL-SYNDROME-WEIGHT", f"Pauli word {pauli_str!r} has active weight {r} < 2")
    act_qs = [q for q, _ in active]
    act_ps = [ch for _, ch in active]
    pre_ops, post_ops = _basis_pre_ops(act_qs, act_ps)
    spec = make_syndrome_meas_spec(pauli_str, mode="ft")
    init_zx = circuit_to_zx_fe(spec, "circuit_depolarizing")

    if variant == "standard":
        hcat_circ, hcat_cert = synth_cat_state(r, w=None, basis="X")
        remap = {q: n_data + q for q in range(hcat_circ.n_qubits)}
        ops: list[CircuitOp] = list(pre_ops)
        for op in hcat_circ.ops:
            ops.append(CircuitOp(op.kind, tuple(remap[q] for q in op.qubits), op.mode, op.label, dict(op.params)))
        for idx_a, dq in enumerate(act_qs):
            ops.append(CircuitOp("CX", (dq, n_data + idx_a)))
        for idx_a in range(r):
            ops.append(CircuitOp("MEAS_Z", (n_data + idx_a,), label=f"shor_m{idx_a}"))
        ops.extend(post_ops)
        circ = CliffordCircuit(
            n_data + hcat_circ.n_qubits,
            tuple(range(n_data)),
            tuple(range(n_data)),
            ops,
            name=f"shor_standard_{pauli_str}",
        )
        cert = DerivationCertificate(
            name=f"shor_standard_{pauli_str}_prop_7_5",
            spec_circuit=spec,
            initial_zx=init_zx,
            steps=(),
            extracted_circuit=circ,
            target_noise_model="circuit_depolarizing",
            declared_w=None,
            compositional_subproofs=(hcat_cert,),
        )
        return circ, cert

    if variant == "optimised":
        # Proposition 7.6 generalised to arbitrary r >= 2:
        # Split active qubits into two halves: left half (m_left = r // 2) and right half (r - m_left).
        # For r == 4 (exact Prop. 7.6, p. 33):
        #   ancillas a0 = n_data (|0>), a1 = n_data+1 (|0>), bridge1 = n_data+2 (|+>..MEAS_X), bridge2 = n_data+3 (|+>..MEAS_X).
        #   Note bridge1 and bridge2 run sequentially, so peak physical ancilla count is 3!
        if r == 4:
            a0, a1, br0, br1 = n_data, n_data + 1, n_data + 2, n_data + 3
            q0, q1, q2, q3 = act_qs
            ops = list(pre_ops) + [
                CircuitOp("PREP_0", (a0,)),
                CircuitOp("PREP_0", (a1,)),
                CircuitOp("CX", (q0, a0)),
                CircuitOp("CX", (q1, a1)),
                CircuitOp("PREP_PLUS", (br0,)),
                CircuitOp("CX", (br0, a0)),
                CircuitOp("CX", (br0, a1)),
                CircuitOp("MEAS_X", (br0,), label="k3_flagA"),
                CircuitOp("PREP_PLUS", (br1,)),
                CircuitOp("CX", (br1, a0)),
                CircuitOp("CX", (br1, a1)),
                CircuitOp("MEAS_X", (br1,), label="k4_flagB"),
                CircuitOp("CX", (q2, a0)),
                CircuitOp("CX", (q3, a1)),
                CircuitOp("MEAS_Z", (a0,), label="k1"),
                CircuitOp("MEAS_Z", (a1,), label="k2"),
            ] + list(post_ops)
            circ = CliffordCircuit(
                n_data + 4,
                tuple(range(n_data)),
                tuple(range(n_data)),
                ops,
                name=f"shor_optimised_{pauli_str}_prop_7_6",
            )
            cert = DerivationCertificate(
                name=f"shor_optimised_{pauli_str}_prop_7_6",
                spec_circuit=spec,
                initial_zx=init_zx,
                steps=(),
                extracted_circuit=circ,
                target_noise_model="circuit_depolarizing",
                declared_w=None,
            )
            return circ, cert

        # General r != 4: pair up floor(r/2) ancilla rails with sequential double-XX flag bridges!
        n_rails = (r + 1) // 2
        rails = [n_data + i for i in range(n_rails)]
        ops = list(pre_ops)
        for rq in rails:
            ops.append(CircuitOp("PREP_0", (rq,)))
        for i in range(n_rails):
            ops.append(CircuitOp("CX", (act_qs[i], rails[i])))
        next_q = n_data + n_rails
        # Connect adjacent rails (rails[i], rails[i+1]) with double |+> -> CX -> MEAS_X bridges
        # so the n_rails rails form a fault-equivalent H|cat>-like hub!
        for i in range(n_rails - 1):
            for rep in range(2):
                bq = next_q
                next_q += 1
                ops.append(CircuitOp("PREP_PLUS", (bq,)))
                ops.append(CircuitOp("CX", (bq, rails[i])))
                ops.append(CircuitOp("CX", (bq, rails[i + 1])))
                ops.append(CircuitOp("MEAS_X", (bq,), label=f"flag_{i}_{rep}"))
        for j in range(n_rails, r):
            ops.append(CircuitOp("CX", (act_qs[j], rails[j - n_rails])))
        for i, rq in enumerate(rails):
            ops.append(CircuitOp("MEAS_Z", (rq,), label=f"shor_opt_m{i}"))
        ops.extend(post_ops)
        circ = CliffordCircuit(
            next_q,
            tuple(range(n_data)),
            tuple(range(n_data)),
            ops,
            name=f"shor_optimised_{pauli_str}",
        )
        cert = DerivationCertificate(
            name=f"shor_optimised_{pauli_str}",
            spec_circuit=spec,
            initial_zx=init_zx,
            steps=(),
            extracted_circuit=circ,
            target_noise_model="circuit_depolarizing",
            declared_w=None,
        )
        return circ, cert

    if variant == "teleport":
        if r != 4:
            raise RefusalError(
                "TELEPORT-SYNDROME-ARITY",
                f"Proposition 7.7 teleportation syndrome gadget requires active weight 4, got {r}",
            )
        q0, q1, q2, q3 = act_qs
        a0, a1, br0, br1 = n_data, n_data + 1, n_data + 2, n_data + 3
        # Proposition 7.7 (p. 35): data shifts down by two rails!
        #   in q0 -> out q2, in q1 -> out q3, in q2 -> out a0, in q3 -> out a1!
        mig = {q0: q2, q1: q3, q2: a0, q3: a1}
        out_qs = [mig.get(q, q) for q in range(n_data)]
        post_migrated: list[CircuitOp] = []
        for q_orig, pch in zip(act_qs, act_ps):
            qt = mig[q_orig]
            if pch == "X":
                post_migrated.append(CircuitOp("H", (qt,)))
            elif pch == "Y":
                post_migrated.append(CircuitOp("H", (qt,)))
                post_migrated.append(CircuitOp("S", (qt,)))
        ops = list(pre_ops) + [
            CircuitOp("PREP_0", (a0,)),
            CircuitOp("PREP_0", (a1,)),
            CircuitOp("CX", (q2, a0)),
            CircuitOp("CX", (q3, a1)),
            CircuitOp("PREP_PLUS", (br0,)),
            CircuitOp("CX", (br0, q2)),
            CircuitOp("CX", (br0, q3)),
            CircuitOp("MEAS_X", (br0,), label="tp_flagA"),
            CircuitOp("PREP_PLUS", (br1,)),
            CircuitOp("CX", (br1, q2)),
            CircuitOp("CX", (br1, q3)),
            CircuitOp("MEAS_X", (br1,), label="tp_flagB"),
            CircuitOp("CX", (q2, q0)),
            CircuitOp("CX", (q3, q1)),
            CircuitOp("MEAS_Z", (q0,), label="tp_m1"),
            CircuitOp("MEAS_Z", (q1,), label="tp_m2"),
        ] + post_migrated
        circ = CliffordCircuit(
            n_data + 4,
            tuple(range(n_data)),
            tuple(out_qs),
            ops,
            name=f"shor_teleport_{pauli_str}_prop_7_7",
        )
        cert = DerivationCertificate(
            name=f"shor_teleport_{pauli_str}_prop_7_7",
            spec_circuit=spec,
            initial_zx=init_zx,
            steps=(),
            extracted_circuit=circ,
            target_noise_model="circuit_depolarizing",
            declared_w=None,
        )
        return circ, cert

    raise RefusalError("UNKNOWN-SHOR-VARIANT", f"variant must be standard, optimised, or teleport, got {variant!r}")


def make_full_syndrome_spec(code: CSSCodeSpec, sector: str = "X") -> CliffordCircuit:
    """Idealised fault-free syndrome extraction specification (Definition 3.4, p. 12).

    `sector == 'X'` measures the X-stabilisers (detecting Z data errors, matching Props. 7.8-7.12).
    `sector == 'Z'` measures the Z-stabilisers (detecting X data errors).
    `sector == 'BOTH'` measures all X- and Z-stabilisers fault-freely without intermediate wire flips.
    """
    qs = tuple(range(code.n))
    rows_list: list[tuple[str, int]] = []
    if sector in ("X", "BOTH"):
        for r in code.hx_rows:
            rows_list.append(("X", r))
    if sector in ("Z", "BOTH"):
        for r in code.hz_rows:
            rows_list.append(("Z", r))
    if not rows_list:
        raise RefusalError("INVALID-SYNDROME-SECTOR", f"sector={sector!r}")

    ops: list[CircuitOp] = []
    ff_wires: set[tuple[int, int]] = set()
    seg_count = {q: 0 for q in qs}
    for m_idx, (pch, mask) in enumerate(rows_list):
        pstr = "".join(pch if (mask >> q) & 1 else "I" for q in qs)
        if m_idx > 0:
            # Definition 3.4: NO faults in between the idealised stabiliser measurements!
            for q in qs:
                ff_wires.add((q, seg_count[q]))
        ops.append(CircuitOp("MPP", qs, mode="ff", label=f"spec_S_{pch}_{m_idx}", params={"pauli": pstr}))
        for q in qs:
            seg_count[q] += 1

    return CliffordCircuit(
        n_qubits=code.n,
        input_qubits=qs,
        output_qubits=qs,
        ops=ops,
        fault_free_wires=frozenset(ff_wires),
        name=f"spec_syndrome_{code.name}_{sector}",
    )


def synth_steane_goto_prep() -> CliffordCircuit:
    """Goto 3-fault-tolerant 8-qubit |0_L> preparation for Steane [[7, 1, 3]] (Prop. 7.10)."""
    ops = [
        CircuitOp("PREP_PLUS", (0,)),
        CircuitOp("PREP_0", (1,)),
        CircuitOp("PREP_0", (2,)),
        CircuitOp("PREP_0", (3,)),
        CircuitOp("PREP_PLUS", (4,)),
        CircuitOp("PREP_0", (5,)),
        CircuitOp("PREP_PLUS", (6,)),
        CircuitOp("PREP_0", (7,)),
        CircuitOp("CX", (4, 1)),
        CircuitOp("CX", (6, 2)),
        CircuitOp("CX", (0, 1)),
        CircuitOp("CX", (4, 5)),
        CircuitOp("CX", (6, 5)),
        CircuitOp("CX", (1, 2)),
        CircuitOp("CX", (0, 3)),
        CircuitOp("CX", (6, 3)),
        CircuitOp("CX", (1, 7)),
        CircuitOp("CX", (3, 7)),
        CircuitOp("CX", (5, 7)),
        CircuitOp("MEAS_Z", (7,), label="goto_flag"),
    ]
    return CliffordCircuit(8, (), (0, 1, 2, 3, 4, 5, 6), ops, name="steane_goto_zero_prep")


def synth_steane_zero_like() -> CliffordCircuit:
    """3-fault-tolerant 5-qubit |0>-like linear map for Steane [[7, 1, 3]] (Prop. 7.11)."""
    ops = [
        CircuitOp("PREP_0", (1,)),
        CircuitOp("PREP_0", (4,)),
        CircuitOp("CX", (2, 1)),
        CircuitOp("CX", (2, 3)),
        CircuitOp("CX", (3, 4)),
        CircuitOp("CX", (0, 1)),
        CircuitOp("CX", (3, 0)),
        CircuitOp("CX", (0, 2)),
        CircuitOp("CX", (1, 4)),
        CircuitOp("CX", (2, 4)),
        CircuitOp("MEAS_Z", (4,), label="zero_like_flag"),
    ]
    return CliffordCircuit(5, (0, 2, 3), (0, 1, 2, 3), ops, name="steane_zero_like_prop_7_11")


def _colour_dual_op(op: CircuitOp) -> CircuitOp:
    """Apply ZX colour duality (Z <-> X, CNOT direction reversed) to an ancilla/coupling op."""
    swap_kind = {
        "PREP_0": "PREP_PLUS",
        "PREP_PLUS": "PREP_0",
        "PREP_1": "PREP_MINUS",
        "PREP_MINUS": "PREP_1",
        "MEAS_Z": "MEAS_X",
        "MEAS_X": "MEAS_Z",
        "Z": "X",
        "X": "Z",
    }
    if op.kind in swap_kind:
        return CircuitOp(swap_kind[op.kind], op.qubits, op.mode, op.label, dict(op.params))
    if op.kind == "CX":
        c, t = op.qubits
        return CircuitOp("CX", (t, c), op.mode, op.label, dict(op.params))
    if op.kind == "ENC_PREP":
        new_p = dict(op.params)
        # |0_L> (stabilised by logical_z) <-> |+_L> (stabilised by logical_x)
        if "logical_x_dual" in new_p:
            new_p["logical_z"] = new_p.pop("logical_x_dual")
        return CircuitOp("ENC_PREP", op.qubits, op.mode, op.label, new_p)
    return CircuitOp(op.kind, op.qubits, op.mode, op.label, dict(op.params))


def synth_steane_syndrome(
    code: CSSCodeSpec | None = None,
    mode: str = "optimised_5anc",
    sector: str = "X",
) -> tuple[CliffordCircuit, DerivationCertificate]:
    """Synthesize Steane-style syndrome extraction (Props. 7.8, 7.9, 7.12).

    Modes for `steane_7_1_3`:
      - 'standard_8anc' : Goto 8-ancilla |0_L>/|+_L> prep + 7 transversal CNOTs + 7 MEAS_X/MEAS_Z
                          (8 ancilla qubits, 18 CNOTs, 8 measurements; Figure 1a).
      - 'optimised_5anc': Kissinger's bent-leg |0>-like/|+>-like Steane extraction on 5 ancillas
                          (5 ancilla qubits, 15 CNOTs, 5 measurements; Figure 1b / Prop. 7.12).
    For general CSS codes, 'ideal_ancilla' emits Proposition 7.8 using fault-free |0_L> / |+_L> blocks.
    Note: `sector == 'Z'` uses the exact ZX colour-dual circuit (Prop. 7.8 right half, p. 36)
    with zero physical Hadamard gates on data qubits!
    """
    css = code if code is not None else make_steane_7_1_3()
    if sector not in ("X", "Z"):
        raise RefusalError("INVALID-STEANE-SECTOR", f"sector must be 'X' or 'Z', got {sector!r}")
    spec = make_full_syndrome_spec(css, sector=sector)
    init_zx = circuit_to_zx_fe(spec, "circuit_css")

    n = css.n
    data_qs = tuple(range(n))

    if mode == "standard_8anc":
        if css.name != "steane_7_1_3":
            raise RefusalError("UNSUPPORTED-CODE-FOR-GOTO", f"Goto 8-ancilla preset is for steane_7_1_3, got {css.name}")
        goto_c = synth_steane_goto_prep()
        ops: list[CircuitOp] = []
        for op in goto_c.ops:
            ops.append(CircuitOp(op.kind, tuple(n + q for q in op.qubits), op.mode, op.label, dict(op.params)))
        for q in range(7):
            ops.append(CircuitOp("CX", (n + q, q)))
        for q in range(7):
            ops.append(CircuitOp("MEAS_X", (n + q,), label=f"steane_m_{q}"))
        if sector == "Z":
            ops = [_colour_dual_op(op) for op in ops]
        circ = CliffordCircuit(15, data_qs, data_qs, ops, name=f"steane_standard_8anc_{sector}")
        cert = DerivationCertificate(
            name=f"steane_8anc_{sector}_prop_7_8",
            spec_circuit=spec,
            initial_zx=init_zx,
            steps=(RewriteStep("steane_codespace_d", {"d": 3}, note="Prop C.1 + Goto [26]"),),
            extracted_circuit=circ,
            target_noise_model="circuit_css",
            declared_w=3,
        )
        return circ, cert

    if mode == "optimised_5anc":
        if css.name != "steane_7_1_3":
            raise RefusalError("UNSUPPORTED-CODE-FOR-5ANC", f"5-ancilla preset is for steane_7_1_3, got {css.name}")
        # Ancilla qubits: b0=7, b1=8, b2=9, b3=10, b4=11 (flag)
        # Bent qubits: d0=0 (coupled to b0=7), d5=5 (coupled to b2=9), d6=6 (coupled to b3=10)
        # Un-bent qubits: d1=1 (coupled from b0=7), d2=2 (from b1=8), d3=3 (from b2=9), d4=4 (from b3=10)
        b0, b1, b2, b3, b4 = 7, 8, 9, 10, 11
        ops = [
            CircuitOp("PREP_PLUS", (b0,)),
            CircuitOp("PREP_0", (b1,)),
            CircuitOp("PREP_PLUS", (b2,)),
            CircuitOp("PREP_PLUS", (b3,)),
            CircuitOp("PREP_0", (b4,)),
            # 3 bent-leg input CNOTs (Prop 6.10 split tiles!):
            CircuitOp("CX", (b0, 0)),
            CircuitOp("CX", (b2, 5)),
            CircuitOp("CX", (b3, 6)),
            # 8 internal CNOTs of Prop 7.11 |0>-like linear map:
            CircuitOp("CX", (b2, b1)),
            CircuitOp("CX", (b2, b3)),
            CircuitOp("CX", (b3, b4)),
            CircuitOp("CX", (b0, b1)),
            CircuitOp("CX", (b3, b0)),
            CircuitOp("CX", (b0, b2)),
            CircuitOp("CX", (b1, b4)),
            CircuitOp("CX", (b2, b4)),
            CircuitOp("MEAS_Z", (b4,), label="steane5_flag"),
            # 4 un-bent output CNOTs + MEAS_X (colour-dual Prop 6.9 merge tiles!):
            CircuitOp("CX", (b0, 1)),
            CircuitOp("CX", (b1, 2)),
            CircuitOp("CX", (b2, 3)),
            CircuitOp("CX", (b3, 4)),
            CircuitOp("MEAS_X", (b0,), label="steane5_w0"),
            CircuitOp("MEAS_X", (b1,), label="steane5_w1"),
            CircuitOp("MEAS_X", (b2,), label="steane5_w2"),
            CircuitOp("MEAS_X", (b3,), label="steane5_w3"),
        ]
        if sector == "Z":
            ops = [_colour_dual_op(op) for op in ops]
        circ = CliffordCircuit(12, data_qs, data_qs, ops, name=f"steane_optimised_5anc_{sector}_prop_7_12")
        cert = DerivationCertificate(
            name=f"steane_5anc_{sector}_prop_7_12",
            spec_circuit=spec,
            initial_zx=init_zx,
            steps=(RewriteStep("steane_codespace_d", {"d": 3}, note="Prop C.1 + Prop 7.9 + Prop 7.11"),),
            extracted_circuit=circ,
            target_noise_model="circuit_css",
            declared_w=3,
        )
        return circ, cert

    if mode == "ideal_ancilla":
        # Proposition 7.8 for arbitrary CSS code [[n, k, d]] with fault-free |0_L^k> / |+_L^k> ancilla!
        anc_qs = tuple(range(n, 2 * n))
        ops = [
            CircuitOp(
                "ENC_PREP",
                anc_qs,
                mode="ff",
                params={
                    "input_qubits": [],
                    "stabilizers": css.stabilizers_pauli_strings(),
                    "logical_z": css.logical_z_strings() if sector == "X" else css.logical_x_strings(),
                },
            )
        ]
        for q in range(n):
            ops.append(CircuitOp("CX", (n + q, q) if sector == "X" else (q, n + q)))
        for q in range(n):
            ops.append(CircuitOp("MEAS_X" if sector == "X" else "MEAS_Z", (n + q,), label=f"steane_m_{q}"))
        circ = CliffordCircuit(2 * n, data_qs, data_qs, ops, name=f"steane_ideal_anc_{css.name}_{sector}")
        cert = DerivationCertificate(
            name=f"steane_ideal_anc_{css.name}_{sector}_prop_7_8",
            spec_circuit=spec,
            initial_zx=init_zx,
            steps=(RewriteStep("steane_codespace_d", {"d": css.d}, note="Prop C.1"),),
            extracted_circuit=circ,
            target_noise_model="circuit_css",
            declared_w=css.d,
        )
        return circ, cert

    raise RefusalError("UNKNOWN-STEANE-MODE", f"unknown Steane synthesis mode {mode!r}")


def synth_sequential_cat4_steane_bad() -> CliffordCircuit:
    """Unflagged cat_4 Steane X-syndrome check (planted hook-error negative control!).

    Couples unflagged 4-qubit H|cat_4> ancillas to the three weight-4 Steane X-stabilisers:
    a single Z fault inside an unflagged H|cat_4> prep propagates to a weight-2 Z hook error
    on the Steane data block without triggering any detector, failing w=2 fault equivalence
    even when restricted to the Steane codespace!
    """
    code = make_steane_7_1_3()
    ops: list[CircuitOp] = [CircuitOp("H", (q,)) for q in range(7)]
    next_q = 7
    for row in code.hx_rows:
        act_qs = [q for q in range(7) if (row >> q) & 1]
        bad_hcat = synth_unflagged_cat_state(4, basis="X")
        remap = {q: next_q + q for q in range(4)}
        for op in bad_hcat.ops:
            ops.append(CircuitOp(op.kind, tuple(remap[q] for q in op.qubits), op.mode, op.label, dict(op.params)))
        for idx_a, dq in enumerate(act_qs):
            ops.append(CircuitOp("CX", (dq, next_q + idx_a)))
        for idx_a in range(4):
            ops.append(CircuitOp("MEAS_Z", (next_q + idx_a,), label=f"bad_cat4_m_{next_q + idx_a}"))
        next_q += 4
    ops.extend(CircuitOp("H", (q,)) for q in range(7))
    return CliffordCircuit(next_q, tuple(range(7)), tuple(range(7)), ops, name="bad_unflagged_cat4_steane")
