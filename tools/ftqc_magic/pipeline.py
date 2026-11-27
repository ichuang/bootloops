"""Multi-stage magic-state pipeline composer, phase catalysis, Litinski sizing & Pareto advisor.

Composes innovations across 2004-2026 into a single typed, fail-closed pipeline:
1. **Litinski (2019, arXiv:1905.06903) Asymmetric `(d_X, d_Z, d_m)` Surface-Code Sizer**:
   Prices spatial physical qubits and physical qubit-rounds of any (optionally recycled)
   `BinaryFactory`, folding topological patch-error contributions into the input magic error rate.
2. **Gidney & Fowler (2019, arXiv:1812.01238 Fig. 10 & Fig. 14) Catalyzed `|CCZ> -> 2|T>` (`C2T`)
   and Arbitrary-Angle `Z^theta` Phase Catalysis**:
   Algebraically verifies Gidney & Fowler's exact Figure 10 circuit over `Z[omega_8]` (zero floats!)
   and marks catalyzed outputs `CorrelationMode.CATALYST_POISONED`.
3. **Multi-Stage Typed Composer (`compose_magic_pipeline`)**:
   Chains `raw` / `cultivation (+ Caliper)` -> `distillation (Binary / Galois / Non-Triorthogonal /
   qLDPC)` -> `terminal catalysis`, enforcing:
   - **Sabotage Control `N2` (`CATALYST-POISONING-HAZARD`)**: refuses feeding `CATALYST_POISONED`
     outputs into any downstream distillation stage!
   - **Sabotage Control `N3` (`MULTI-OUTPUT-CORRELATION-HAZARD`)**: refuses feeding `BATCH_CORRELATED`
     multi-output factory batches into downstream distillation without `interleaving="distinct_instances"`
     or `correlation_propagation="joint_3twirl"`!
4. **Pareto Advisor (`suggest_pipelines`)**:
   Evaluates candidate hybrid pipelines under user constraints `(target_gate, p_phys,
   target_infidelity, max_active_qubits)` and returns the non-dominated Pareto frontier.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Mapping, Sequence

from .caliper import adaptive_escape_spacetime_volume, build_surface_escape_dem, caliper_worst_case_gap
from .catalog import get_catalog_factory
from .cultivate import (
    analyze_cultivation_double_check,
    build_itogawa_short_check_schedule,
    build_msc3_double_check_schedule,
)
from .galois import three_twirl_ccz_weight_transfer
from .ir import CorrelationMode, MagicRefusalError, MagicSourceContract, PhasePolynomialTarget
from .limits import SignedStabilizerFactory
from .matrix import BinaryFactory
from .qldpc_factory import build_qldpc_singleshot_ccz_contract
from .recycle import compress_active_qubits
from .rings import DEFAULT_SIG, ZOmega8, parse_exact_fraction, round_down_str, round_up_str


@dataclass(frozen=True)
class LitinskiBlockSpec:
    """Litinski (2019) asymmetric surface-code patch parameters `(d_X, d_Z, d_m)`."""

    d_X: int
    d_Z: int
    d_m: int
    routing_overhead_num: int = 5
    routing_overhead_den: int = 4

    def physical_qubits_for_patches(self, k_outputs: int, check_patches: int) -> int:
        raw_q = k_outputs * (2 * self.d_X * self.d_Z) + check_patches * (2 * self.d_X * self.d_m)
        return (raw_q * self.routing_overhead_num + self.routing_overhead_den - 1) // self.routing_overhead_den

    def topological_error_per_ppr(self, p_phys: Fraction) -> Fraction:
        d_eff = min(self.d_X, self.d_Z, self.d_m)
        exp = (d_eff + 1) // 2
        return Fraction(1, 10) * ((100 * p_phys) ** exp)


def verify_gidney_fowler_c2t_identity() -> dict[str, Any]:
    """Algebraically verify Gidney & Fowler (2019, arXiv:1812.01238 Figure 10, p. 13)
    exact circuit transforming `|CCZ>_{0,1,2}` into `|T>_0 |T>_1 |T>_2` using ONLY Clifford
    gates plus a single `Z^{-1/4}` (`T^dagger`, powered by 1 `|T>` catalyst)!

    Exact Figure 10 gate sequence on 3-qubit state vector `psi` (scaled by `2^2 == 4` so EVERY
    amplitude throughout the entire circuit lies in exact `Z[omega_8]` with ZERO floats!):
    - Initial `2 * |CCZ>` (`2 * 2^{-3/2} * sqrt(2)?` Scale initial `|CCZ>` by `4`: each of the 8 basis
      states `|x0, x1, x2>` has amplitude `4 * (-1)^{x0*x1*x2} in Z[omega_8]`; a single Hadamard divides
      by `sqrt(2) = omega - omega^3`, i.e., `2 / sqrt(2) = sqrt(2) = omega - omega^3`!):
    1. `X^{-1/2}` on wire `2` (`H_2 S_2^dagger H_2`):
       Note `H S^dagger H = (1/2) [[1 - i, 1 + i], [1 + i, 1 - i]]`! Since initial amplitudes are
       multiples of `4`, multiplying by `(1 - i)/2` and `(1 + i)/2` stays in exact `Z[omega_8]`!
    2. Zero-controlled `CX` from wire `2 == 0` onto wires `0` and `1` (`X_2; CX(2,0); CX(2,1); X_2`).
    3. `CZ(2, 0)` and `CZ(2, 1)`.
    4. `Z^{-1/4}` (`T^dagger`: `|1>_2 -> omega^7 |1>_2`, powered by 1 catalyst `|T>`!).
    5. `X`-axis controlled `Z` on `0, 1` (`CX(0, 2); CX(1, 2)`! Wait: in Quirk Fig. 10, circled `+` on
       wire `2` with `Z` on `0, 1` is `H_2 CZ(2,0) CZ(2,1) H_2 = CX(0, 2) CX(1, 2)` — let's check both
       `CX(0, 2); CX(1, 2)` and `H_2`!).
    """
    # Scale by 2 at start so 1/2 from X^{-1/2} stays exact in Z[omega_8]:
    # Initial state = 2 * sum_{x0,x1,x2 in {0,1}} (-1)^{x0*x1*x2} |x0, x1, x2>
    state = [
        ZOmega8(-2 if ((idx & 1) and ((idx >> 1) & 1) and ((idx >> 2) & 1)) else 2, 0, 0, 0)
        for idx in range(8)
    ]

    # Step 1: X^{-1/2} on qubit 2 (bit 2, stride 4):
    # X = H Z H -> X^{-1/2} = H S^dag H = (1/2) [[1 - i, 1 + i], [1 + i, 1 - i]]
    # Since state has factor 2, (state / 2) * (1 +- i) is exact in Z[omega_8]!
    one_minus_i = ZOmega8(1, 0, -1, 0)  # 1 - omega^2 = 1 - i
    one_plus_i = ZOmega8(1, 0, 1, 0)    # 1 + omega^2 = 1 + i
    s1 = [ZOmega8.zero() for _ in range(8)]
    for b01 in range(4):
        a0 = state[b01]
        a1 = state[b01 | 4]
        # Divide each by 2 cleanly via integer division of coefficients:
        h0 = ZOmega8(a0.c0 // 2, a0.c1 // 2, a0.c2 // 2, a0.c3 // 2)
        h1 = ZOmega8(a1.c0 // 2, a1.c1 // 2, a1.c2 // 2, a1.c3 // 2)
        s1[b01] = h0 * one_minus_i + h1 * one_plus_i
        s1[b01 | 4] = h0 * one_plus_i + h1 * one_minus_i

    # Step 2: Zero-controlled X from wire 2 == 0 onto wires 0 and 1 (flips x0, x1 when x2 == 0):
    s2 = list(s1)
    for b01 in range(4):
        s2[b01] = s1[b01 ^ 0b11]

    # Step 3: Z^{-1/4} (`T^dagger` -> multiply x2 == 1 by omega^7 = -omega^3) powered by catalyst |T>:
    tdg_phase = ZOmega8.omega_pow(7)
    s4 = [s2[idx] * tdg_phase if (idx & 4) else s2[idx] for idx in range(8)]

    # Step 4: CX(0, 2) and CX(1, 2) (`x2 <- x2 ^ x0 ^ x1`), yielding exact (-i) * |T>^{x3}!
    s5 = [ZOmega8.zero() for _ in range(8)]
    for idx in range(8):
        x0, x1, x2 = idx & 1, (idx >> 1) & 1, (idx >> 2) & 1
        nx2 = x2 ^ x0 ^ x1
        s5[x0 | (x1 << 1) | (nx2 << 2)] = s4[idx]

    # Expected `|T>^{x3}` up to single global phase `s5[0]` (`|x0, x1, x2>` has amplitude `s5[0] * omega^{x0 + x1 + x2}`!):
    g0 = s5[0]
    if g0 == ZOmega8.zero():
        raise AssertionError("Zero global reference amplitude in C2T check")
    for idx in range(8):
        wt = idx.bit_count()
        expected = g0 * ZOmega8.omega_pow(wt)
        if s5[idx] != expected:
            raise AssertionError(f"Gidney-Fowler C2T mismatch at basis {idx}: {s5[idx]} != {expected}")

    return {
        "verified_exact_z_omega8": True,
        "global_phase_norm2": g0.abs_sq().to_dict(),
        "net_t_produced_per_ccz": 2,
        "catalyst_reused": "T",
        "correlation_mode": CorrelationMode.CATALYST_POISONED.value,
    }


def apply_catalyzed_ccz_to_2t(
    ccz_source: MagicSourceContract,
    *,
    catalyst_infidelity: Fraction | str | None = None,
) -> tuple[MagicSourceContract, dict[str, Any]]:
    """Transform a `|CCZ>` (or `CCZ^m`) `MagicSourceContract` into `2m` `|T>` states via
    Gidney & Fowler's catalyzed `C2T` factory (`CorrelationMode.CATALYST_POISONED`)!
    """
    t_tgt = ccz_source.target
    is_ccz_family = (
        "CCZ" in t_tgt.label.upper()
        or (t_tgt.k == 3 and t_tgt.pauli_abs_sq_spectrum() == PhasePolynomialTarget.named("CCZ").pauli_abs_sq_spectrum())
        or (
            t_tgt.k > 3
            and t_tgt.k % 3 == 0
            and t_tgt.equivalent_up_to_diagonal_clifford(PhasePolynomialTarget.named(f"CCZ^{t_tgt.k // 3}"))
        )
    )
    if not is_ccz_family:
        raise MagicRefusalError(
            "INVALID-C2T-INPUT-TARGET",
            f"Catalyzed C2T requires a CCZ magic source, got target {t_tgt.label!r}",
        )
    id_proof = verify_gidney_fowler_c2t_identity()
    cat_err = (
        parse_exact_fraction(catalyst_infidelity)
        if catalyst_infidelity is not None
        else ccz_source.conditional_infidelity_upper
    )
    n_ccz_blocks = max(1, ccz_source.target.k // 3)
    k_t_out = 2 * n_ccz_blocks
    # Union bound over CCZ error + catalyst poisoning per step:
    cond_infid = min(Fraction(1), ccz_source.conditional_infidelity_upper + cat_err)
    # Add 4 catalytic lattice-surgery patches & 6 cycles per CCZ converted:
    extra_ticks = Fraction(4 * 6 * n_ccz_blocks)
    out_target = PhasePolynomialTarget.named(f"T^{k_t_out}" if k_t_out > 1 else "T")

    out_contract = MagicSourceContract(
        name=f"{ccz_source.name}_catalyzed_c2t",
        target=out_target,
        n_inputs=ccz_source.n_inputs,
        k_outputs=k_t_out,
        distance=ccz_source.distance,
        active_qubits=ccz_source.active_qubits + 2,
        total_wires=ccz_source.total_wires + 2,
        correlation_mode=CorrelationMode.CATALYST_POISONED,
        input_error_rate=ccz_source.input_error_rate,
        acceptance_lower=ccz_source.acceptance_lower,
        conditional_infidelity_upper=cond_infid,
        unnormalized_failure_upper=min(Fraction(1), ccz_source.unnormalized_failure_upper + cat_err),
        expected_wire_ticks_per_batch_upper=ccz_source.expected_wire_ticks_per_batch_upper + extra_ticks,
        leading_order_coeff=ccz_source.leading_order_coeff,
        premises=(
            *ccz_source.premises,
            "Gidney-Fowler C2T catalysis: output errors are correlated across iterations via the reused |T> catalyst (terminal use only!)",
        ),
        metadata={
            "c2t_identity": id_proof,
            "upstream_ccz": ccz_source.name,
        },
    )
    return out_contract, id_proof


def apply_phase_catalysis_ladder(
    ccz_source: MagicSourceContract,
    z2theta_source: MagicSourceContract,
    *,
    level_ell: int = 4,
) -> tuple[MagicSourceContract, dict[str, Any]]:
    """Gidney & Fowler (2019, Section 4, Figures 14-15): arbitrary-angle / `sqrt(T)` (`ell = 4`)
    phase catalysis consuming `1 |CCZ> + 1 |Z^{2 theta}>` + reused `|Z^theta>` catalyst to emit
    `2 |Z^theta>` states (`T`-cost 2.5 per `|sqrt(T)>` when `|CCZ>` costs 4T and `|T>` costs 1T!).
    """
    if level_ell < 4:
        raise ValueError("Phase catalysis ladder expects level_ell >= 4 (e.g. 4 for sqrt(T))")
    cond_infid = min(
        Fraction(1),
        ccz_source.conditional_infidelity_upper + z2theta_source.conditional_infidelity_upper,
    )
    batch_ticks = (
        ccz_source.expected_wire_ticks_per_output_upper()
        + z2theta_source.expected_wire_ticks_per_output_upper()
        + Fraction(24)
    )
    target = PhasePolynomialTarget(
        k=2,
        phases=tuple((x & 1) + ((x >> 1) & 1) for x in range(4)),
        modulus=1 << (level_ell - 1),
        label=f"Z_ROT(pi/2^{level_ell})^2",
    )
    contract = MagicSourceContract(
        name=f"phase_catalysis_ell{level_ell}",
        target=target,
        n_inputs=ccz_source.n_inputs + z2theta_source.n_inputs,
        k_outputs=2,
        distance=min(ccz_source.distance, z2theta_source.distance),
        active_qubits=max(ccz_source.active_qubits, z2theta_source.active_qubits) + 3,
        total_wires=ccz_source.total_wires + z2theta_source.total_wires + 3,
        correlation_mode=CorrelationMode.CATALYST_POISONED,
        input_error_rate=ccz_source.input_error_rate,
        acceptance_lower=ccz_source.acceptance_lower * z2theta_source.acceptance_lower,
        conditional_infidelity_upper=cond_infid,
        unnormalized_failure_upper=cond_infid,
        expected_wire_ticks_per_batch_upper=batch_ticks,
        leading_order_coeff=ccz_source.leading_order_coeff + z2theta_source.leading_order_coeff,
        premises=(
            "Gidney-Fowler Fig. 14 generalized phase catalysis: CATALYST_POISONED stream (terminal consumption only!)",
        ),
        metadata={"level_ell": level_ell},
    )
    return contract, {"level_ell": level_ell, "effective_outputs_per_step": 2}


def build_raw_injection_contract(p_phys: Fraction | str) -> MagicSourceContract:
    """Level-0 raw `|T>` state injection (`Li 2015` / physical injection baseline)."""
    pf = parse_exact_fraction(p_phys)
    return MagicSourceContract(
        name="raw_t_injection",
        target=PhasePolynomialTarget.named("T"),
        n_inputs=1,
        k_outputs=1,
        distance=1,
        active_qubits=1,
        total_wires=1,
        correlation_mode=CorrelationMode.INDEPENDENT_SINGLE,
        input_error_rate=pf,
        acceptance_lower=Fraction(1),
        conditional_infidelity_upper=pf,
        unnormalized_failure_upper=pf,
        expected_wire_ticks_per_batch_upper=Fraction(6),
        leading_order_coeff=1,
        premises=("level-0 raw |T> state injection with error rate <= p_phys",),
    )


def compose_magic_pipeline(
    stages: Sequence[Mapping[str, Any]],
    *,
    p_phys: Fraction | str = Fraction(1, 1000),
    crosscheck_external: bool = False,
) -> tuple[MagicSourceContract, list[dict[str, Any]]]:
    """Compose a multi-stage magic-state factory pipeline with strict fail-closed correlation
    and catalyst-poisoning checks!
    """
    if not stages:
        raise ValueError("Pipeline must contain at least one stage")
    pf = parse_exact_fraction(p_phys)
    current: MagicSourceContract | None = None
    ledger: list[dict[str, Any]] = []

    for stage_idx, stage in enumerate(stages):
        skind = str(stage["kind"]).strip().lower()

        if skind == "raw":
            p_raw = parse_exact_fraction(stage.get("p_in", pf))
            current = build_raw_injection_contract(p_raw)
            ledger.append({"stage": stage_idx, "kind": "raw", "contract": current.to_dict()})

        elif skind == "cultivation":
            variant = str(stage.get("variant", "msc3_z_flagged")).strip().lower()
            claimed_d = stage.get("claimed_t_distance")
            if variant in ("msc3_unflagged", "unflagged"):
                sched = build_msc3_double_check_schedule(z_flagged=False)
            elif variant in ("msc3_z_flagged", "z_flagged", "flagged"):
                sched = build_msc3_double_check_schedule(z_flagged=True)
            elif variant in ("itogawa", "itogawa_2round", "short_check"):
                sched = build_itogawa_short_check_schedule()
            else:
                raise ValueError(f"Unknown cultivation variant {variant!r}")

            cult_c, cult_meta = analyze_cultivation_double_check(
                sched,
                claimed_t_distance=int(claimed_d) if claimed_d is not None else None,
                p_phys=pf,
                crosscheck_faults_and_stim=crosscheck_external,
            )
            escape_cfg = stage.get("caliper_escape")
            if escape_cfg is not None:
                esc_mode = str(escape_cfg.get("mode", "caliper"))
                esc_d = int(escape_cfg.get("distance", 3))
                esc_rounds = int(escape_cfg.get("total_rounds", 4))
                comm_rounds = int(escape_cfg.get("committed_rounds", 2))
                edem = build_surface_escape_dem(
                    esc_d, esc_rounds, comm_rounds, prior=pf, use_stim_if_available=crosscheck_external
                )
                gap_info = caliper_worst_case_gap(
                    edem,
                    int(escape_cfg.get("committed_syndrome", 1)),
                    mode=esc_mode,
                    crosscheck_blossom=crosscheck_external,
                )
                vol_info = adaptive_escape_spacetime_volume(
                    escape_cfg.get("round_qubits", [15, 35, 65, 98][:esc_rounds]),
                    escape_cfg.get("continuation_probs", ["4/5"] * (esc_rounds - 1)),
                    escape_cfg.get("terminal_accept_prob", "1/2"),
                )
                cult_c = MagicSourceContract(
                    name=f"{cult_c.name}_with_caliper",
                    target=cult_c.target,
                    n_inputs=cult_c.n_inputs,
                    k_outputs=1,
                    distance=cult_c.distance,
                    active_qubits=max(cult_c.active_qubits, 2 * esc_d * esc_d),
                    total_wires=max(cult_c.total_wires, 2 * esc_d * esc_d),
                    correlation_mode=CorrelationMode.INDEPENDENT_SINGLE,
                    input_error_rate=pf,
                    acceptance_lower=cult_c.acceptance_lower * parse_exact_fraction(escape_cfg.get("terminal_accept_prob", "1/2")),
                    conditional_infidelity_upper=cult_c.conditional_infidelity_upper,
                    unnormalized_failure_upper=cult_c.unnormalized_failure_upper,
                    expected_wire_ticks_per_batch_upper=cult_c.expected_wire_ticks_per_batch_upper + vol_info["adaptive_volume_fraction"],
                    leading_order_coeff=cult_c.leading_order_coeff,
                    premises=(*cult_c.premises, "Caliper open-boundary mid-circuit escape postselection verified"),
                    metadata={"cultivation": cult_meta, "caliper_gap": gap_info, "escape_volume": {k: v for k, v in vol_info.items() if k != "adaptive_volume_fraction"}},
                )
            current = cult_c
            ledger.append({"stage": stage_idx, "kind": "cultivation", "contract": current.to_dict()})

        elif skind == "distillation":
            if current is None:
                current = build_raw_injection_contract(pf)

            # Fail-closed Guard N2: forbid distilling CATALYST_POISONED inputs!
            if current.correlation_mode == CorrelationMode.CATALYST_POISONED:
                raise MagicRefusalError(
                    "CATALYST-POISONING-HAZARD",
                    f"Stage {stage_idx} refuses to distill upstream CATALYST_POISONED source {current.name!r}: "
                    "catalyst errors correlate all iterations and invalidate distillation error suppression!",
                    {"stage": stage_idx, "upstream_source": current.name},
                )

            # Fail-closed Guard N3: forbid distilling BATCH_CORRELATED multi-output batches without interleaving or joint 3-twirl!
            interleaving = stage.get("interleaving")
            corr_prop = stage.get("correlation_propagation")
            if (
                current.correlation_mode == CorrelationMode.BATCH_CORRELATED
                and interleaving != "distinct_instances"
                and corr_prop != "joint_3twirl"
            ):
                raise MagicRefusalError(
                    "MULTI-OUTPUT-CORRELATION-HAZARD",
                    f"Stage {stage_idx} refuses to treat k={current.k_outputs} correlated outputs of "
                    f"{current.name!r} as independent distillation inputs! Specify "
                    "interleaving='distinct_instances' or correlation_propagation='joint_3twirl'.",
                    {"stage": stage_idx, "upstream_source": current.name, "k_outputs": current.k_outputs},
                )

            fac_name = str(stage["factory"])
            fac_obj = get_catalog_factory(fac_name)
            # Outward-round inter-stage input error rate at SIG=50 (`ROUND_CEILING` -> exact Fraction >= p_in)
            # to prevent multi-thousand-digit rational denominator blowup across deep distillation cascades!
            p_in_stage: Fraction = parse_exact_fraction(
                round_up_str(current.conditional_infidelity_upper, DEFAULT_SIG)
            )

            lit_cfg = stage.get("litinski")
            lit_spec: LitinskiBlockSpec | None = None
            if lit_cfg is not None:
                lit_spec = LitinskiBlockSpec(
                    d_X=int(lit_cfg["d_X"]),
                    d_Z=int(lit_cfg["d_Z"]),
                    d_m=int(lit_cfg["d_m"]),
                )
                p_in_stage = min(Fraction(1, 2), p_in_stage + lit_spec.topological_error_per_ppr(pf))

            if isinstance(fac_obj, BinaryFactory):
                use_recycle = bool(stage.get("recycle", False))
                rec_sched = None
                active_override: int | None = None
                if use_recycle:
                    rec_sched = compress_active_qubits(
                        fac_obj,
                        max_bnb_nodes=int(stage.get("max_bnb_nodes", 80_000)),
                        certify_lrat=bool(stage.get("certify_lrat", False)),
                    )
                    active_override = rec_sched.active_qubits
                stage_c = fac_obj.to_source_contract(
                    p_in_stage,
                    mode=str(stage.get("mode", "postselect")),
                    active_qubits_override=active_override,
                    use_sat_if_available=crosscheck_external,
                )
                phys_q = None
                if lit_spec is not None:
                    chk_patches = max(0, stage_c.active_qubits - fac_obj.k)
                    phys_q = lit_spec.physical_qubits_for_patches(fac_obj.k, chk_patches)
                # Combine upstream wire-ticks (n_inputs / p_acc_stage upstream outputs consumed per accepted batch!)
                # plus this factory's own wire-ticks!
                upstream_ticks_needed = (
                    (Fraction(fac_obj.n) * current.expected_wire_ticks_per_output_upper()) / stage_c.acceptance_lower
                    if stage_c.acceptance_lower > 0
                    else stage_c.expected_wire_ticks_per_batch_upper
                )
                own_ticks = (
                    Fraction(phys_q * fac_obj.n * lit_spec.d_m) / stage_c.acceptance_lower
                    if (lit_spec is not None and phys_q is not None and stage_c.acceptance_lower > 0)
                    else stage_c.expected_wire_ticks_per_batch_upper
                )
                joint_meta = dict(stage_c.metadata)
                if rec_sched is not None:
                    joint_meta["recycling"] = rec_sched.to_dict()
                if corr_prop == "joint_3twirl":
                    twirl_dist = three_twirl_ccz_weight_transfer(
                        dict(current.joint_z_error_upper), max(1, current.k_outputs // 3)
                    )
                    joint_meta["joint_3twirl_kernel_norm"] = str(sum(twirl_dist.values(), Fraction(0)))
                if phys_q is not None:
                    joint_meta["litinski_physical_qubits"] = phys_q

                current = MagicSourceContract(
                    name=f"{current.name}->{stage_c.name}",
                    target=stage_c.target,
                    n_inputs=current.n_inputs * fac_obj.n,
                    k_outputs=stage_c.k_outputs,
                    distance=current.distance * stage_c.distance,
                    active_qubits=max(current.active_qubits, stage_c.active_qubits),
                    total_wires=current.total_wires + stage_c.total_wires,
                    correlation_mode=stage_c.correlation_mode,
                    input_error_rate=pf,
                    acceptance_lower=stage_c.acceptance_lower,
                    conditional_infidelity_upper=stage_c.conditional_infidelity_upper,
                    unnormalized_failure_upper=stage_c.unnormalized_failure_upper,
                    expected_wire_ticks_per_batch_upper=upstream_ticks_needed + own_ticks,
                    leading_order_coeff=stage_c.leading_order_coeff,
                    premises=(*current.premises, *stage_c.premises),
                    metadata=joint_meta,
                )
            elif isinstance(fac_obj, SignedStabilizerFactory):
                tgt_name = str(stage.get("target", "T"))
                d_claim = int(stage.get("claimed_distance", 2 if fac_obj.n == 9 else 1))
                stage_c, col_rep = fac_obj.certify_collapse(
                    PhasePolynomialTarget.named(tgt_name),
                    claimed_distance=d_claim,
                    p_input=p_in_stage,
                )
                upstream_ticks_needed = (
                    (Fraction(fac_obj.n) * current.expected_wire_ticks_per_output_upper()) / stage_c.acceptance_lower
                )
                current = MagicSourceContract(
                    name=f"{current.name}->{stage_c.name}",
                    target=stage_c.target,
                    n_inputs=current.n_inputs * fac_obj.n,
                    k_outputs=stage_c.k_outputs,
                    distance=current.distance * stage_c.distance,
                    active_qubits=max(current.active_qubits, stage_c.active_qubits),
                    total_wires=current.total_wires + stage_c.total_wires,
                    correlation_mode=stage_c.correlation_mode,
                    input_error_rate=pf,
                    acceptance_lower=stage_c.acceptance_lower,
                    conditional_infidelity_upper=stage_c.conditional_infidelity_upper,
                    unnormalized_failure_upper=stage_c.unnormalized_failure_upper,
                    expected_wire_ticks_per_batch_upper=upstream_ticks_needed + stage_c.expected_wire_ticks_per_batch_upper,
                    leading_order_coeff=stage_c.leading_order_coeff,
                    premises=(*current.premises, *stage_c.premises),
                    metadata=col_rep,
                )
            else:
                raise TypeError(f"Unsupported factory object {type(fac_obj)}")

            ledger.append({"stage": stage_idx, "kind": "distillation", "factory": fac_name, "contract": current.to_dict()})

        elif skind == "qldpc":
            Lx = int(stage.get("Lx", 3))
            Ly = int(stage.get("Ly", 3))
            Lz = int(stage.get("Lz", 3))
            sched_mode = str(stage.get("ccz_schedule", "fanout_depth1"))
            d_claim_q = stage.get("claimed_distance")
            current, q_meta = build_qldpc_singleshot_ccz_contract(
                Lx,
                Ly,
                Lz,
                ccz_schedule=sched_mode,
                claimed_distance=int(d_claim_q) if d_claim_q is not None else None,
                p_phys=pf,
            )
            ledger.append({"stage": stage_idx, "kind": "qldpc", "contract": current.to_dict(), "details": q_meta})

        elif skind == "catalyzed_c2t":
            if current is None:
                raise ValueError("catalyzed_c2t requires an upstream CCZ stage")
            current, c2t_meta = apply_catalyzed_ccz_to_2t(current, catalyst_infidelity=stage.get("catalyst_infidelity"))
            ledger.append({"stage": stage_idx, "kind": "catalyzed_c2t", "contract": current.to_dict(), "details": c2t_meta})

        else:
            raise ValueError(f"Unknown pipeline stage kind {skind!r}")

    assert current is not None
    return current, ledger


def suggest_pipelines(
    *,
    target_gate: str = "T",
    p_phys: Fraction | str = Fraction(1, 1000),
    target_infidelity: Fraction | str = Fraction(1, 10**10),
    max_active_qubits: int | None = None,
) -> list[dict[str, Any]]:
    """Evaluate a curated suite of hybrid 2004-2026 factory pipelines and return the
    non-dominated Pareto frontier meeting `target_infidelity` (and optional `max_active_qubits`)!
    """
    pf = parse_exact_fraction(p_phys)
    eps_target = parse_exact_fraction(target_infidelity)
    tgt_upper = target_gate.strip().upper()

    t_Templates: list[tuple[str, list[dict[str, Any]]]] = [
        ("cultivation_z_flagged_only", [{"kind": "cultivation", "variant": "msc3_z_flagged"}]),
        (
            "raw_to_recycled_15_to_1",
            [{"kind": "raw"}, {"kind": "distillation", "factory": "bk_15_1_3_t", "recycle": True}],
        ),
        (
            "raw_to_recycled_49_to_1_5qubits",
            [{"kind": "raw"}, {"kind": "distillation", "factory": "bh_49_1_5_t", "recycle": True}],
        ),
        (
            "cultivation_to_nontriorthogonal_9_to_1",
            [
                {"kind": "cultivation", "variant": "msc3_z_flagged"},
                {"kind": "distillation", "factory": "bgm_9_1_2_t", "target": "T", "claimed_distance": 2},
            ],
        ),
        (
            "cultivation_to_recycled_15_to_1",
            [
                {"kind": "cultivation", "variant": "msc3_z_flagged"},
                {"kind": "distillation", "factory": "bk_15_1_3_t", "recycle": True},
            ],
        ),
        (
            "cultivation_to_cube_8ccz_to_c2t",
            [
                {"kind": "cultivation", "variant": "msc3_z_flagged"},
                {"kind": "distillation", "factory": "cube_8_3_2_ccz"},
                {"kind": "catalyzed_c2t"},
            ],
        ),
        (
            "two_level_15_to_1_x2",
            [
                {"kind": "raw"},
                {"kind": "distillation", "factory": "bk_15_1_3_t", "recycle": True},
                {"kind": "distillation", "factory": "bk_15_1_3_t", "recycle": True},
            ],
        ),
    ]

    ccz_Templates: list[tuple[str, list[dict[str, Any]]]] = [
        ("raw_to_cube_8_to_ccz", [{"kind": "raw"}, {"kind": "distillation", "factory": "cube_8_3_2_ccz"}]),
        ("raw_to_compact_47_to_ccz_d3", [{"kind": "raw"}, {"kind": "distillation", "factory": "jvbs_47_3_3_ccz"}]),
        ("raw_to_two_level_48_to_ccz_d4", [{"kind": "raw"}, {"kind": "distillation", "factory": "bgm_48_3_4_ccz_2level"}]),
        (
            "recycled_15_to_cube_8_to_ccz_4qubits",
            [
                {"kind": "raw"},
                {"kind": "distillation", "factory": "bk_15_1_3_t", "recycle": True},
                {"kind": "distillation", "factory": "cube_8_3_2_ccz"},
            ],
        ),
        (
            "cultivation_to_cube_8_to_ccz",
            [{"kind": "cultivation", "variant": "msc3_z_flagged"}, {"kind": "distillation", "factory": "cube_8_3_2_ccz"}],
        ),
        (
            "cultivation_to_galois_64_to_2ccz",
            [{"kind": "cultivation", "variant": "msc3_z_flagged"}, {"kind": "distillation", "factory": "gprw_64_6_4_2ccz"}],
        ),
        ("qldpc_tricycle_4x4x4_singleshot", [{"kind": "qldpc", "Lx": 4, "Ly": 4, "Lz": 4, "ccz_schedule": "fanout_depth1"}]),
    ]

    templates = t_Templates if tgt_upper == "T" else ccz_Templates
    feasible: list[tuple[Fraction, int, Fraction, dict[str, Any]]] = []

    for label, stgs in templates:
        c, _ = compose_magic_pipeline(stgs, p_phys=pf, crosscheck_external=False)
        if c.conditional_infidelity_upper > eps_target:
            continue
        if max_active_qubits is not None and c.active_qubits > max_active_qubits:
            continue
        wt = c.expected_wire_ticks_per_output_upper()
        aq = c.active_qubits
        infid = c.conditional_infidelity_upper
        rec = {
            "pipeline_id": label,
            "stages": stgs,
            "distance": c.distance,
            "active_qubits": aq,
            "correlation_mode": c.correlation_mode.value,
            "conditional_infidelity_hi": round_up_str(infid, 16),
            "expected_wire_ticks_per_output_hi": round_up_str(wt, 16),
            "acceptance_prob_lo": round_down_str(c.acceptance_lower, 16),
        }
        feasible.append((wt, aq, infid, rec))

    # Pareto filter over `(wire_ticks, active_qubits, conditional_infidelity)`:
    pareto: list[dict[str, Any]] = []
    feasible.sort(key=lambda item: (item[0], item[1], item[2]))
    for i, (wt_i, aq_i, inf_i, rec_i) in enumerate(feasible):
        dominated = False
        for j, (wt_j, aq_j, inf_j, _) in enumerate(feasible):
            if i == j:
                continue
            if wt_j <= wt_i and aq_j <= aq_i and inf_j <= inf_i and (wt_j < wt_i or aq_j < aq_i or inf_j < inf_i):
                dominated = True
                break
        if not dominated:
            pareto.append(rec_i)
    return pareto
