"""Live bridge to sibling BootLoops FTQC packages (`ftqc_audit`, `ftqc_zx`, `ftqc_extractors`).

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction
import importlib.util
import sys
from pathlib import Path
from typing import Any

from ..ir_kickmix import ANNOTATION_OPS, KickMixCircuit, Qubit
from ..kickmix_sim import z_pow_synthesis_cost


def _ensure_tools_on_syspath() -> None:
    tools_dir = str(Path(__file__).resolve().parents[2])
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)


def ftqc_audit_available() -> bool:
    _ensure_tools_on_syspath()
    return importlib.util.find_spec("ftqc_faults") is not None


def ftqc_zx_available() -> bool:
    _ensure_tools_on_syspath()
    return importlib.util.find_spec("ftqc_faults") is not None


def ftqc_compose_available() -> bool:
    _ensure_tools_on_syspath()
    return importlib.util.find_spec("ftqc_compose") is not None


def kmx_clifford_to_ftqc_audit_ticks(
    circuit: KickMixCircuit,
    *,
    prep_public_plus: bool = True,
) -> list[list[tuple[Any, ...]]]:
    """Lower a Clifford + `HMR` + `R` + single-qubit `Z if b` `KickMixCircuit` into
    `ftqc_audit.Circuit` ticks (`R`, `RX`, `H`, `S`, `SDG`, `CX`, `MX`, `Z_IF`).

    Raises `ValueError` if non-Clifford gates (`CCX`, `CCZ`, `CSWAP`, `T`) or multi-qubit
    classical feedback (`CZ ... if b`) appear without prior magic-state / gadget expansion.
    """
    pub = circuit.public_qubits() if circuit.registers else frozenset(range(circuit.num_qubits))
    ticks: list[list[tuple[Any, ...]]] = []

    # Tick 0: prepare all qubits (`!RX` for public superpositions or `!R` for scratch `|0>`)
    prep_tick: list[tuple[Any, ...]] = []
    for q in range(circuit.num_qubits):
        if q in pub and prep_public_plus:
            prep_tick.append(("!RX", q))
        else:
            prep_tick.append(("!R", q))
    if prep_tick:
        ticks.append(prep_tick)

    # Each bit holds a constant XOR an optional absolute MX record index.
    bit_values: dict[int, tuple[bool, int | None]] = {}
    meas_counter = 0

    for op in circuit.ops:
        kind = op.kind
        if kind in ANNOTATION_OPS:
            continue
        if kind in ("CCX", "CCZ", "CSWAP"):
            raise ValueError(f"Non-Clifford {kind} requires magic-state gadget expansion before ftqc_audit export")
        if kind in ("PUSH_CONDITION", "POP_CONDITION"):
            raise ValueError("ftqc_audit tick export supports direct `Z q if b` feedback; flatten condition stack first")

        if op.condition is not None:
            if kind != "Z":
                raise ValueError(f"ftqc_audit supports 1-qubit `Z_IF` feedback, got `{kind} ... if b{op.condition}`")
            inverted, record = bit_values.get(op.condition, (False, None))
            q = op.qubits[0]
            if inverted:
                ticks.extend([[("S", q)], [("S", q)]])
            if record is not None:
                ticks.append([("Z_IF", q, record)])
            continue

        if kind == "X":
            q = op.qubits[0]
            ticks.append([("H", q)])
            ticks.append([("S", q)])
            ticks.append([("S", q)])
            ticks.append([("H", q)])
        elif kind == "Z":
            q = op.qubits[0]
            ticks.append([("S", q)])
            ticks.append([("S", q)])
        elif kind == "Z_POW":
            assert op.angle is not None
            q = op.qubits[0]
            if op.angle == 0:
                continue
            if op.angle == Fraction(1, 2):
                ticks.append([("S", q)])
            elif op.angle == Fraction(-1, 2) or op.angle == Fraction(3, 2):
                ticks.append([("SDG", q)])
            elif op.angle == Fraction(1, 1) or op.angle == Fraction(-1, 1):
                ticks.append([("S", q)])
                ticks.append([("S", q)])
            else:
                raise ValueError(f"Non-Clifford Z_POW({op.angle}) cannot be exported directly to ftqc_audit")
        elif kind == "CX":
            ticks.append([("CX", op.qubits[0], op.qubits[1])])
        elif kind == "CZ":
            c1, c2 = op.qubits
            ticks.append([("H", c2)])
            ticks.append([("CX", c1, c2)])
            ticks.append([("H", c2)])
        elif kind == "SWAP":
            q1, q2 = op.qubits
            ticks.append([("CX", q1, q2)])
            ticks.append([("CX", q2, q1)])
            ticks.append([("CX", q1, q2)])
        elif kind == "R":
            ticks.append([("R", op.qubits[0])])
        elif kind == "HMR":
            q = op.qubits[0]
            b = op.bits[0]
            ticks.append([("MX", q)])
            bit_values[b] = (False, meas_counter)
            meas_counter += 1
            ticks.append([("R", q)])
        elif kind in ("BIT_STORE0", "BIT_STORE1"):
            bit_values[op.bits[0]] = (kind == "BIT_STORE1", None)
        elif kind == "BIT_INVERT":
            b = op.bits[0]
            inverted, record = bit_values.get(b, (False, None))
            bit_values[b] = (not inverted, record)
        elif kind == "NEG":
            continue  # An unconditional global phase has no stabilizer observable.
        else:
            raise ValueError(f"Unsupported ftqc_audit operation {kind!r}")

    return ticks


def audit_clifford_hmr_with_ftqc_audit(
    circuit: KickMixCircuit,
    *,
    name: str = "hmr_clifford_gadget",
) -> dict[str, Any]:
    """Build and verify a live `ftqc_audit` certificate for a Clifford+HMR `KickMixCircuit`."""
    _ensure_tools_on_syspath()
    from ftqc_faults import audit as ftqc_audit  # type: ignore[import-not-found]

    ticks = kmx_clifford_to_ftqc_audit_ticks(circuit, prep_public_plus=True)
    pub = sorted(circuit.public_qubits()) if circuit.registers else list(range(circuit.num_qubits))
    audit_circ = ftqc_audit.Circuit(ticks, n=circuit.num_qubits, outputs=pub, name=name)
    rec, _model = ftqc_audit.build(audit_circ, max_weight=1)
    ftqc_audit.verify(rec, full=True)
    return {
        "ftqc_audit_version": ftqc_audit.VERSION,
        "num_ticks": len(audit_circ.ticks),
        "num_measurements": len(audit_circ.meas),
        "outputs": audit_circ.outputs,
        "summary": ftqc_audit.summary(rec),
        "verified": True,
    }


def price_circuit_with_zx_contract(
    circuit: KickMixCircuit,
    zx_contract: dict[str, Any],
    *,
    t_factory_wire_ticks: int = 5000,
    t_factory_escape_hi: str | Fraction = "1e-12",
) -> dict[str, Any]:
    """Compute outward-rounded physical FTQC spacetime volume and logical failure budget
    by combining logical gate/HMR counts from `circuit` with a certified `ftqc_zx-contract-v1`
    EC/Clifford gadget contract and an `ftqc_compose` logical-T half-diamond bound (rational or Decimal str).
    """
    if zx_contract.get("schema") != "ftqc_zx-contract-v1":
        raise ValueError(f"Expected schema 'ftqc_zx-contract-v1', got {zx_contract.get('schema')!r}")
    if zx_contract.get("verdict") != "CONTRACT-CERTIFIED":
        raise ValueError(f"Refusing uncertified ZX contract verdict {zx_contract.get('verdict')!r}")

    ctx = Context(prec=28, rounding=ROUND_CEILING)
    stats = circuit.operation_statistics()
    t_info = z_pow_synthesis_cost(circuit)
    total_t = int(t_info["total_t_upper"])
    rx_depth = circuit.reaction_depth()
    cliff_rounds = stats.num_clifford_2q + rx_depth + 1

    wt_per_ec = Decimal(str(zx_contract["enclosures"]["expected_wire_ticks_per_shot_hi"]))
    esc_per_ec = Decimal(str(zx_contract["enclosures"]["conditional_escape_hi"]))
    if isinstance(t_factory_escape_hi, Fraction):
        esc_per_t = ctx.divide(Decimal(t_factory_escape_hi.numerator), Decimal(t_factory_escape_hi.denominator))
    elif "/" in str(t_factory_escape_hi):
        n_s, d_s = str(t_factory_escape_hi).split("/", 1)
        esc_per_t = ctx.divide(Decimal(int(n_s)), Decimal(int(d_s)))
    else:
        esc_per_t = Decimal(str(t_factory_escape_hi))

    if type(t_factory_wire_ticks) is not int or t_factory_wire_ticks < 0:
        raise ValueError("Factory wire-ticks must be a nonnegative integer")
    if not wt_per_ec.is_finite() or wt_per_ec < 0:
        raise ValueError("Expected EC wire-ticks must be finite and nonnegative")
    for rate in (esc_per_ec, esc_per_t):
        if not rate.is_finite() or not 0 <= rate <= 1:
            raise ValueError("Escape bounds must lie in [0,1]")

    ec_volume = ctx.multiply(Decimal(cliff_rounds + total_t), wt_per_ec)
    factory_volume = ctx.multiply(Decimal(total_t), Decimal(t_factory_wire_ticks))
    total_volume = ctx.add(ec_volume, factory_volume)

    total_escape = ctx.add(
        ctx.multiply(Decimal(cliff_rounds + total_t), esc_per_ec),
        ctx.multiply(Decimal(total_t), esc_per_t),
    )
    if total_escape > Decimal(1):
        total_escape = Decimal(1)

    return {
        "zx_contract_circuit": zx_contract.get("circuit_name"),
        "codespace": zx_contract.get("codespace"),
        "logical_t_upper": total_t,
        "reaction_depth": rx_depth,
        "ec_gadget_invocations": cliff_rounds + total_t,
        "ec_wire_ticks_upper": str(ec_volume),
        "t_factory_wire_ticks_upper": str(factory_volume),
        "total_wire_ticks_upper": str(total_volume),
        "total_logical_escape_upper": str(total_escape),
    }


def estimate_reaction_decoder_latency(
    circuit_or_depth: KickMixCircuit | int,
    schedule_metrics: dict[str, Any],
    *,
    decoder_stall_ticks_per_reaction: int | Fraction | str = 0,
    shuttle_ticks_per_reaction: int = 0,
    installed_action_per_reaction: dict[str, Any] | None = None,
    surviving_logicals: int = 1,
    bare_code_k: int | None = None,
    layout_name: str = "planar_local",
) -> dict[str, Any]:
    """Compute exact / outward-rounded reaction-depth latency and decoder-stall spacetime budget.

    Distinguishes bare code rate `k/n` from `physical_qubits_per_surviving_logical`
    during active surgery/measurement schedules, and charges data-qubit idling across
    all measured classical decoder stall ticks. Zero binary floating-point values are used.
    """
    if isinstance(circuit_or_depth, KickMixCircuit):
        rx_depth = int(circuit_or_depth.reaction_depth())
    elif type(circuit_or_depth) is int and circuit_or_depth >= 0:
        rx_depth = circuit_or_depth
    else:
        raise ValueError("circuit_or_depth must be a KickMixCircuit or nonnegative int")

    for req_k in ("num_data_qubits", "peak_ancilla_qubits", "rounds_per_reaction", "ticks_per_round"):
        val = schedule_metrics.get(req_k)
        if type(val) is not int or val < 0:
            raise ValueError(f"schedule_metrics[{req_k!r}] must be a nonnegative int")

    n_data = schedule_metrics["num_data_qubits"]
    n_anc = schedule_metrics["peak_ancilla_qubits"]
    rounds_per_rx = schedule_metrics["rounds_per_reaction"]
    ticks_per_rnd = schedule_metrics["ticks_per_round"]
    couplers = schedule_metrics.get("installed_couplers", 0)

    if n_data < 1:
        raise ValueError("num_data_qubits must be >= 1")
    if type(couplers) is not int or couplers < 0:
        raise ValueError("installed_couplers must be a nonnegative int")
    if type(shuttle_ticks_per_reaction) is not int or shuttle_ticks_per_reaction < 0:
        raise ValueError("shuttle_ticks_per_reaction must be a nonnegative int")
    if type(surviving_logicals) is not int or surviving_logicals < 1:
        raise ValueError("surviving_logicals must be a positive int")
    if bare_code_k is not None and (type(bare_code_k) is not int or bare_code_k < 0):
        raise ValueError("bare_code_k must be a nonnegative int or None")

    if type(decoder_stall_ticks_per_reaction) is int:
        stall_per_rx = Fraction(decoder_stall_ticks_per_reaction, 1)
    elif isinstance(decoder_stall_ticks_per_reaction, Fraction):
        stall_per_rx = decoder_stall_ticks_per_reaction
    elif isinstance(decoder_stall_ticks_per_reaction, str):
        stall_per_rx = Fraction(decoder_stall_ticks_per_reaction)
    else:
        raise ValueError("decoder_stall_ticks_per_reaction must be int, Fraction, or exact string")

    if stall_per_rx < 0:
        raise ValueError("decoder_stall_ticks_per_reaction must be nonnegative")

    if installed_action_per_reaction is not None:
        from ftqc_gadgets.installed_action import analyze
        route = analyze(installed_action_per_reaction)
        if route["n_data"] != n_data or route["width"] > n_data + n_anc:
            raise ValueError("installed action exceeds the reaction allocation")
        if shuttle_ticks_per_reaction:
            raise ValueError("supply the installed route or manual shuttle ticks")
        if route["installed_couplers"] > couplers:
            raise ValueError("installed action exceeds the reaction coupler budget")
        shuttle_ticks_per_reaction = route["ticks"]

    ctx = Context(prec=28, rounding=ROUND_CEILING)
    peak_total = n_data + n_anc
    active_rounds = rx_depth * rounds_per_rx
    active_wall_ticks = active_rounds * ticks_per_rnd
    shuttle_wall_ticks = rx_depth * shuttle_ticks_per_reaction

    total_stall_frac = Fraction(rx_depth, 1) * stall_per_rx
    stall_ticks_ceil = (total_stall_frac.numerator + total_stall_frac.denominator - 1) // total_stall_frac.denominator
    stall_ticks_dec = ctx.divide(Decimal(total_stall_frac.numerator), Decimal(total_stall_frac.denominator))

    total_wall_ticks_ceil = active_wall_ticks + shuttle_wall_ticks + stall_ticks_ceil
    active_qubit_ticks = (active_wall_ticks + shuttle_wall_ticks) * peak_total
    stall_idle_qubit_ticks_ceil = stall_ticks_ceil * n_data
    total_qubit_ticks_ceil = active_qubit_ticks + stall_idle_qubit_ticks_ceil

    return {
        "layout": str(layout_name),
        "reaction_depth": rx_depth,
        "rounds_per_reaction": rounds_per_rx,
        "active_rounds": active_rounds,
        "ticks_per_round": ticks_per_rnd,
        "active_wall_ticks": active_wall_ticks,
        "shuttle_wall_ticks": shuttle_wall_ticks,
        "decoder_stall_ticks_exact": f"{total_stall_frac.numerator}/{total_stall_frac.denominator}",
        "decoder_stall_ticks_upper": str(stall_ticks_dec),
        "decoder_stall_ticks_ceil": stall_ticks_ceil,
        "total_wall_ticks_ceil": total_wall_ticks_ceil,
        "installed_couplers": couplers,
        "peak_data_qubits": n_data,
        "peak_ancilla_qubits": n_anc,
        "peak_total_qubits": peak_total,
        "surviving_logicals": surviving_logicals,
        "bare_code_rate": f"{bare_code_k}/{n_data}" if bare_code_k is not None else None,
        "physical_qubits_per_surviving_logical": f"{peak_total}/{surviving_logicals}",
        "active_qubit_ticks": active_qubit_ticks,
        "stall_idle_qubit_ticks_ceil": stall_idle_qubit_ticks_ceil,
        "total_qubit_ticks_ceil": total_qubit_ticks_ceil,
    }
