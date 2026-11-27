"""Replay a producer before exposing its ordered, acceptance-conditioned output."""

from fractions import Fraction
from math import isqrt

from .certify import evaluate_source, verify_receipt
from .ir import PhasePolynomialTarget


def sqrt_upper(q):
    scale = 10**18
    n = (q.numerator * scale**2 + q.denominator - 1) // q.denominator
    r = isqrt(n)
    return Fraction(r + (r * r < n), scale)


def read_delivery(spec, receipt):
    identity = verify_receipt(spec, receipt)
    kind, c, _, _ = evaluate_source(spec)
    mode = str(spec.get("mode", "postselect")).lower().strip()
    correcting = mode in ("ec", "correct", "correcting")
    if kind == "factory" and not correcting and mode != "postselect":
        raise ValueError("unknown factory acceptance rule")
    t = c.target
    groups = []
    for name, width in (("T", 1), ("CCZ", 3)):
        if t.k % width:
            continue
        count = t.k // width
        expected = PhasePolynomialTarget.named(f"{name}^{count}")
        if (t.modulus, t.phases) == (expected.modulus, expected.phases):
            groups = [list(range(i, i + width)) for i in range(0, t.k, width)]
            break
    if not groups:
        name, groups = "phase-polynomial", [list(range(t.k))]
    exact = kind == "factory" and (
        correcting or c.metadata.get("exact_probability_dp") is True
    )
    acc = c.acceptance_lower
    if not 0 < acc <= 1:
        raise ValueError("delivery needs a positive acceptance lower bound")
    joint = conditional = None
    delta = sqrt_upper(c.conditional_infidelity_upper)
    if exact:
        masses = dict(c.joint_z_error_upper)
        masses[0] = acc - sum(masses.values())
        if any(not 0 <= m < 1 << t.k or v < 0 for m, v in masses.items()):
            raise ValueError("invalid accepted mask distribution")
        delta = 1 - masses[0] / acc
        if delta != c.conditional_infidelity_upper:
            raise ValueError("joint law disagrees with factory infidelity")
        joint = {str(m): str(v) for m, v in sorted(masses.items())}
        conditional = {m: str(Fraction(v) / acc) for m, v in joint.items()}
    return dict(
        schema="ftqc-magic-delivery-1",
        producer=identity,
        state=name,
        modulus=t.modulus,
        phases=list(t.phases),
        groups=groups,
        wire_order="little-endian",
        frame="canonical-positive",
        correlation=c.correlation_mode.value,
        distance=c.distance,
        leading_count=c.leading_order_coeff,
        acceptance_rule=mode if kind == "factory" else kind,
        acceptance_lower=str(acc),
        acceptance_exact=bool(exact),
        conditional_infidelity_upper=str(c.conditional_infidelity_upper),
        conditional_trace_distance_upper=str(delta),
        accepted_bad_mass_upper=str(c.unnormalized_failure_upper),
        law="exact-accepted-Z" if exact else "state-bound-only",
        accepted_z_mass=joint,
        conditional_z_mass=conditional,
        premises=list(c.premises),
        active_qubits=c.active_qubits,
        expected_wire_ticks_per_output_upper=str(
            c.expected_wire_ticks_per_output_upper()
        ),
    )
