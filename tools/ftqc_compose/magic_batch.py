"""Ideal diagonal injection with a shared source-acceptance event and ordered ports."""

from fractions import Fraction
from itertools import combinations

from ftqc_magic.delivery import read_delivery
from ftqc_magic.ir import PhasePolynomialTarget
from .supports import fields


def injection(name):
    """Replay the CX(data -> resource), MZ and diagonal Clifford on each transcript."""
    target = PhasePolynomialTarget.named(name)
    if name not in ("T", "CCZ"):
        raise ValueError("unsupported injection target")
    table = target.phases
    k = target.k
    dim = 1 << k
    branches = []
    for measured in range(dim):
        correction = [(table[x] - table[x ^ measured]) % 8 for x in range(dim)]
        constant = correction[0]
        linear = [(correction[1 << i] - constant) % 8 for i in range(k)]
        quadratic = [
            (
                i,
                j,
                (correction[(1 << i) | (1 << j)] - constant - linear[i] - linear[j])
                % 8,
            )
            for i, j in combinations(range(k), 2)
        ]
        if any(c % 2 for c in linear) or any(c not in (0, 4) for _, _, c in quadratic):
            raise ValueError("non-Clifford feed-forward")
        for x in range(dim):
            phase = (
                constant
                + sum(c for i, c in enumerate(linear) if x >> i & 1)
                + sum(c for i, j, c in quadratic if x >> i & 1 and x >> j & 1)
            ) % 8
            # Check every branch, computational basis amplitude and resource Z.
            for mask in range(dim):
                actual = (
                    table[x ^ measured]
                    + phase
                    + 4 * ((x ^ measured) & mask).bit_count()
                ) % 8
                expected = (
                    table[x]
                    + 4 * (x & mask).bit_count()
                    + 4 * (measured & mask).bit_count()
                ) % 8
                if actual != expected:
                    raise ValueError("injection residual or branch phase mismatch")
        branches.append(
            dict(
                measured=measured,
                probability=str(Fraction(1, dim)),
                global_phase=constant,
                single_wire_S_powers=[v // 2 for v in linear],
                cz=[[i, j] for i, j, c in quadratic if c],
            )
        )
    return dict(
        target=name,
        branches=branches,
        operations=[["CX", f"data.{i}", f"resource.{i}"] for i in range(k)]
        + [["MZ", f"resource.{i}"] for i in range(k)],
        residual="same ordered Z mask on data",
    )


def consume_batch(spec):
    fields(
        spec,
        "kind",
        "factory",
        "factory_receipt",
        "batch",
        "ports",
        *(k for k in ("catalyst", "twirl") if k in spec),
    )
    src = read_delivery(spec["factory"], spec["factory_receipt"])
    if not isinstance(spec["batch"], str) or not spec["batch"].strip():
        raise ValueError("missing batch identity")
    ports = spec["ports"]
    if len(ports) != len(src["groups"]):
        raise ValueError("need one consumer per batch member")
    seen = set()
    for i, (port, wires) in enumerate(zip(ports, src["groups"], strict=True)):
        fields(port, "member", "wires", "data")
        if (
            type(port["member"]) is not int
            or port["member"] != i
            or port["wires"] != wires
            or len(port["data"]) != len(wires)
        ):
            raise ValueError("batch member or wire order mismatch")
        for q in port["data"]:
            if not isinstance(q, str) or not q or q in seen:
                raise ValueError("duplicate or missing data port")
            seen.add(q)
    terminal = src["correlation"] == "CATALYST_POISONED"
    if terminal:
        c = spec.get("catalyst", {})
        if (
            set(c) != {"id", "history", "retired_after"}
            or not isinstance(c["id"], str)
            or not c["id"]
            or c["history"] != []
            or c["retired_after"] != spec["batch"]
        ):
            raise ValueError(
                "terminal catalysis needs a fresh identified catalyst retired after this batch"
            )
    elif "catalyst" in spec:
        raise ValueError("unexpected catalyst")
    inst = injection(src["state"])
    # The instrument above proves that every mask reaches the named data
    # ports unchanged. Keep a single common atom for all members it affects.
    atoms = None
    accepted = src["accepted_z_mass"]
    conditional = src["conditional_z_mass"]
    if "twirl" in spec:
        if (
            spec["twirl"] != "joint_3twirl"
            or src["state"] != "CCZ"
            or conditional is None
        ):
            raise ValueError("joint CCZ twirl needs a verified full Z-mask law")
        from .magic_correlations import ccz_twirl

        accepted = ccz_twirl(accepted, len(ports))
        conditional = {
            m: str(Fraction(v) / Fraction(src["acceptance_lower"]))
            for m, v in accepted.items()
        }
    if conditional is not None:
        atoms = []
        for key, p in conditional.items():
            mask = int(key)
            hits = [
                q
                for port in ports
                for w, q in zip(port["wires"], port["data"], strict=True)
                if mask >> w & 1
            ]
            atoms.append(
                dict(
                    batch=spec["batch"],
                    mask=mask,
                    data_Z=hits,
                    conditional_mass=p,
                    accepted_mass=accepted[key],
                )
            )
    return dict(
        status="proved",
        source=src,
        batch=spec["batch"],
        ports=ports,
        twirl=spec.get("twirl"),
        instrument=inst,
        conditional_joint_atoms=atoms,
        conditional_joint_trace_distance_upper=src["conditional_trace_distance_upper"],
        catalyst=spec.get("catalyst"),
        scope="ideal logical consumers; joint source independent of data, references and subsequent noise",
        conditioning="one acceptance event for the whole batch",
    )
