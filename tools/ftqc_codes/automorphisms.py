"""Signed Clifford actions on a specified CSS memory and ordered logical basis.

A step is either a push permutation (old site -> new site), or one disjoint
layer of local or paired Clifford gates. Layers are applied in list order.
"""
from . import basis, certify, codes, interchange, symmetry

SCHEMA = "ftqc_codes-action-1"
LOCAL = {"H", "S", "SDG", "X", "Y", "Z"}
PAIRED = {"CX", "CZ", "SWAP"}


def validate(n, steps):
    if not isinstance(steps, list):
        raise codes.SpecError("steps must be a list")
    for step in steps:
        if isinstance(step, dict) and set(step) == {"perm"}:
            p = step["perm"]
            if not isinstance(p, list) or any(type(q) is not int for q in p):
                raise codes.SpecError("expected integer permutation")
            interchange.validate_permutation(p, n, "action")
            continue
        if not isinstance(step, dict) or set(step) != {"layer"} or not isinstance(step["layer"], list):
            raise codes.SpecError("expected perm or layer")
        used = set()
        for op in step["layer"]:
            if not isinstance(op, list) or not op or not isinstance(op[0], str):
                raise codes.SpecError("invalid Clifford operation")
            gate, *sites = op
            if (gate not in LOCAL | PAIRED or len(sites) != 1 + (gate in PAIRED)
                    or any(type(q) is not int or not 0 <= q < n for q in sites)
                    or len(set(sites)) != len(sites) or used.intersection(sites)):
                raise codes.SpecError(f"invalid or overlapping operation {op!r}")
            used.update(sites)


def conjugate(pauli, steps):
    """Apply validated steps to a Hermitian Pauli (x,z,negative-sign bit)."""
    x, z, sign = pauli
    for step in steps:
        if "perm" in step:
            x = interchange.transport_word(x, step["perm"])
            z = interchange.transport_word(z, step["perm"])
            continue
        for gate, q, *rest in step["layer"]:
            a, b = (x >> q) & 1, (z >> q) & 1
            if gate == "H":
                sign ^= a & b
                x ^= (a ^ b) << q
                z ^= (a ^ b) << q
            elif gate in ("S", "SDG"):
                sign ^= a & (b ^ (gate == "SDG"))
                z ^= a << q
            elif gate in ("X", "Y", "Z"):
                sign ^= (a if gate in ("Y", "Z") else 0) ^ (b if gate in ("X", "Y") else 0)
            else:
                t = rest[0]
                c, d = (x >> t) & 1, (z >> t) & 1
                if gate == "CX":
                    sign ^= a & d & (1 ^ b ^ c)
                    x ^= a << t
                    z ^= d << q
                elif gate == "CZ":
                    sign ^= a & c & (b ^ d)
                    z ^= (a << t) ^ (c << q)
                else:
                    pair = (1 << q) | (1 << t)
                    if a ^ c:
                        x ^= pair
                    if b ^ d:
                        z ^= pair
    return x, z, sign


def build(code, paired, steps):
    basis.verify_basis_receipt(code, paired)
    validate(code.n, steps)
    checks = basis._code_signed_checks(code)
    group = basis._SignedStabilizerGroup(code.n, checks)
    for check in checks:
        if not group.contains_exact(conjugate(check, steps)):
            raise certify.Refusal("physical Clifford does not preserve signed stabilizers")
    logical = symmetry.logical_action(code, paired, lambda p: conjugate(p, steps))
    return dict(schema=SCHEMA, basis=paired, steps=steps, **logical)


def verify(code, receipt):
    if not isinstance(receipt, dict) or receipt.get("schema") != SCHEMA:
        raise certify.Refusal("expected a code action receipt")
    if build(code, receipt.get("basis"), receipt.get("steps")) != receipt:
        raise certify.Refusal("signed code action replay mismatch")
    return receipt


def logical_image(action, pauli):
    """Conjugate a signed logical Pauli using the certified generator images."""
    x, z, sign = pauli
    ops = []
    for mask, name in ((x, "logical_x"), (z, "logical_z")):
        for j, image in enumerate(action[name]):
            if mask >> j & 1:
                ops.append((image["x_mask"], image["z_mask"], image["sign"]))
    return basis._mul_chain_signed(ops, sign, (x & z).bit_count())

