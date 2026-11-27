"""Exact CSS source isometries for simultaneous injection (Liu et al., sections III–IV)."""

from ftqc_codes.gf2 import Echelon, systematic
from .algebra import bits, parity, transpose
from .block_reading import memory, row_solver


def site_list(sites, n):
    if (
        not isinstance(sites, list)
        or any(type(q) is not int or not 0 <= q < n for q in sites)
        or len(set(sites)) != len(sites)
    ):
        raise ValueError("expected distinct physical site indices")
    return sites


def mask(sites):
    return sum(1 << q for q in sites)


def logicals(s):
    n, hx, hz = memory(s["memory"])
    pairs = s["logicals"]
    if not isinstance(pairs, list) or not pairs:
        raise ValueError("supply at least one signed CSS logical pair")
    for p in pairs:
        if (
            not isinstance(p, list)
            or len(p) != 4
            or any(type(v) is not int for v in p)
            or min(p[:2]) <= 0
            or max(p[:2]).bit_length() > n
            or not set(p[2:]) <= {0, 1}
        ):
            raise ValueError(
                "logical pair requires [X mask, Z mask, X sign bit, Z sign bit]"
            )
        if any(parity(p[0], v) for v in hz) or any(parity(p[1], v) for v in hx):
            raise ValueError("representative anticommutes with a memory check")
    if any(
        parity(p[0], q[1]) != (i == j)
        for i, p in enumerate(pairs)
        for j, q in enumerate(pairs)
    ):
        raise ValueError("representatives are not paired logical coordinates")
    return n, hx, hz, pairs


def clean(word, basis, pivots):
    for row, q in zip(basis, pivots):
        if word >> q & 1:
            word ^= row
    return word


def carrier_encoder(xs, sites):
    """Reduce the X images by physical CNOT/SWAP conjugations and reverse the circuit."""
    current = list(xs)
    undo = []

    def cx(a, b):
        for i, v in enumerate(current):
            if v >> a & 1:
                current[i] ^= 1 << b
        undo.append(["CX", a, b])

    for i, q in enumerate(sites):
        if not current[i] >> q & 1:
            other = next((r for r in sites[i + 1 :] if current[i] >> r & 1), None)
            if other is None:
                raise ValueError("carrier action is singular")
            cx(q, other)
            cx(other, q)
            cx(q, other)
        for r in list(bits(current[i] ^ (1 << q))):
            cx(q, r)
    return list(reversed(undo))


def layout(s, n, hx, hz, pairs):
    mode = s["method"]
    if mode == "general":
        if "preparation" in s or len(pairs) != n - Echelon(hx).rank - Echelon(hz).rank:
            raise ValueError(
                "general injection requires the complete ordered logical basis"
            )
        xb, xp = systematic(hx, n, range(n))
        zb, zp = systematic(hz, n, [q for q in range(n) if q not in xp])
        sites = [q for q in range(n) if q not in xp and q not in zp]
        cleaned = [
            [clean(x, xb, xp), clean(z, zb, zp), sx, sz] for x, z, sx, sz in pairs
        ]
        gates = carrier_encoder([p[0] & mask(sites) for p in cleaned], sites)
        return sites, zp, xp, [], cleaned, gates
    if mode != "selected":
        raise ValueError("unknown injection method")
    prep = s["preparation"]
    if set(prep) != {"carriers", "plus", "zero", "bells"}:
        raise ValueError("provide a complete preparation")
    sites, plus, zero = (site_list(prep[k], n) for k in ("carriers", "plus", "zero"))
    bells = prep["bells"]
    if any(len(site_list(b, n)) != 2 for b in bells):
        raise ValueError("a Bell pair has two endpoints")
    allsites = sites + plus + zero + [q for b in bells for q in b]
    if len(sites) != len(pairs) or sorted(allsites) != list(range(n)):
        raise ValueError("preparation must partition the memory")
    return sites, plus, zero, bells, pairs, []


def subsystem_error(word, single_sites, bells):
    if word & mask(single_sites):
        return "incompatible single-qubit preparation"
    if any(parity(word, mask(pair)) for pair in bells):
        return "logical acts on only one Bell endpoint"
    return None


def frame(checks, logical_rows, n):
    independent = Echelon(checks).picked
    measurements = [checks[i] for i in independent] + logical_rows
    solve = row_solver(transpose(measurements, n))
    return [[idx, solve(1 << i)] for i, idx in enumerate(independent)]


def analyze(s):
    n, hx, hz, pairs = logicals(s)
    sites, plus, zero, bells, reps, gates = layout(s, n, hx, hz, pairs)
    # Prepend single-site sign corrections in carrier coordinates.
    signs = [
        op
        for q, (_, _, sx, sz) in zip(sites, reps)
        for op, b in ((["Z", q], sx), (["X", q], sz))
        if b
    ]
    gates = signs + gates

    def pull(word, axis, sign):
        for gate, *qs in reversed(gates):
            if gate == "CX":
                c, t = qs if axis == "X" else qs[::-1]
                if word >> c & 1:
                    word ^= 1 << t
            elif gate != axis:
                sign ^= (word >> qs[0]) & 1
        return word, sign

    fixed = {}
    for axis, checks, incompatible in (("X", hx, zero), ("Z", hz, plus)):
        fixed[axis] = [
            i
            for i, w in enumerate(checks)
            if not w & mask(sites) and not subsystem_error(w, incompatible, bells)
        ]
    for i, (q, p) in enumerate(zip(sites, reps)):
        for axis, w, sign, incompatible in (
            ("X", p[0], p[2], zero),
            ("Z", p[1], p[3], plus),
        ):
            fault = subsystem_error(w, incompatible, bells)
            reduced, phase = pull(w & mask(sites), axis, sign)
            if fault or (reduced, phase) != (1 << q, 0):
                return dict(
                    status="counterexample",
                    logical=i,
                    axis=axis,
                    reason=fault or "incorrect signed carrier action",
                )
    return dict(
        status="proved",
        scope="ideal source isometry on the specified logical subsystem",
        carriers=sites,
        plus=plus,
        zero=zero,
        bells=bells,
        carrier_gates=gates,
        representatives=reps,
        fixed_checks=fixed,
        syndrome_frame={
            "Z": frame(hx, [p[0] for p in reps], n),
            "X": frame(hz, [p[1] for p in reps], n),
        },
    )
