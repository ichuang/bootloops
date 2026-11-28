"""Reference gadgets: positive and negative controls with known answers.

cat:N[,verify]            cat state by a CX chain; `verify` adds one ancilla that
                          measures Z_0 Z_{N-1} (catches the weight-2 hook of the chain)
css_prep:CODE             non-fault-tolerant |0_L> of a ftqc_codes CSS code: |+> on an
                          information set, CX fan-out (CODE is a ftqc_codes spec, e.g. golay)
surface:D,R[,BASIS[,bad]] rotated surface code memory, R rounds of syndrome extraction,
                          BASIS z (default) or x; `bad` swaps the CX orders of the X and Z
                          checks so that hook errors line up with the logical operators
"""
import os
import sys

from .circuit import Circuit, CircuitError


def with_idles(layers, n, first=None, last=None):
    """Insert I on every qubit that is alive (between its first and last use) and
    untouched in a tick."""
    used = [set(q for op in t for q in op[1:]) for t in layers]
    lo, hi = {}, {}
    for k, s in enumerate(used):
        for q in s:
            lo.setdefault(q, k)
            hi[q] = k
    out = []
    for k, t in enumerate(layers):
        row = list(t)
        for q in range(n):
            if q in lo and lo[q] < k < hi[q] and q not in used[k]:
                row.append(("I", q))
        out.append(row)
    return out


def cat(n, verify=False):
    n = int(n)
    if n < 2:
        raise CircuitError("cat needs n >= 2")
    layers = [[("RX", 0)] + [("R", q) for q in range(1, n)]]
    layers += [[("CX", q, q + 1)] for q in range(n - 1)]
    total = n
    if verify:
        a = n
        total = n + 1
        layers[0].append(("R", a))
        layers += [[("CX", 0, a)], [("CX", n - 1, a)], [("M", a)]]
    c = Circuit(with_idles(layers, total), total, outputs=list(range(n)),
                name=f"cat{n}" + ("_verified" if verify else ""))
    return c


def _ftqc_codes():
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if here not in sys.path:
        sys.path.insert(0, here)
    import ftqc_codes
    from ftqc_codes import gf2
    return ftqc_codes, gf2


def css_prep(spec):
    """|0_L> of a CSS code: the +1 eigenstate of the X stabilizers and of all Z-type
    operators in ker(H_X) (stabilizers and logical Z).  Circuit: put the rows of H_X in
    systematic form, prepare |+> on the pivots and |0> elsewhere, fan out with CX."""
    ftqc_codes, gf2 = _ftqc_codes()
    code = ftqc_codes.load(spec) if isinstance(spec, str) else spec
    n = code.n
    rows, piv = gf2.systematic(code.hx, n, list(range(n)))
    layers = [[("RX", q) if q in piv else ("R", q) for q in range(n)]]
    edges = [(p, t) for r, p in zip(rows, piv) for t in gf2.supp(r) if t != p]
    # greedy layering, one gate per qubit per tick
    busy = {}
    cx = []
    for c, t in edges:
        k = 0
        while (k, c) in busy or (k, t) in busy:
            k += 1
        busy[k, c] = busy[k, t] = True
        while len(cx) <= k:
            cx.append([])
        cx[k].append(("CX", c, t))
    return Circuit(with_idles(layers + cx, n), n, outputs=list(range(n)),
                   name=f"zero_prep({code.name})")


def surface(d, rounds, basis="z", bad=False):
    """Rotated surface code memory experiment.  Data qubits r*d + c; one ancilla per
    check; each round: prepare ancillas, four CX ticks, measure ancillas."""
    ftqc_codes, gf2 = _ftqc_codes()
    d, rounds = int(d), int(rounds)
    basis = basis.lower()
    code = ftqc_codes.load(f"rotated_surface:{d}")
    n = d * d
    row0 = gf2.from_supp(range(d))                       # first row of the grid
    # a row is an X-type logical iff it commutes with every Z check
    row_is_lx = all(gf2.dot(row0, h) == 0 for h in code.hz)
    # corners of a plaquette in the two standard orders
    Zshape, Nshape = ("NW", "NE", "SW", "SE"), ("NW", "SW", "NE", "SE")
    # An X-check ancilla fault spreads X to the last two data qubits of its order:
    # Z-shape ends on a horizontal pair, N-shape on a vertical pair.  The good choice
    # ends perpendicular to the logical operator of the same type.
    x_order = Nshape if row_is_lx else Zshape            # logical X horizontal -> end vertical
    z_order = Zshape if row_is_lx else Nshape            # logical Z is the other direction
    if bad:
        x_order, z_order = z_order, x_order
    from ftqc_codes.codes import rotated_plaquettes
    checks = [(axis, cells) for axis, _, cells in rotated_plaquettes(d)]
    anc = {i: n + i for i in range(len(checks))}
    total = n + len(checks)
    data_prep, data_meas = ("R", "M") if basis == "z" else ("RX", "MX")
    layers = [[(data_prep, q) for q in range(n)]]
    for _ in range(rounds):
        layers.append([("RX" if t == "X" else "R", anc[i]) for i, (t, _) in enumerate(checks)])
        for step in range(4):
            tick = []
            for i, (t, cells) in enumerate(checks):
                key = (x_order if t == "X" else z_order)[step]
                if key in cells:
                    tick.append(("CX", anc[i], cells[key]) if t == "X" else ("CX", cells[key], anc[i]))
            layers.append(tick)
        layers.append([("MX" if t == "X" else "M", anc[i]) for i, (t, _) in enumerate(checks)])
    layers.append([(data_meas, q) for q in range(n)])
    circ = Circuit(with_idles(layers, total), total, name="tmp")
    # logical observable: the final data measurements along a logical operator of the basis
    nm = len(circ.meas)
    final = {q: nm - n + q for q in range(n)}
    col0 = gf2.from_supp(r * d for r in range(d))
    if basis == "z":
        string = col0 if row_is_lx else row0             # a Z-type logical
    else:
        string = row0 if row_is_lx else col0
    obs = [[final[q] for q in gf2.supp(string)]]
    # detectors, declared explicitly: checks of the memory basis in the first round,
    # every check against itself one round earlier, and at the end the data parity of
    # each basis-type check against its last ancilla outcome
    C = len(checks)
    kind = "Z" if basis == "z" else "X"
    dets = [[i] for i, (t, _) in enumerate(checks) if t == kind]
    for r in range(1, rounds):
        dets += [[(r - 1) * C + i, r * C + i] for i in range(C)]
    dets += [[(rounds - 1) * C + i] + [final[q] for q in cells.values()]
             for i, (t, cells) in enumerate(checks) if t == kind]
    return Circuit(with_idles(layers, total), total, detectors=dets, observables=obs, outputs=[],
                   name=f"surface_d{d}_r{rounds}_{basis}" + ("_bad" if bad else ""))


def load(spec):
    """'cat:4,verify', 'css_prep:golay', 'surface:3,2', 'surface:3,2,x,bad' or a JSON path."""
    if spec.endswith(".json"):
        from . import circuit
        return circuit.load(spec)
    fam, _, arg = spec.partition(":")
    parts = [p for p in arg.split(",") if p]
    try:
        if fam == "cat":
            return cat(parts[0], "verify" in parts[1:])
        if fam == "css_prep":
            return css_prep(arg)
        if fam == "surface":
            return surface(parts[0], parts[1], parts[2] if len(parts) > 2 else "z", "bad" in parts)
    except IndexError:
        raise CircuitError(f"missing arguments for {fam}") from None
    raise CircuitError(f"unknown gadget {fam!r}; see ftqc_faults.audit.library")
