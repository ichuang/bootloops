"""Construct finite homomorphic and mapping-cone measurements in a local Pauli basis.

Rows of A label auxiliary X checks on vertex coordinates. F assigns each
vertex a data Z mask after local basis rotation. Products in ker(A) are the
unmasked outcomes. The instrument checker verifies their action on the code.
"""

from ftqc_codes.gf2 import Echelon, nullspace, supp, systematic
from ftqc_gadgets.surgery import pauli


def combine(rows, mask):
    out = 0
    for i in supp(mask):
        out ^= rows[i]
    return out


def independent_section(rows, width):
    """A right inverse on independent equations, in original row order."""
    basis = Echelon(rows)
    echelon = {}
    for i in basis.picked:
        v, coeff = rows[i], 1 << i
        for pivot in sorted(echelon, reverse=True):
            if v >> pivot & 1:
                r, c = echelon[pivot]
                v ^= r
                coeff ^= c
        echelon[v.bit_length() - 1] = v, coeff
    masks = [0] * width
    for pivot in sorted(echelon):
        v, c = echelon[pivot]
        masks[pivot] = c ^ combine(masks, v ^ (1 << pivot))
    return masks


def construct(n, checks, axes, ports, auxiliary, logicals, method, rounds=1):
    """logicals contains (vertex parity mask, public sign bit) pairs."""
    if type(n) != int or n < 1 or len(axes) != n or set(axes) - set("XYZ"):
        raise ValueError("expected a local measurement basis per data site")
    v = len(ports)
    if not v or any(type(p) != int or not 0 <= p < 1 << n for p in ports):
        raise ValueError("invalid vertex port")
    if any(type(a) != int or not 0 < a < 1 << v for a in auxiliary):
        raise ValueError("invalid auxiliary X check")
    if type(rounds) != int or rounds < 1 or method not in ("homomorphic", "cone"):
        raise ValueError("unknown construction or repetition count")
    for check in checks:
        pauli(check, n)
    targets = []
    for mask, sign in logicals:
        if (
            type(mask) != int
            or not 0 < mask < 1 << v
            or type(sign) != int
            or sign not in (0, 1)
        ):
            raise ValueError("invalid signed logical selector")
        if any((a & mask).bit_count() % 2 for a in auxiliary):
            raise ValueError("public parity depends on a private gauge outcome")
        data = combine(ports, mask)
        targets.append(
            ("-" if sign else "+")
            + "".join(axes[q] if data >> q & 1 else "I" for q in range(n))
        )
    ticks, detectors = [], []

    def op(name, *sites):
        ticks.append([[name, *sites]])

    used = {q for port in ports for q in supp(port)}

    def rotate(inverse):
        for q in sorted(used):
            if axes[q] == "Y" and not inverse:
                op("SDG", q)
            if axes[q] != "Z":
                op("H", q)
            if axes[q] == "Y" and inverse:
                op("S", q)

    rotate(False)
    records = 0
    prior = None
    first = []
    for _ in range(rounds):
        vertex = list(range(n, n + v))
        edges = list(range(n + v, n + v + len(auxiliary)))
        if method == "homomorphic":
            rows, pivots = systematic(auxiliary, v, list(range(v)))
            for j, q in enumerate(vertex):
                op("RX" if j in pivots else "R", q)
            for row, pivot in zip(rows, pivots):
                for j in supp(row ^ (1 << pivot)):
                    op("CX", vertex[pivot], vertex[j])
        else:
            for e in edges:
                op("RX", e)
            for q in vertex:
                op("R", q)
        parities = []
        vertex_records = []
        for j, q in enumerate(vertex):
            for site in supp(ports[j]):
                op("CX", site, q)
            if method == "cone":
                for i, row in enumerate(auxiliary):
                    if row >> j & 1:
                        op("CX", edges[i], q)
            op("M", q)
            vertex_records.append(records)
            records += 1
        for mask, _sign in logicals:
            parities.append([vertex_records[j] for j in supp(mask)])
        if prior is None:
            first = parities
        else:
            detectors.extend(a + b for a, b in zip(prior, parities))
        prior = parities
        if method == "cone":
            start = records
            for e in edges:
                op("MX", e)
                records += 1
            transpose = [
                sum(((a >> j) & 1) << i for i, a in enumerate(auxiliary))
                for j in range(v)
            ]
            for relation in nullspace(transpose, len(edges)):
                detectors.append([start + i for i in supp(relation)])
            section = independent_section(auxiliary, v)
            for i in range(len(edges)):
                mask = sum((c >> i & 1) << j for j, c in enumerate(section))
                for q in supp(combine(ports, mask)):
                    op("Z_IF", q, start + i)
    rotate(True)
    return dict(
        n=n,
        checks=list(checks),
        targets=targets,
        ticks=ticks,
        readouts=[dict(records=r, sign=s) for r, (_, s) in zip(first, logicals)],
        detectors=detectors,
        terminal_syndrome=True,
    )


def graph(n, checks, word, edges, method="cone", rounds=1):
    """Connected gauging graph with one vertex for each nonidentity factor."""
    target = pauli(word, n)
    sites = supp(target.x | target.z)
    size = len(sites)
    rows = []
    for a, b in edges:
        if (
            type(a) != int
            or type(b) != int
            or not 0 <= a < size
            or not 0 <= b < size
            or a == b
        ):
            raise ValueError("invalid graph edge")
        rows.append((1 << a) | (1 << b))
    if not size or Echelon(rows).rank != size - 1 or len(set(rows)) != len(rows):
        raise ValueError("expected a connected simple graph on the logical support")
    axes = word[1:].replace("I", "Z")
    return construct(
        n,
        checks,
        axes,
        [1 << q for q in sites],
        rows,
        [((1 << size) - 1, target.s)],
        method,
        rounds,
    )
