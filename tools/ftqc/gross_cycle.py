"""Gate layers for one extraction of the BB144 unprimed-Z merged checks."""

from .algebra import bits
from .gross import cell, gross_z, shifted_z

# The first 162 wires hold data. Meters follow the row order in gross_z().
X0, Z0, PORT = 162, 235, 307
NWIRES = 319
A = ((3, 0), (0, 1), (0, 2))
B = ((0, 3), (1, 0), (2, 0))


def auxiliary_layers(ports, labels):
    # Fix the colored subgraph and both early interactions in Fig. 14(bottom).
    # Complete the coloring on the remaining edges of Fig. 9 deterministically.
    edges = []
    for j, h in enumerate(ports):
        for a in bits(h >> 144):
            edges.append((PORT + j, 144 + a))
    gauge = X0 + 72
    aux = lambda i, j: 144 + labels.index(cell(i, j))
    sites = list(bits(shifted_z()))

    def port(half, i, j):
        return PORT + sites.index(72 * half + cell(i, j))

    colored = {}

    def fix(meter, triples):
        for i, j, color in triples:
            colored[meter, aux(i, j)] = color

    fix(gauge, [(3, 2, 0), (3, 1, 1), (0, 4, 2), (1, 4, 3)])
    edges.extend(e for e in colored)
    fix(port(1, 2, 4), [(0, 4, 1), (1, 4, 2)])
    fix(port(0, 3, 4), [(3, 2, 1), (0, 4, 3)])
    fix(port(1, 3, 4), [(3, 1, 0), (1, 4, 1)])
    fix(port(0, 3, 3), [(3, 2, 2), (3, 1, 3)])
    if not colored.keys() <= set(edges):
        raise ValueError("Figure 14 support mismatch")
    # The original green edge moved to early leaves a green gap at both ends.
    gap = {port(1, 3, 4), aux(3, 1)}

    def options(edge):
        m, q = edge
        used = {v for (a, b), v in colored.items() if a == m or b == q}
        if m in gap or q in gap:
            used.add(2)
        return sorted({1, 2, 3} - used)

    def complete():
        rest = [e for e in edges if e not in colored]
        if not rest:
            return True
        edge = min(rest, key=lambda e: len(options(e)))
        for color in options(edge):
            colored[edge] = color
            if complete():
                return True
            del colored[edge]
        return False

    if not complete():
        raise ValueError("no proper completion of port coloring")
    layers = [[] for _ in range(4)]
    for (m, q), color in sorted(colored.items()):
        layers[color].append(("CX", m, q) if m == gauge else ("CX", q, m))
    return layers


def cycle():
    memory, (n, hx, hz), ports, labels = gross_z()
    meters = list(range(X0, NWIRES))
    layers = [[("RX" if q < Z0 else "R", q) for q in meters]]
    layers.extend(auxiliary_layers(ports, labels))
    layers.append([("CX", q, PORT + j) for j, q in enumerate(bits(shifted_z()))])
    # Fig. 13a in cyclic order II through VIII. The prep/readout barriers and
    # placement of port couplings are explicit compilation choices.
    x_order = [None, (0, A[1]), (1, B[1]), (1, B[0]), (1, B[2]), (0, A[0]), (0, A[2])]
    z_order = [(1, A[0]), (1, A[2]), (0, B[0]), (0, B[1]), (0, B[2]), (1, A[1]), None]
    for j, (xt, zt) in enumerate(zip(x_order, z_order)):
        layer = []
        for lab in range(72):
            i, k = divmod(lab, 6)
            if xt:
                half, (u, v) = xt
                layer.append(("CX", X0 + lab, 72 * half + cell(i + u, k + v)))
            if zt:
                half, (u, v) = zt
                layer.append(("CX", 72 * half + cell(i - u, k - v), Z0 + lab))
        if j == 0:
            layer.extend(("CX", X0 + lab, 144 + i) for i, lab in enumerate(labels))
        layers.append(layer)
    layers.append([("MX" if q < Z0 else "M", q) for q in meters])
    # All data wires are live throughout. Count an idle at each waiting meter
    # between its preparation and readout, even if its gate has been moved.
    for layer in layers:
        used = {q for op in layer for q in op[1:]}
        if len(used) != sum(len(op) - 1 for op in layer):
            raise ValueError("tick collision")
        layer.extend(("I", q) for q in range(NWIRES) if q not in used)
    validate_cycle(layers, hx, hz)
    return layers


def validate_cycle(layers, hx, hz):
    # CSS CNOT transfer in reverse: each measured meter must pull back to its
    # own prepared Pauli and the intended data Pauli. Other same-basis
    # prepared meters are harmless. Opposite-basis factors are forbidden.
    for is_x, rows, start in ((True, hx, X0), (False, hz, Z0)):
        for j, row in enumerate(rows):
            mask = 1 << (start + j)
            for layer in reversed(layers):
                for op in layer:
                    if op[0] != "CX":
                        continue
                    c, t = op[1:]
                    src, dst = (c, t) if is_x else (t, c)
                    if mask >> src & 1:
                        mask ^= 1 << dst
            forbidden = mask >> Z0 if is_x else (mask >> X0) & ((1 << (Z0 - X0)) - 1)
            if mask & ((1 << X0) - 1) != row or forbidden:
                raise ValueError(f"incorrect ideal check at meter {start + j}")
