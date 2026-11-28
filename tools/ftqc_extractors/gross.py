"""The mono-layer unprimed Z measurement of Cross et al., arXiv:2407.18393v4."""

from functools import reduce
from operator import xor

from ftqc_codes.codes import load

from .algebra import bits, css_dimension, echelon, kernel, parity, reduce_row, transpose

# Table 1 and Eqs. (63)-(64). A memory site is L(x^i y^j) or R(x^i y^j).
LEFT_Z = ((2, 3), (2, 5), (3, 0), (3, 3), (3, 4), (3, 5))
RIGHT_Z = ((1, 3), (1, 5), (2, 0), (2, 4), (3, 2), (3, 4))
# Four red X-check-labelled ancilla data qubits incident to gauge check 6, Fig. 9.
GAUGE_6 = ((0, 4), (1, 4), (3, 2), (3, 1))


def cell(i, j):
    return 6 * (i % 12) + (j % 6)


def shifted_z(i=0, j=0):
    return sum(
        1 << (half * 72 + cell(x + i, y + j))
        for half, terms in enumerate((LEFT_Z, RIGHT_Z))
        for x, y in terms
    )


def gross_z():
    """Return the BB memory, merged CSS code and the Z port-to-check labels."""
    memory = load("bb144")
    target = shifted_z()
    incident = [i for i, row in enumerate(memory.hx) if row & target]
    if len(incident) != 18 or any(
        (memory.hx[i] & target).bit_count() != 2 for i in incident
    ):
        raise ValueError("incorrect BB144 coordinates for the unprimed Z support")
    anc = {label: 1 << (memory.n + j) for j, label in enumerate(incident)}
    x = [row | anc.get(i, 0) for i, row in enumerate(memory.hx)]
    gauge = sum(anc[cell(i, j)] for i, j in GAUGE_6)
    x.append(gauge)
    ports = [
        (1 << q) | sum(anc[i] for i in incident if memory.hx[i] >> q & 1)
        for q in bits(target)
    ]
    z = memory.hz + ports
    return memory, (memory.n + len(anc), x, z), ports, incident


def logical_action():
    """Check that projection measures just Z and preserves the surviving Pauli algebra."""
    memory, merged, ports, incident = gross_z()
    n, x, z = merged
    data = (1 << memory.n) - 1
    target = shifted_z()
    old_x, old_z = echelon(memory.hx), echelon(memory.hz)
    if (
        css_dimension((memory.n, memory.hx, memory.hz)) != 12
        or css_dimension(merged) != 11
    ):
        raise ValueError("unexpected number of encoded qubits")
    if (
        reduce(xor, ports, 0) != target
        or any(parity(h, target) for h in memory.hx)
        or not reduce_row(target, old_z)
    ):
        raise ValueError("port outcomes do not measure the selected logical Z")

    # Shorten merged Z stabilizers to data by eliminating their ancilla pivots.
    # Project merged X stabilizers to the original X group. Thus no additional
    # data-only Z logical becomes trivial, and the X equivalences are unchanged.
    merged_z = echelon(z)
    data_z = echelon(row for row in merged_z.values() if row <= data)
    expected_z = echelon([*old_z.values(), target])
    if len(data_z) != len(expected_z) or any(
        reduce_row(row, expected_z) for row in data_z.values()
    ):
        raise ValueError("additional Z projection on the data")
    if echelon(row & data for row in x) != old_x:
        raise ValueError("X stabilizer projection mismatch")
    # X representatives commuting with the measured Z must lift through the
    # port constraints. Project the entire merged X centralizer onto data.
    projected_x = echelon(row & data for row in kernel(z, n))
    surviving_x = kernel([*memory.hz, target], memory.n)
    if len(projected_x) != len(surviving_x) or any(
        reduce_row(row, projected_x) for row in surviving_x
    ):
        raise ValueError("a surviving X logical has no lift")
    if any(parity(row, check) for row in kernel(memory.hx, memory.n) for check in x):
        raise ValueError("a data Z logical is disturbed")
    paths, cycles = port_frame(ports, memory.n, len(incident))
    return dict(
        split_frame_paths=[list(bits(path)) for path in paths],
        split_x_parities=[list(bits(c)) for c in cycles],
        memory_n=memory.n,
        memory_k=12,
        merged_n=n,
        merged_k=11,
        measured_z=list(bits(target)),
        ancilla_x_check_labels=incident,
        gauge_check=6,
        gauge_support=list(bits(x[-1])),
        x_rank=len(echelon(x)),
        z_rank=len(merged_z),
    )


def port_frame(ports, first_ancilla, edge_count):
    """Recover a vertex potential from edge X outcomes along a spanning tree."""
    adjacency = [[] for _ in ports]
    incidence = [row >> first_ancilla for row in ports]
    endpoints = transpose(incidence, edge_count)
    for edge, pair in enumerate(endpoints):
        if pair.bit_count() != 2:
            raise ValueError("port incidence requires two endpoints per ancilla")
        u, v = bits(pair)
        adjacency[u].append((v, edge))
        adjacency[v].append((u, edge))
    paths = {0: 0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v, edge in adjacency[u]:
            if v not in paths:
                paths[v] = paths[u] ^ (1 << edge)
                stack.append(v)
    if len(paths) != len(ports):
        raise ValueError("disconnected port graph")
    paths = [paths[i] for i in range(len(ports))]
    cycles = []
    for edge, pair in enumerate(endpoints):
        u, v = bits(pair)
        cycle = paths[u] ^ paths[v] ^ (1 << edge)
        if cycle:
            cycles.append(cycle)
    # A path's boundary must be the vertex plus the reference vertex.
    # These identities verify the frame on all physically allowed X records.
    for vertex, path in enumerate(paths):
        if reduce(xor, (endpoints[e] for e in bits(path)), 0) != (1 ^ (1 << vertex)):
            raise ValueError("incorrect split correction")
    if len(echelon(cycles)) != edge_count - len(ports) + 1 or any(
        parity(c, row) for c in cycles for row in incidence
    ):
        raise ValueError("incomplete split parity checks")
    return paths, cycles


def measurement():
    """Export ideal joint Pauli measurements and their signed merge/split readout."""
    from ftqc_codes.codes import CSSCode

    memory, (n, x, z), ports, labels = gross_z()
    sites = list(bits(shifted_z()))
    paths, cycles = port_frame(ports, memory.n, len(labels))
    ancillas = list(range(memory.n, n))
    return dict(
        model="ideal-pauli-measurement",
        prepare_plus=ancillas,
        measure_z_ports=[list(bits(h)) for h in ports],
        logical_z_record_parity=list(range(len(ports))),
        measure_x_ancillas=ancillas,
        x_record_constraints=[list(bits(c)) for c in cycles],
        data_z_frame=[
            dict(site=q, x_record_parity=list(bits(p))) for q, p in zip(sites, paths)
        ],
        merged_code=CSSCode("BB144 unprimed Z", n, x, z).to_json(),
    )
