"""Independent full-pulse Clifford probes of the adaptive factory."""

from .bicycle_install import install
from ftqc_extractors.bicycle import layout


def physical_probe(study, seed=1):
    import stim

    p = install(study, True)
    m, xs, zs = layout()
    ref = p["width"]

    def pauli(x, z):
        g = stim.PauliString(ref + 1)
        for i in range(ref + 1):
            g[i] = (0, 1, 3, 2)[((x >> i) & 1) + 2 * ((z >> i) & 1)]
        return g

    checks = [pauli(x, 0) for x in m.hx] + [pauli(0, z) for z in m.hz]
    checks.extend(pauli(0, z ^ (1 << ref) if i == 7 else z) for i, z in enumerate(zs))
    checks.append(pauli(xs[7] | (1 << ref), 0))
    checks.extend(pauli(0, 1 << i) for i in range(144, ref))
    sim = stim.TableauSimulator(seed=seed)
    sim.set_inverse_tableau(
        stim.Tableau.from_stabilizers(checks, allow_redundant=True).inverse()
    )
    records = {}
    for step in p["steps"]:
        idx = sum(records[k] << i for i, k in enumerate(step["controls"]))
        g = p["templates"][step["arms"][idx]]
        circ = stim.Circuit()
        for tick in g["ticks"]:
            for op, *qs in tick:
                if op == "T":
                    circ.append("RX", qs)
                    circ.append("S", qs)
                elif op == "Z_IF":
                    circ.append(
                        "CZ", [stim.target_rec(qs[1] - circ.num_measurements), qs[0]]
                    )
                else:
                    circ.append({"SDG": "S_DAG"}.get(op, op), qs)
        sim.do(circ)
        nr = circ.num_measurements
        stream = sim.current_measurement_record()[-nr:] if nr else []
        if any(sum(stream[j] for j in ds) % 2 for ds in g["detectors"]):
            raise ValueError("clean Clifford probe triggered Gross detector")
        if step["record"]:
            (read,) = g["readouts"]
            records[step["record"]] = (
                read["sign"] + sum(stream[j] for j in read["records"])
            ) % 2
    # Independently contract five logical qubits with |+i> inputs and
    # the observed physical transcript. Each parity branch has unit phase.
    amps = [1 + 0j] * 32
    for i in range(15):
        a, b, c, d = [records[f"r{i}.{bit}"] for bit in "abcd"]
        for state in range(32):
            parity = (state & (1 | ((i + 1) << 1))).bit_count() % 2
            anc = parity ^ c
            raw = anc ^ a
            amp = (1j**raw) * (-1) ** (raw * b + anc * d)
            if a ^ c:
                amp *= (-1j) ** anc
            if parity and (b ^ d ^ ((a ^ c) & (1 ^ c))):
                amp = -amp
            amps[state] *= amp
    out = [0j, 0j]
    syndrome = sum(records[f"check.{j}"] << j for j in range(1, 5))
    for state, amp in enumerate(amps):
        out[state & 1] += (
            amp * (1j if state & 1 else 1) * (-1) ** ((state & syndrome).bit_count())
        )
    u, v = out
    norm = abs(u) ** 2 + abs(v) ** 2
    if not norm:
        raise ValueError("impossible postselection record")
    physical_y = pauli(xs[0], zs[0])
    physical_y.sign = (-1) ** (((xs[0] & zs[0]).bit_count() - 1) // 2)
    measured = [
        sim.peek_observable_expectation(g)
        for g in (pauli(xs[0], 0), physical_y, pauli(0, zs[0]))
    ]
    predicted = [
        2 * (u.conjugate() * v).real / norm,
        2 * (u.conjugate() * v).imag / norm,
        (abs(u) ** 2 - abs(v) ** 2) / norm,
    ]
    if any(abs(a - b) > 1e-8 for a, b in zip(measured, predicted, strict=True)):
        raise ValueError(
            f"pulse probe {measured} differs from logical contraction {predicted}"
        )
    if any(
        sim.peek_observable_expectation(g) != 1
        for g in (pauli(xs[7] | (1 << ref), 0), pauli(0, zs[7] | (1 << ref)))
    ):
        raise ValueError("spectator reference disturbed")
    return dict(
        seed=seed,
        input="all +i Clifford control",
        output_XYZ=measured,
        final_check_bits=syndrome >> 1,
        spectator_Bell=True,
    )
