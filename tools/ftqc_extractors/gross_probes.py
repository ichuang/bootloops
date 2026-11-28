"""Complementary encoded CSS experiments for the gross Z instrument."""

from functools import reduce
from operator import xor

from ftqc_codes.gf2 import Echelon
from ftqc_faults.audit.circuit import Circuit
from .algebra import bits, echelon, kernel, reduce_row
from .gross import gross_z, port_frame, shifted_z
from .gross_cycle import cycle, X0, Z0, NWIRES


def quotient(rows, stabilizers):
    span = echelon(stabilizers)
    out = []
    for row in rows:
        rem = reduce_row(row, span)
        if rem:
            span[rem.bit_length() - 1] = rem
            out.append(row)
    return out


def encoder(memory, basis):
    # Uniform superposition over the X-check group gives encoded all-zero.
    # Its CSS dual gives encoded all-plus. All encoder locations are ideal.
    group = Echelon(memory.hx if basis == "Z" else memory.hz)
    pivots = set(group.pivots)
    layers = [
        [
            ("!RX" if ((q in pivots) == (basis == "Z")) else "!R", q)
            for q in range(memory.n)
        ]
    ]
    for p, h in zip(group.pivots, group.basis):
        for q in bits(h ^ (1 << p)):
            layers.append([("!CX", p, q) if basis == "Z" else ("!CX", q, p)])
    return layers


def experiment(rounds, basis):
    if type(rounds) is not int or not 1 <= rounds <= 32 or basis not in ("X", "Z"):
        raise ValueError("choose 1 to 32 syndrome rounds and an X or Z readout basis")
    mem, (_, hx, hz), ports, labels = gross_z()
    layers = encoder(mem, basis)
    idle = [("I", q) for q in range(mem.n)]
    edge = list(range(mem.n, X0))
    layers.append([("RX", q) for q in edge] + idle)
    body = cycle()
    for _ in range(rounds):
        layers.extend(body)
    layers.extend(
        [
            [("MX", q) for q in edge] + idle,
            [("!M" if basis == "Z" else "!MX", q) for q in range(mem.n)],
        ]
    )
    width = NWIRES - X0
    split = width * rounds
    end = split + len(edge)
    is_z = basis == "Z"
    first = Z0 - X0 if is_z else 0
    size = len(hz) if is_z else len(hx)
    # First port outcomes carry quantum information. Old Z checks and all
    # extended X checks (including gauge 6) start at known + signs.
    known = len(mem.hz) if is_z else size
    detectors = [[first + i] for i in range(known)]
    for r in range(1, rounds):
        detectors.extend(
            [width * (r - 1) + first + j, width * r + first + j] for j in range(size)
        )
    paths, cycles = port_frame(ports, mem.n, len(edge))
    sites = list(bits(shifted_z()))

    def terminal(row):
        rec = [end + q for q in bits(row)]
        if not is_z:
            path = reduce(xor, (p for q, p in zip(sites, paths) if row >> q & 1), 0)
            rec.extend(split + i for i in bits(path))
        return sorted(rec)

    if not is_z:
        detectors.extend([split + i for i in bits(c)] for c in cycles)
    detectors.extend(terminal(h) for h in (mem.hz if is_z else mem.hx))
    if is_z:
        reps = quotient([shifted_z(), *kernel(mem.hx, mem.n)], mem.hz)
        obs = [list(range(first + known, first + size))]
    else:
        reps = quotient(kernel([*mem.hz, shifted_z()], mem.n), mem.hx)
        obs = []
    obs.extend(terminal(h) for h in reps)
    if len(reps) != (12 if is_z else 11):
        raise ValueError("incorrect logical quotient")
    return Circuit(
        layers,
        detectors=detectors,
        observables=obs,
        outputs=[],
        name=f"BB144 unprimed Z, {rounds} rounds, {basis} probe",
    )
