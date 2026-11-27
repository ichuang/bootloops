"""Finite seed translations, product connections, and movement-layer checks."""

from fractions import Fraction
from itertools import combinations
from ftqc_codes.gf2 import Echelon, supp
from ftqc_gadgets.surgery_geometry import restricted_cone, transpose
from ftqc_gadgets.surgery_merge import word
from .surgery_products import shift_word, lifted_product


def translated_settings(n, hx, hz, lift, seed, shifts, rounds):
    """Verify the code automorphism and independence of every requested orbit member."""
    x, z = Echelon(hx), Echelon(hz)
    if any(
        shift_word(r, n, lift) not in space
        for rows, space in ((hx, x), (hz, z))
        for r in rows
    ):
        raise ValueError("cyclic shift does not preserve code")
    seen = Echelon(hz)
    result = []
    for shift in shifts:
        mask = shift_word(seed, n, lift, shift)
        if any((r & mask).bit_count() % 2 for r in hx) or not seen.add(mask):
            raise ValueError("translated seed is dependent or not logical")
        # Keep port order equivariant even when the orbit wraps around.
        sites = [q // lift * lift + (q % lift + shift) % lift for q in supp(seed)]
        result.append(
            dict(
                method="gauging",
                n=n,
                hx=hx,
                hz=hz,
                axis="Z",
                sites=sites,
                edges=[[i, i + 1] for i in range(len(sites) - 1)],
                targets=[word(n, mask, "Z")],
                rounds=rounds,
            )
        )
    return result


def product_connection(a, b, other_a, other_b, lift, rounds):
    """Gu et al. Lemma V.3. Expand f1=H'Z transpose and both chain squares."""
    n, hx, hz = lifted_product(a, b, lift)
    m, fx, fz = lifted_product(other_a, other_b, lift)
    if (n, len(hx), len(hz)) != (m, len(fx), len(fz)):
        raise ValueError("incompatible product patches")
    from ftqc_codes.gf2 import dot

    if any(dot(x, zp) != dot(xp, z) for x, xp in zip(hx, fx) for z, zp in zip(hz, fz)):
        raise ValueError("product connection square fails")
    aux = transpose(hz, n)
    cone = restricted_cone(n, hx, hz, fz, aux)
    if not cone["measured"]:
        raise ValueError("connection induces no logical measurement")
    return dict(
        n=n,
        hx=hx,
        hz=hz,
        method="linear",
        ports=fz,
        aux=aux,
        axis="Z",
        targets=[word(n, z, "Z") for z in cone["measured"]],
        rounds=rounds,
    )


def transport(initial, layers, ticks, gate_range):
    """Exact synchronous straight-line trap motion, followed by one circuit tick.

    Coordinates, movement and gate durations are integer device units. All
    occupied traps participate in collision checks, even while waiting.
    """
    if len(layers) != len(ticks) or type(gate_range) != int or gate_range < 0:
        raise ValueError("expected one timed layout per gate tick")

    def point(p):
        if len(p) != 2 or any(type(v) != int for v in p):
            raise ValueError("invalid trap position")
        return tuple(p)

    pos = [point(p) for p in initial]
    if len(set(pos)) != len(pos):
        raise ValueError("initial trap collision")
    elapsed = movement = 0
    for layer, row in zip(layers, ticks):
        dst = pos.copy()
        moved = set()
        for q, p in layer["moves"]:
            if type(q) != int or not 0 <= q < len(pos) or q in moved:
                raise ValueError("invalid moving trap")
            moved.add(q)
            dst[q] = point(p)
        dt, gate_time = layer["move_time"], layer["gate_time"]
        if (
            any(type(v) != int or v < 0 for v in (dt, gate_time))
            or (dt == 0 and pos != dst)
            or (row and not gate_time)
        ):
            raise ValueError("invalid layer duration")
        for i, j in combinations(range(len(pos)), 2):
            if i not in moved and j not in moved:
                continue
            roots = set()
            for k in range(2):
                delta = pos[i][k] - pos[j][k]
                speed = dst[i][k] - dst[j][k] - delta
                if speed:
                    roots.add(Fraction(-delta, speed))
                elif delta:
                    break
            else:
                if len(roots) == 1 and 0 <= next(iter(roots)) <= 1:
                    raise ValueError("collision during transport")
        for name, *args in row:
            wires = args[:1] if name.lstrip("!") == "Z_IF" else args
            if any(q < 0 or q >= len(pos) for q in wires):
                raise ValueError("missing trap")
            if name.lstrip("!") == "CX":
                a, b = (dst[q] for q in wires)
                if sum((i - j) ** 2 for i, j in zip(a, b)) > gate_range**2:
                    raise ValueError("gate outside interaction radius")
        pos = dst
        elapsed += dt + gate_time
        movement += dt
    return dict(elapsed=elapsed, movement=movement, final=pos)
