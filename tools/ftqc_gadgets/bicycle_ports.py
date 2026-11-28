"""Physical parity meters and Tanner-edge SWAP layers for the BB factory."""

from collections import deque
from functools import lru_cache
from ftqc_codes.gf2 import supp
from ftqc_extractors.gross import cell, gross_z
from ftqc_extractors.block_map import word
from .surgery import Instrument


def meter(x, z, anc, sign=0):
    """Measure a signed tensor Pauli with a single parity meter; restore local bases."""
    sites = supp(x | z)
    if not sites or anc in sites:
        raise ValueError("invalid meter")
    ys = supp(x & z)
    hs = supp(x)
    ticks = [[["R", anc]]]
    if ys:
        ticks.append([["SDG", q] for q in ys])
    if hs:
        ticks.append([["H", q] for q in hs])
    ticks.extend([["CX", q, anc]] for q in sites)
    ticks.append([["M", anc]])
    if hs:
        ticks.append([["H", q] for q in hs])
    if ys:
        ticks.append([["S", q] for q in ys])
    return dict(ticks=ticks, readouts=[dict(records=[0], sign=sign)], detectors=[])


@lru_cache(maxsize=128)
def verify_meter(n, checks, x, z, sign):
    gadget = meter(x, z, n, sign)
    Instrument(
        dict(
            n=n,
            checks=list(checks),
            targets=[word(n, x, z, sign)],
            terminal_syndrome=False,
            **gadget,
        )
    )
    return True


def translated(q, delta):
    return (q // 72) * 72 + cell(q % 72 // 6 + delta[0], q % 6 + delta[1])


@lru_cache(maxsize=1)
def generators():
    m = gross_z()[0]
    generators = {}
    for lrows, rrows, lbase, rbase in [(m.hx, m.hz, 144, 216), (m.hz, m.hx, 216, 144)]:

        def offsets(rows, half):
            return [q % 72 for q in supp(rows[0]) if q // 72 == half]

        for first in offsets(lrows, 0):
            for second in offsets(lrows, 0):
                if first == second:
                    continue
                delta = ((second // 6 - first // 6) % 12, (second % 6 - first % 6) % 6)
                if delta in generators:
                    continue
                right = next(
                    (
                        q
                        for q in offsets(rrows, 1)
                        if translated(q, delta) in offsets(rrows, 1)
                    ),
                    None,
                )
                if right is None:
                    continue
                legs = [[], []]
                for half, rows, base, offset in (
                    (0, lrows, lbase, first),
                    (1, rrows, rbase, right),
                ):
                    for i, row in enumerate(rows):
                        q = half * 72 + cell(i // 6 + offset // 6, i % 6 + offset % 6)
                        t = translated(q, delta)
                        if not ((row >> q) & (row >> t) & 1):
                            raise ValueError("SWAP outside memory Tanner graph")
                        legs[0].append((q, base + i))
                        legs[1].append((base + i, t))
                ticks = []
                for pairs in legs:
                    if len({q for pair in pairs for q in pair}) != 288:
                        raise ValueError("SWAP collision")
                    for reverse in (False, True, False):
                        ticks.append(
                            [["CX", *(reversed(p) if reverse else p)] for p in pairs]
                        )
                # Independently propagate each original data label through both SWAPs.
                labels = list(range(288))
                for pairs in legs:
                    for a, b in pairs:
                        labels[a], labels[b] = labels[b], labels[a]
                if any(labels[translated(q, delta)] != q for q in range(144)):
                    raise ValueError("shift permutation mismatch")
                # Propagate all data X and Z generators through the emitted CXs.
                # A bit in row i identifies an input Pauli with output on site i.
                initial = [1 << q for q in range(144)] + [0] * 144
                xs, zs = initial.copy(), initial.copy()
                for tick in ticks:
                    for _, control, target in tick:
                        xs[target] ^= xs[control]
                        zs[control] ^= zs[target]
                expected = [0] * 288
                for q in range(144):
                    expected[translated(q, delta)] = 1 << q
                if xs != expected or zs != expected:
                    raise ValueError("emitted CX route changes the signed data action")
                generators[delta] = ticks
    return generators


def route(dx, dy):
    paths = {(0, 0): []}
    queue = deque(paths)
    while queue:
        at = queue.popleft()
        if at == (dx % 12, dy % 6):
            return paths[at]
        for d in generators():
            to = ((at[0] + d[0]) % 12, (at[1] + d[1]) % 6)
            if to not in paths:
                paths[to] = paths[at] + [d]
                queue.append(to)
    raise ValueError("unroutable translation")
