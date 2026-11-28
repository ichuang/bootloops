"""Four-piece round-robin CCZ (Yoder, Takagi & Chuang, Fig. 5 / Eq. 34).

Physical sites are block-major. At a cut, repeat full constant-check rounds
until consecutive syndromes agree. Keep each cut's X correction for terminal
phase-error decoding. No classical bit in this procedure sees the fault label.
"""

from collections import Counter
from functools import lru_cache
from itertools import product
from ftqc_codes.codes import steane
from ftqc_codes.gf2 import dot, supp, Echelon
from ftqc_decoders.staged import Recovery, syndrome
from ftqc_faults.staged import cat_faults, conjugate_ccz, z_expansion
from .staged import Counterexample, endpoint
from .staged_switch import span


def example():
    code = steane()
    line = next(
        v for v in range(128) if v.bit_count() == 3 and not syndrome(v, code.hx)
    )
    active = supp(line)
    pieces = [[] for _ in range(4)]
    for abc in product(range(3), repeat=3):
        special = [i for i, v in enumerate(abc) if v == 2]
        if len(special) in (0, 3):
            group = 0
        elif len(special) == 1:
            group = 1 + special[0]
        else:
            group = 1 + next(i for i in range(3) if i not in special)
        pieces[group].append([7 * i + active[v] for i, v in enumerate(abc)])
    return dict(
        kind="staged-gadget",
        operation="pieceable-ccz",
        blocks=["A", "B", "C"],
        pieces=pieces,
        cuts=["repeat-Z-until-agreement"] * 3,
    )


def check(spec, rule=None):
    if (
        (set(spec) - {"record_history"})
        != {"kind", "operation", "blocks", "pieces", "cuts"}
        or len(spec["blocks"]) != 3
        or len(spec["pieces"]) != 4
        or spec["cuts"] != ["repeat-Z-until-agreement"] * 3
    ):
        raise ValueError("expected four CCZ pieces and three constant-check cuts")
    keep_history = spec.get("record_history", True)
    if type(keep_history) is not bool:
        raise ValueError("record_history must be Boolean")
    pieces = spec["pieces"]
    if any(
        not isinstance(piece, list) or not 1 <= len(piece) <= 27 for piece in pieces
    ):
        raise ValueError("expected 1..27 gates per piece")
    for piece in pieces:
        for gate in piece:
            if (
                not isinstance(gate, list)
                or len(gate) != 3
                or any(type(q) is not int or q // 7 != b for b, q in enumerate(gate))
            ):
                raise ValueError("CCZ requires one ordered physical site per block")
    small = steane()
    # All amplitudes in a logical basis column must acquire the same sign.
    words = span(small.hx)
    phase_masks = [sum(1 << q for q in gate) for piece in pieces for gate in piece]
    for a, b, c in product(words + [w ^ 127 for w in words], repeat=3):
        word = a | b << 7 | c << 14
        logical = [v.bit_count() % 2 for v in (a, b, c)]
        phase = sum(word & mask == mask for mask in phase_masks) % 2
        if phase != all(logical):
            raise Counterexample(
                "logical CCZ phase", word=word, logical=logical, sign=phase
            )
    ports = endpoint("steane", spec["blocks"])
    return dict(
        status="proved",
        inputs=ports,
        outputs=ports,
        action="CCZ",
        ideal=dict(
            coset_words=(2 * len(words)) ** 3, piece_lengths=list(map(len, pieces))
        ),
        faults=faults(pieces, small, keep_history, rule),
    )


def events(pieces, checks):
    """Only a first fault can occur on the fault-free execution prefix."""
    for q in range(21):
        for x, z in ((1, 0), (0, 1), (1, 1)):
            yield (), ("incoming", q, x, z), x << q, z << q, 0
    for i, piece in enumerate(pieces):
        for j, gate in enumerate(piece):
            tick = (i, j)
            for x, z in product(range(8), repeat=2):
                if x or z:
                    mx, mz = [
                        sum(((mask >> k) & 1) << q for k, q in enumerate(gate))
                        for mask in (x, z)
                    ]
                    yield tick, ("CCZ", x, z), mx, mz, 0
            for q in range(21):
                if q not in gate:
                    for x, z in ((1, 0), (0, 1), (1, 1)):
                        yield tick, ("idle", q, x, z), x << q, z << q, 0
        if i < len(pieces) - 1:
            for r, j in product(range(2), range(len(checks))):
                for mode, x, z, flip in cat_faults(21, checks[j], "Z"):
                    yield (i, r, j), mode, x, z, flip


def faults(pieces, small, keep_history, rule):
    checks = [h << (7 * b) for b in range(3) for h in small.hz]
    xchecks = [h << (7 * b) for b in range(3) for h in small.hx]
    decoder = Recovery(xchecks, checks, rule)
    z_group = Echelon(checks)
    correction = {syndrome(1 << q, small.hz): 1 << q for q in range(7)}
    correction[0] = 0

    @lru_cache(maxsize=None)
    def replay(at, fx, fz, flip):
        x, z = (fx, fz) if at == () else (0, 0)
        edges = set()
        history = []
        for i, piece in enumerate(pieces):
            for j, gate in enumerate(piece):
                if x:
                    z = conjugate_ccz(x, z, edges, gate)
                if at == (i, j):
                    x ^= fx
                    z ^= fz
            if i == len(pieces) - 1:
                break
            previous = None
            for r in range(4):
                reading = 0
                for j, h in enumerate(checks):
                    bit = dot(x, h)
                    if at == (i, r, j):
                        x ^= fx
                        z ^= fz
                        bit ^= flip
                    reading |= bit << j
                if reading == previous:
                    break
                previous = reading
            else:
                raise Counterexample("constant-check rounds failed to settle", tick=at)
            history.append(reading)
            for b in range(3):
                x ^= correction[(reading >> (3 * b)) & 7] << (7 * b)
        # Distinct syndrome projections are classical branches; within one
        # branch add amplitudes for Pauli terms identical on the output code.
        terms = Counter()
        for mask, amplitude in z_expansion(z, edges).items():
            terms[z_group.reduce(mask)] += amplitude
        return tuple(history), x, tuple(v for v, k in terms.items() if k), len(edges)

    alternatives = Counter()
    sites = set()
    maximum = 0
    for at, mode, x, z, flip in [(None, ("clean",), 0, 0, 0), *events(pieces, checks)]:
        history, out_x, terms, nedges = replay(at, x, z, flip)
        maximum = max(maximum, nedges)
        for out_z in terms:
            witness = decoder.observe(
                history if keep_history else (), out_x, out_z, dict(tick=at, mode=mode)
            )
            if witness:
                raise Counterexample("terminal recovery collision", **witness)
        alternatives[mode[0]] += 1
        if mode[0] not in ("clean", "incoming"):
            sites.add((at, mode[0], None if mode[0] == "CCZ" else mode[1]))
    replay.cache_clear()
    return dict(
        order=1,
        alternatives=dict(alternatives),
        noisy_locations=len(sites),
        maximum_residual_cz_terms=maximum,
        recovery=decoder.rows(),
        terminal="ideal full syndrome and recovery",
        ideal_resources=[
            "verified cat seeds",
            "classical processing and intermediate X correction",
        ],
        schedule="native serial CCZ; other data idle; parallel cat couplings per check",
    )
