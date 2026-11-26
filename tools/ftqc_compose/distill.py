"""Correcting 15-to-1 magic-state preparation (He et al., Proposition 7.9).

Integer masks use labels 1..15 for the nonzero columns of the punctured
Reed-Muller code. Phase exponents are in units of pi/4.
"""

from itertools import combinations, product

from ftqc_codes.gf2 import dot, nullspace, rank

from .supports import Family


def syndrome(error, hx, hz):
    x, z = error
    return tuple(dot(z, h) for h in hx) + tuple(dot(x, h) for h in hz)


def phase_sum(exponents):
    """Exact sum of eighth roots in the independent basis 1,w,w²,w³."""
    coeff = [0] * 4
    for exponent in exponents:
        k = exponent % 8
        coeff[k % 4] += 1 if k < 4 else -1
    return tuple(coeff)


def teleport_residual(bit, outcome, flip, phase):
    # CX from the data to X^flip Z^phase |T>, read the resource in Z, then
    # apply S to the data for outcome one. Remove the desired data T phase.
    original_bit = bit ^ outcome ^ flip
    return (original_bit + 4 * phase * original_bit + 2 * outcome * bit - bit) % 8


def code_check():
    n = 15
    full = (1 << n) - 1
    hx = [sum(1 << (j - 1) for j in range(1, 16) if j & (1 << i)) for i in range(4)]
    hz = nullspace([*hx, full], n)
    if (rank(hx), rank(hz)) != (4, 10) or any(dot(a, b) for a in hx for b in hz):
        raise ValueError("invalid punctured Reed-Muller CSS checks")
    linear = [0]
    for row in hx:
        linear += [x ^ row for x in linear]
    # In this logical frame the odd coset encodes zero: transversal T has
    # phases (w^7,1), implementing logical T up to a constant phase.
    cosets = [[w ^ full for w in linear], linear]
    if any(
        w.bit_count() % 8 != (7 + b) % 8
        for b, words in enumerate(cosets)
        for w in words
    ):
        raise ValueError("transversal T has the wrong logical phase")
    identity = (0, 0)
    errors = [identity] + [
        (x << i, z << i) for i in range(n) for x, z in ((1, 0), (0, 1), (1, 1))
    ]
    decoder = {syndrome(e, hx, hz): e for e in errors}
    if len(decoder) != 46:
        raise ValueError("single-site minimum-weight correction is ambiguous")
    # The identity and the three Paulis span a one-site error, including
    # off-diagonal pairs. Syndrome measurement removes unequal syndromes.
    pauli_pairs = 0
    for i in range(n):
        local = [identity, (1 << i, 0), (0, 1 << i), (1 << i, 1 << i)]
        for left, right in product(local, repeat=2):
            sl, sr = syndrome(left, hx, hz), syndrome(right, hx, hz)
            if sl == sr and (decoder[sl] != left or decoder[sr] != right):
                raise ValueError("surviving syndrome branch was not corrected")
            pauli_pairs += 1

    # Independently evaluate the decoded one-site injection Kraus operators.
    # After removing ideal transversal T the branch is diagonal at site i,
    # hence lies in span(I,Z_i). Only these two syndromes can occur. Following
    # correction, the diagonal logical entries are coset averages. Check both
    # entries match the same scalar exactly, retaining phases for cross terms.
    checked = 0
    for i, c, flip, phase, correct_z in product(
        range(n), range(2), range(2), range(2), range(2)
    ):
        residual = [teleport_residual(b, c, flip, phase) for b in (0, 1)]
        if flip == phase == 0 and residual != [c, c]:
            raise ValueError("fault-free injection disagrees with T")
        if flip == 0 and phase == 1 and (residual[1] - residual[0]) % 8 != 4:
            raise ValueError("resource Z did not induce a data Z")
        expected = phase_sum([residual[0], residual[1] + 4 * correct_z] * 8)
        for words in cosets:
            exponents = []
            for word in words:
                bit = (word >> i) & 1
                exponents.append(residual[bit] + 4 * correct_z * bit)
            if phase_sum(exponents) != expected:
                raise ValueError("corrected Kraus operator depends on the encoded bit")
        checked += 1
    # Two Z faults on distinct resources always mimic a third single Z.
    # The residual commutes with the checks and anticommutes with logical X.
    pair_count = 0
    witness = None
    for i, j in combinations(range(n), 2):
        error = (1 << i) | (1 << j)
        correction = decoder.get(syndrome((0, error), hx, hz))
        if not correction or correction[0]:
            raise ValueError("missing weight-one Z correction")
        residual = error ^ correction[1]
        if (
            residual.bit_count() != 3
            or any(syndrome((0, residual), hx, hz))
            or not dot(full, residual)
        ):
            raise ValueError("pair witness does not produce logical Z")
        if witness is None:
            witness = {
                "input_labels": [i + 1, j + 1],
                "decoder_Z_label": correction[1].bit_length(),
                "logical": "Z",
            }
        pair_count += 1
    return {
        "code": [15, 1, 3],
        "check_ranks": [4, 10],
        "logical_T_phase_exponents": [7, 0],
        "verified_single_syndromes": len(decoder) - 1,
        "single_site_pauli_pairs": pauli_pairs,
        "decoded_injection_kraus_checks": checked,
        "harmful_input_Z_pairs": pair_count,
        "pair_witness": witness,
    }


def output_check(vertices, inputs, outputs, outer):
    sources = [v for v in vertices if v["ideal"] == "prepare-T"]
    sinks = [v for v in vertices if v["ideal"] == "RM15-inject-correct-decode"]
    if inputs or len(vertices) != 16 or len(sources) != 15 or len(sinks) != 1:
        raise ValueError(
            "distillation requires fifteen native magic preparations and one ideal distiller"
        )
    sink = sinks[0]
    if (
        set(sink["inputs"]) != {v["outputs"][0] for v in sources}
        or outputs != sink["outputs"]
    ):
        raise ValueError(
            "distiller wiring does not match its magic-state sources and output"
        )
    # Arbitrary logical faults at the final operation must also be covered in
    # the outer argument. Its singleton vanishes only at physical substitution.
    ids = [v["id"] for v in sources]
    necessary = Family(outer.sites, [*combinations(ids, 2), [sink["id"]]])
    if not outer.covers(necessary):
        raise ValueError("logical family misses a resource pair or distiller fault")
    return {
        "method": "Proposition 7.9, finite CSS and injection verification",
        "ideal_output": "T|+>",
        "output_metric": "trace distance (half trace norm)",
        "circuit_assumption": "fault-free preparation of encoded plus, local injections, syndrome correction and decoding",
        "code_check": code_check(),
    }
