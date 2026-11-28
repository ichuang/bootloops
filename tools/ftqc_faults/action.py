"""Single-location residual spread for disjoint physical Clifford ticks."""
from itertools import product

from ftqc_codes import automorphisms as aut


def single_fault_spread(width, n_data, ticks):
    """Enumerate post-gate Paulis and XYZ idles. Spectators can propagate into data."""
    data = (1 << n_data) - 1
    xs = [1 << q for q in range(width)]
    zs = [1 << (width + q) for q in range(width)]
    histogram = {}
    for tick in reversed(ticks):
        used = {q for op in tick for q in op[1:]}
        locations = [op[1:] for op in tick] + [[q] for q in range(width) if q not in used]
        for loc in locations:
            for errors in product(range(4), repeat=len(loc)):
                if not any(errors):
                    continue
                mask = 0
                for q, e in zip(loc, errors):
                    if e & 1:
                        mask ^= xs[q]
                    if e & 2:
                        mask ^= zs[q]
                weight = ((mask | mask >> width) & data).bit_count()
                histogram[weight] = histogram.get(weight, 0) + 1
        # Prepend the tick to the suffix tableau, ignoring Pauli signs for weight.
        for op in tick:
            gate, q, *other = op
            if gate == "H":
                xs[q], zs[q] = zs[q], xs[q]
            elif gate in ("S", "SDG"):
                xs[q] ^= zs[q]
            elif other:
                t = other[0]
                if gate == "CX":
                    xs[q] ^= xs[t]
                    zs[t] ^= zs[q]
                elif gate == "CZ":
                    xs[q] ^= zs[t]
                    xs[t] ^= zs[q]
                elif gate == "SWAP":
                    xs[q], xs[t] = xs[t], xs[q]
                    zs[q], zs[t] = zs[t], zs[q]
    return {"post_gate_and_idle_faults": sum(histogram.values()),
            "data_weight_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "max_data_weight": max(histogram, default=0)}


def check_data_action(width, n, ticks, steps):
    """Prove the same signed data Pauli images, for arbitrary entangled inputs."""
    layers = [{"layer": tick} for tick in ticks]
    aut.validate(width, layers)
    for q in range(n):
        for p in ((1 << q, 0, 0), (0, 1 << q, 0)):
            if aut.conjugate(p, layers) != aut.conjugate(p, steps):
                raise ValueError(f"installed circuit changes signed data image at site {q}")

