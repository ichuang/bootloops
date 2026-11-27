"""Statistical summaries conditional on a shared batch acceptance decision."""

import math


def interval(count, total):
    if total == 0:
        return None
    z = 1.959963984540054
    p = count / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total**2)) / denom
    return dict(estimate=p, wilson95=[max(0, center - half), min(1, center + half)])


def summarize(sample, k, physical_qubits, depth):
    n, a = sample["shots"], sample["accepted"]
    h = sample["accepted_mask_counts"]
    if (
        type(n) is not int
        or n < 1
        or type(a) is not int
        or not 0 <= a <= n
        or sample["outputs"] != k
        or not isinstance(h, list)
        or len(h) != 1 << k
        or any(type(c) is not int or c < 0 for c in h)
        or sum(h) != a
    ):
        raise ValueError("inconsistent shot denominators or mask histogram")
    failures = a - h[0]
    correlated = sum(v for m, v in enumerate(h) if m.bit_count() > 1)
    return dict(
        evidence="empirical closed logical Y; individual Wilson 95% intervals",
        shots=n,
        accepted=a,
        discard=interval(n - a, n),
        any_error=interval(failures, a),
        correlated_error=interval(correlated, a),
        correlated_share=interval(correlated, failures),
        batch_error_divided_by_outputs=failures / (k * a) if a else None,
        marginals=[
            interval(sum(v for m, v in enumerate(h) if m >> i & 1), a) for i in range(k)
        ],
        accepted_states_per_attempt=k * a / n,
        physical_footprint=physical_qubits,
        two_qubit_layers=depth,
        qubit_layers_per_accepted_state=(
            physical_qubits * depth * n / (k * a) if a else None
        ),
    )
