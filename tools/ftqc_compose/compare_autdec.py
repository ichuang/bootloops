"""Compare native and AutDEC decoders on identical weight-zero and one Steane faults.

Run this optional reference in an AutDEC environment with PYTHONPATH unset.
"""
from pathlib import Path
import sys
from unittest.mock import patch


def main():
    import numpy as np
    from autdec import code_capacity
    sys.path.append(str(Path(__file__).resolve().parents[1]))
    from ftqc_codes import basis, codes, interchange
    from ftqc_decoders.dem import DecodingProblem
    from ftqc_decoders.effective_distance import run_decoder

    c = codes.load("steane")
    paired = basis.derive_paired_basis(c)
    # A non-involutive column map in the Steane Hamming check coordinates.
    cols = [sum(((row >> q) & 1) << i for i, row in enumerate(c.hx)) for q in range(c.n)]
    perm = [cols.index((v >> 1) | ((v & 1) << 2)) for v in cols]
    cfg = dict(algorithm="automorphism_ensemble", permutations=[perm])
    def dense(rows):
        return np.array([[(v >> q) & 1 for q in range(c.n)] for v in rows])
    problems = [DecodingProblem.from_css_code(c, sector=axis) for axis in ("Z", "X")]
    faults = [(0, 0)] + [(z << q, x << q) for q in range(c.n) for z, x in ((0,1),(1,0),(1,1))]
    for bits in faults:
        for problem, error in zip(problems, bits):
            result = run_decoder(problem, problem.syndrome_of(error), cfg)
            if not result.converged or problem.observable_of(error ^ result.error_mask):
                raise ValueError("native ensemble failed a single-qubit control")
    lx = [basis.parse_pauli(s["x"], c.n)[0] for s in paired["slots"]]
    lz = [basis.parse_pauli(s["z"], c.n)[1] for s in paired["slots"]]
    autos = [np.eye(c.n, dtype=int), dense([interchange.transport_word(1 << q, perm) for q in range(c.n)])]
    stream = (dense(bits).reshape(-1) for bits in faults)
    with patch.object(code_capacity, "depolarizing_channel", side_effect=stream):
        counts = code_capacity.AutDEC_code_capacity(dense(c.hx), dense(c.hz), dense(lx), dense(lz),
                0.001, len(faults), "BPOSD", {"max_iter": 10, "osd_order": 0}, autos, autos)
    if counts[0] != len(faults):
        raise ValueError("AutDEC did not run the full common fault set")
    print(f"{len(faults)} common patterns; native ensemble failures: 0; "
          f"AutDEC base failures: {counts[1]}; AutDEC ensemble failures: {counts[2]}")


if __name__ == "__main__":
    main()

