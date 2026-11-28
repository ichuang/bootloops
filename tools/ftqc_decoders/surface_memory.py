"""Decode the measured Pauli sector of a closed rotated-patch memory experiment."""

import argparse
import json

from ftqc_faults.audit.library import surface
from ftqc_faults.audit.record_model import RecordModel
from .dem import DecodingProblem, transpose_bitmasks
from .effective_distance import run_decoder


def benchmark(distance, rounds, basis, bad=False):
    basis = basis.upper()
    if basis not in ("X", "Z") or rounds < 1:
        raise ValueError("expected positive syndrome rounds and X or Z readout")
    circ = surface(distance, rounds, basis, bad)
    records = RecordModel(circ)
    selected = [
        i
        for i, det in enumerate(circ.detectors)
        if all(circ.meas[j][3] == basis for j in det)
    ]

    def project(det):
        return sum(((det >> i) & 1) << j for j, i in enumerate(selected))

    modes = []
    for site, group in enumerate(records.groups):
        for label in group:
            det, logical = records.replay([(site, label)])
            modes.append((project(det), logical))
    signatures = sorted(set(modes) - {(0, 0)})
    rows = transpose_bitmasks([d for d, _ in signatures], len(selected))
    logical_row = sum(bit << i for i, (_, bit) in enumerate(signatures))
    model = DecodingProblem.from_matrices(rows, [logical_row], len(signatures))
    if not model.is_graphlike():
        raise ValueError("sector projection contains a non-graphlike signature")
    predictions = {}
    for syndrome, _ in signatures + [(0, 0)]:
        if syndrome not in predictions:
            match = run_decoder(model, syndrome, "mwpm")
            if not match.converged or model.syndrome_of(match.error_mask) != syndrome:
                raise ValueError("incomplete matching")
            predictions[syndrome] = match.observable_mask
    failures = sum(predictions[d] != bit for d, bit in modes)
    return dict(
        distance=distance,
        rounds=rounds,
        basis=basis,
        hook_order="bad" if bad else "standard",
        circuit_sha256=circ.sha256(),
        data=distance**2,
        physical_qubits=circ.n,
        ticks=len(circ.ticks),
        locations=len(records.groups),
        single_modes=len(modes),
        failed_single_modes=failures,
        detector_sector=basis,
        graph_edges=len(signatures),
        decoder="equal-weight sector MWPM",
        boundary="noisy product preparation and destructive readout",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=int, default=3)
    parser.add_argument("--rounds", type=int)
    parser.add_argument("--basis", choices=["X", "Z"], default="Z")
    parser.add_argument("--bad", action="store_true")
    args = parser.parse_args()
    rounds = args.distance if args.rounds is None else args.rounds
    print(json.dumps(benchmark(args.distance, rounds, args.basis, args.bad), indent=2))


if __name__ == "__main__":
    main()
