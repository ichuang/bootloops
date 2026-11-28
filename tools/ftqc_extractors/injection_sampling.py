"""Decoder parameters in Liu et al.'s released injection study."""

import hashlib
import importlib.metadata
from collections import Counter
from functools import lru_cache


def build_decoder(circuit):
    from stimbposd.bp_osd import BPOSD

    return BPOSD(
        circuit.detector_error_model(),
        max_bp_iters=1000,
        bp_method="ms",
        osd_order=10,
        osd_method="osd_cs",
        ms_scaling_factor=0,
    )


def sample(path, mask_size, shots, seed):
    import numpy as np
    import stim

    if min(shots, mask_size) < 1:
        raise ValueError("positive shot count and postselection width required")
    data = path.read_bytes()
    circuit = stim.Circuit(data.decode())
    if mask_size > circuit.num_detectors:
        raise ValueError("postselection exceeds detector count")
    decoder = build_decoder(circuit)
    k = circuit.num_observables

    @lru_cache(maxsize=65536)
    def predict(row):
        synd = np.unpackbits(np.frombuffer(row, dtype=np.uint8), bitorder="little")[
            : circuit.num_detectors
        ]
        return sum(int(v) << i for i, v in enumerate(decoder.decode(synd)))

    sampler = circuit.compile_detector_sampler(seed=seed)
    counts = Counter()
    for start in range(0, shots, 4096):
        det, obs = sampler.sample(
            min(shots - start, 4096), separate_observables=True, bit_packed=True
        )
        full, tail = divmod(mask_size, 8)
        keep = ~np.any(det[:, :full], axis=1)
        if tail:
            keep &= (det[:, full] & ((1 << tail) - 1)) == 0
        for row, out in zip(det[keep], obs[keep]):
            counts[
                int.from_bytes(out.tobytes(), "little") ^ predict(row.tobytes())
            ] += 1
    return dict(
        circuit_sha256=hashlib.sha256(data).hexdigest(),
        postselected_detectors=mask_size,
        seed=seed,
        shots=shots,
        accepted=sum(counts.values()),
        outputs=k,
        accepted_mask_counts=[counts[i] for i in range(1 << k)],
        sampler_batch_size=4096,
        versions={
            p: importlib.metadata.version(p)
            for p in ("stim", "stimbposd", "ldpc", "numpy")
        },
        decoder="BP+OSD min-sum; max_iter=1000, osd_cs order 10, ms_scaling_factor=0",
    )


def main():
    import argparse
    import json
    from pathlib import Path
    from ftqc_extractors.injection_circuit import paths, run

    ap = argparse.ArgumentParser(
        description="Sample the pinned logical-Y injection experiment with BP+OSD."
    )
    ap.add_argument("case")
    ap.add_argument("--shots", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("-o", type=Path, required=True)
    a = ap.parse_args()
    _, path, _, _ = paths(a.case)
    report = run({"case": a.case})
    a.o.write_text(
        json.dumps(
            sample(path, report["postselect_width"], a.shots, a.seed),
            sort_keys=True,
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
