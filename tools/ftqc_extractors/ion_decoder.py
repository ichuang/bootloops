"""Syndrome differences and a joint BP+OSD decoder on the released H2 circuits."""

from itertools import combinations
from ftqc_codes.gf2 import Echelon, supp


def detector_masks(wires, initial, deterministic, sequential):
    index = [[i for i, q in enumerate(wires) if s == q] for s in range(27)]
    gauges = sorted(i for q in (20, 21, 22) for i in index[q])
    masks = []
    basis = Echelon()

    def add(ids, frames=()):
        for frame in frames:
            trial = sorted(set(ids) ^ set(frame))
            mask = sum(1 << i for i in trial)
            if deterministic(trial):
                if basis.add(mask):
                    masks.append(mask)
                return True
        return False

    for i in initial:
        add([i], [()])
    for q in (*range(5), *range(15, 20), 24, 25, 26):
        if q >= 24 and sequential:
            continue
        for a, b in zip(index[q], index[q][1:]):
            available = [j for j in gauges if a < j < b]
            frames = (
                p for r in range(len(available) + 1) for p in combinations(available, r)
            )
            add([a, b], frames)
    # Final stabilizers of the bicycle code, in the destructive data basis.
    offsets = (1, 2, 4) if sequential else (2, 1, 0)
    a, b, off = offsets
    for i in range(5):
        sites = [5 + i, 5 + (i + a) % 5, 10 + (i + off) % 5, 10 + (i + off + b) % 5]
        data = [index[q][-1] for q in sites]
        # Physical QIR check-meter numbering need not follow polynomial order.
        for meter in (*range(5), *range(15, 20)):
            if add(data + [index[meter][-1]], [()]):
                break
    return masks


def predictor(circuit, detectors, observables):
    import numpy as np
    from scipy.sparse import csc_matrix
    import stim
    from ldpc import BpOsdDecoder
    from ftqc_decoders.dem import DecodingProblem

    for mask in detectors:
        circuit.append(
            "DETECTOR",
            [stim.target_rec(j - circuit.num_measurements) for j in supp(mask)],
        )
    for i, ids in enumerate(observables):
        circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(j - circuit.num_measurements) for j in ids],
            i,
        )
    model = DecodingProblem.from_stim_dem_text(
        str(circuit.detector_error_model()), merge_parallel=False
    )
    rs, cs = [], []
    for col, mask in enumerate(model.mech_detectors):
        for row in supp(mask):
            rs.append(row)
            cs.append(col)
    H = csc_matrix(
        (np.ones(len(rs), dtype=np.uint8), (rs, cs)),
        shape=(len(detectors), len(model.priors)),
    )
    dec = BpOsdDecoder(
        H,
        error_channel=[float(p) for p in model.priors],
        max_iter=20,
        osd_method="OSD_0",
        osd_order=0,
    )
    cache = {0: 0}

    def predict(delta):
        syndrome = sum(
            ((delta & m).bit_count() % 2) << i for i, m in enumerate(detectors)
        )
        if syndrome not in cache:
            e = dec.decode(
                np.fromiter(
                    ((syndrome >> i) & 1 for i in range(len(detectors))), dtype=np.uint8
                )
            )
            out = 0
            for j in np.flatnonzero(e):
                out ^= model.mech_observables[j]
            cache[syndrome] = out
        return cache[syndrome]

    return predict
