"""Replay the straight-line H2-2 QIR and shot histograms from Zenodo 23165109."""

from collections import Counter
from hashlib import sha256
import json, math, re
from pathlib import Path


def parse(text, p=0):
    import stim

    circuit = stim.Circuit()
    measured = []
    outputs = []
    registers = {}
    ops = Counter()
    array_size = None
    if re.search(r"^\s*(?:br|switch|invoke)\b", text, re.M):
        raise ValueError("expected a straight-line circuit")
    for line in text.splitlines():
        if not re.search(r"\bcall\b", line):
            continue
        match = re.search(r"call void @__quantum__(\w+)\(", line)
        if not match:
            raise ValueError("unsupported QIR call")
        op = match[1]
        qs = [int(i) for i in re.findall(r"i64 (\d+) to %Qubit\*", line)]
        ids = [int(i) for i in re.findall(r"i64 (\d+) to %Result\*", line)]
        if op == "rt__initialize":
            continue
        if op == "rt__array_record_output":
            array_size = int(re.search(r"output\(i64 (\d+)", line)[1])
            continue
        if op == "rt__result_record_output":
            outputs.append(registers[ids[0]])
            continue
        ops[op] += 1
        if op == "qis__m__body":
            registers[ids[0]] = len(measured)
            measured.append(qs[0])
            circuit.append("M", qs, p)
            continue
        if op == "qis__reset__body":
            circuit.append("R", qs)
        elif op in ("qis__rx__body", "qis__rz__body", "qis__rzz__body"):
            theta = float(re.search(r"double ([^,]+)", line)[1]) / (math.pi / 2)
            count = round(theta)
            if abs(count - theta) > 1e-10:
                raise ValueError("expected Clifford native rotation")
            table = {
                "qis__rx__body": ["I", "SQRT_X", "X", "SQRT_X_DAG"],
                "qis__rz__body": ["I", "S", "Z", "S_DAG"],
                "qis__rzz__body": ["I", "SQRT_ZZ", "Z", "SQRT_ZZ_DAG"],
            }
            circuit.append(table[op][count % 4], qs)
        else:
            raise ValueError("unknown QIR opcode " + op)
        if p:
            circuit.append("DEPOLARIZE2" if len(qs) == 2 else "DEPOLARIZE1", qs, p)
    if len(outputs) != array_size or sorted(outputs) != list(range(len(measured))):
        raise ValueError("expected permutation of every acquired bit")
    return circuit, measured, outputs, dict(ops)


def case_files(folder, name):
    qi = folder / (name + "-input.json")
    hist = folder / (name + "-output.json")
    data = json.loads(hist.read_text())["Results"][0]
    circuit, wires, order, ops = parse(qi.read_text())

    def mask(row):
        if len(row) != len(order) or any(v not in (0, 1) for v in row):
            raise ValueError("invalid shot")
        return sum(v << order[j] for j, v in enumerate(row))

    observed = Counter({mask(r["Outcome"]): r["Count"] for r in data["Histogram"]})
    if observed != Counter(mask(row) for row in data["Shots"]):
        raise ValueError("shot count mismatch")
    return (
        circuit,
        wires,
        order,
        observed,
        dict(
            source_sha256={
                f.name: sha256(f.read_bytes()).hexdigest() for f in (qi, hist)
            },
            operations=ops,
            physical_qubits=circuit.num_qubits,
            acquisitions=len(wires),
            shots=sum(observed.values()),
        ),
    )


def logical_parities(name, wires):
    # QIR numbers q0..q9 of Fig. 1 as physical wires 5..14.
    def last(sites):
        return [len(wires) - 1 - wires[::-1].index(q) for q in sites]

    def gadget(sites, occurrence):
        ids = [[i for i, q in enumerate(wires) if q == j] for j in sites]
        return [group[occurrence] for group in ids]

    if name.startswith("SEQ"):
        # Second gadget after the three kappa qubits have been released and reset.
        data = last([7, 11, 12])
        fusion = gadget([20], 0)
        return dict(data=data + fusion, gadget=gadget([24, 25, 26], -1) + fusion)
    read = last([7, 8, 12])
    spectator = last([6, 12, 14])
    return dict(data=read, spectator=spectator, gadget=gadget([24, 25, 26], 0))


def replay(folder, name, p, shots, seed):
    import stim
    from scipy.stats import beta

    ideal, wires, order, observed, report = case_files(folder, name)
    probe, _, _, _ = parse((folder / (name + "-input.json")).read_text(), p)
    parities = logical_parities(name, wires)

    def bits(ids):
        return sum(1 << j for j in ids)

    det = []
    measurements = ideal.num_measurements

    def deterministic(ids):
        t = ideal.copy()
        t.append("DETECTOR", [stim.target_rec(j - measurements) for j in ids])
        try:
            t.detector_error_model()
        except ValueError:
            return False
        return True

    # Accept trivial initial eigenvalues for deterministic first-cycle checks.
    ref = ideal.reference_sample()
    for j in range(1, 11):
        if deterministic([j]):
            det.append(j)
    sim = probe.compile_sampler(seed=seed).sample(shots, bit_packed=True)
    counts = Counter(int.from_bytes(row.tobytes(), "little") for row in sim)
    from .ion_decoder import detector_masks, predictor

    active = {k: ids for k, ids in parities.items() if deterministic(ids)}
    comparison = sorted(set(parities["gadget"]) ^ set(parities["data"]))
    if deterministic(comparison):
        active["gadget_data_disagreement"] = comparison
    detectors = detector_masks(wires, det, deterministic, name.startswith("SEQ"))
    predict = predictor(probe.copy(), detectors, list(active.values()))
    reference = bits([i for i, value in enumerate(ref) if value])

    def stats(hist, ids, logical=None):
        keep = fail = 0
        expected = sum(int(ref[i]) for i in ids) % 2
        sig = bits(ids)
        for value, n in hist.items():
            if any((value >> j & 1) != ref[j] for j in det):
                continue
            keep += n
            correction = (
                0 if logical is None else (predict(value ^ reference) >> logical) & 1
            )
            fail += (((value & sig).bit_count() % 2) ^ expected ^ correction) * n
        return dict(
            accepted=keep,
            errors=fail,
            jeffreys95=[
                float(x)
                for x in beta.ppf([0.025, 0.975], fail + 0.5, keep - fail + 0.5)
            ],
        )

    results = {}
    for label, ids in dict(
        parities, **{k: v for k, v in active.items() if k not in parities}
    ).items():
        if deterministic(ids):
            j = list(active).index(label)
            results[label] = dict(
                acquisition_indices=ids,
                raw_data=stats(observed, ids),
                simulated_raw=stats(counts, ids),
                decoded_data=stats(observed, ids, j),
                simulated_decoded=stats(counts, ids, j),
            )
        else:
            results[label] = dict(
                acquisition_indices=ids,
                status="nondeterministic without decoding or fusion frame",
            )
    return dict(
        **report,
        case=name,
        seed=seed,
        p=p,
        simulated_shots=shots,
        stim=stim.__version__,
        decoder="BP+OSD-0, 20 iterations",
        detector_masks=detectors,
        deterministic_preparation_checks=det,
        observables=results,
    )


def main():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", type=Path)
    ap.add_argument("--p", type=float, default=0.001)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--shots", type=int, default=4000)
    args = ap.parse_args()
    results = [
        replay(args.folder, q.name[:-11], args.p, args.shots, args.seed)
        for q in sorted(args.folder.glob("*-input.json"))
        if not q.name.startswith("MEM")
    ]
    if not results:
        raise ValueError("no surgery data found")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
