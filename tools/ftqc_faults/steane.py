"""Compare two physical single-location engines on the scheduled Steane extractor."""

import hashlib
import json
from itertools import product
from types import SimpleNamespace

from ftqc_audit.circuit import Circuit
from ftqc_audit.library import css_prep
from ftqc_audit.record_model import RecordModel
from ftqc_zx import (check_fault_equiv, make_steane_7_1_3, schedule_circuit_ticks,
                     synth_steane_syndrome, synthesize_flag_decoder, verify_derivation)
from ftqc_zx.contracts import propagate_single_fault

NATIVE = {"CX": "CX", "IDLE": "I", "PREP_PLUS": "RX", "PREP_0": "R",
          "MEAS_X": "MX", "MEAS_Z": "M"}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def closed_probe(schedule, code, rows, basis):
    # Use the numbering in the actual ZX code. All boundary locations are ideal.
    prep = css_prep(SimpleNamespace(n=code.n, hx=code.hx_rows, name="steane"))
    ticks = [[("!" + op.name, *op.qs) for op in tick] for tick in prep.ticks]
    ticks.extend([[(NATIVE[op.kind], *op.qubits) for op in tick] for tick in schedule.tick_ops])
    ticks.append([("!M" if basis == "Z" else "!MX", i) for i in range(code.n)])
    measurements = [op for op in schedule.locations if op.meas_index >= 0]
    indices = {op.label: i for i, op in enumerate(measurements)}
    # Audit records are chronological, while ZX uses the synthesis record order.
    raw = [indices[op.label] for op in sorted(measurements, key=lambda op: op.meas_index)]
    flag, = rows.flag_meas_labels
    detectors = [[indices[flag]]]
    for row in rows.h_prime_rows:
        detectors.append([indices[label] for label, bit in zip(rows.syndrome_meas_labels, row) if bit])
    circ = Circuit(ticks, detectors=detectors, observables=[], outputs=[], name="steane-probe-" + basis)
    return RecordModel(circ, response_records=raw + list(range(len(raw), len(raw) + code.n)))


def compare_axis(sector):
    code = make_steane_7_1_3()
    circuit, proof = synth_steane_syndrome(code, mode="optimised_5anc", sector=sector)
    equiv = check_fault_equiv(circuit, proof.spec_circuit, w=3, noise_model_lhs="circuit_css", codespace=code)
    if not equiv.equivalent or not equiv.spackle_crosschecked or verify_derivation(proof)["verdict"] != "PROOF-VALID":
        raise ValueError("Steane ZX equivalence or rewrite proof failed")
    decoder = synthesize_flag_decoder(circuit, code, extracted_sector=sector)
    schedule = schedule_circuit_ticks(circuit)
    probes = [closed_probe(schedule, code, decoder, basis) for basis in ("Z", "X")]
    rows = []
    m = schedule.meas_count
    for loc, *groups in zip(schedule.locations, *(p.groups for p in probes), strict=True):
        is_meas = loc.meas_index >= 0
        modes = ["flip"] if is_meas else ["".join(ps) for ps in product("IXYZ", repeat=len(loc.qubits)) if any(p != "I" for p in ps)]
        if any(set(group) != set(modes) for group in groups):
            raise ValueError(f"unmatched physical modes at {loc}")
        for label in modes:
            paulis = ("Z" if loc.kind == "MEAS_X" else "X") if is_meas else label
            x = sum(1 << q for q, p in zip(loc.qubits, paulis, strict=True) if p in "XY")
            z = sum(1 << q for q, p in zip(loc.qubits, paulis, strict=True) if p in "ZY")
            out_x, out_z, bits = propagate_single_fault(schedule, circuit.output_qubits, loc.loc_id, x, z)
            # ZX raw bit 0 is the flag, followed by the four syndrome records.
            detectors = bits & 1
            for i, row in enumerate(decoder.h_prime_rows, 1):
                parity = sum(bit << j for j, bit in enumerate(row))
                detectors |= ((parity & (bits >> 1)).bit_count() % 2) << i
            for probe, output in zip(probes, (out_x, out_z)):
                actual = probe.replay([(loc.loc_id, label)])
                if actual != (detectors, bits | (output << m)):
                    raise ValueError(f"single-fault disagreement: sector={sector}, tick={loc.tick}, site={loc.loc_id}, mode={label}")
            rows.append([loc.tick, loc.kind, loc.qubits, label, detectors, bits, out_x, out_z])
    return {"sector": sector, "equivalence": equiv.to_dict(), "derivation": "PROOF-VALID",
            "ticks": schedule.ticks, "locations": schedule.total_locations, "cnots": schedule.cnot_count,
            "idle_locations": schedule.idle_location_count, "single_fault_modes": len(rows),
            "response_sha256": digest(rows), "probe_sha256": [p.circuit.sha256() for p in probes],
            "detectors": probes[0].detectors, "record_labels": [op.label for op in circuit.ops if op.kind.startswith("MEAS_")]}


def receipt():
    return {"schema": "ftqc_faults-steane-1", "input": "ideal Steane logical zero",
            "schedule": "ZX ASAP with live-wire idles", "frame": "raw branches before recovery",
            "physical_model": "one noisy location; joint Pauli after gate, preparation or idle; readout bit flip",
            "axes": [compare_axis(sector) for sector in ("X", "Z")]}
