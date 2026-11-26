"""Accepted diagonal-Pauli contract for the shipped eight-qubit HGP instrument."""

import json
from collections import Counter
from functools import cache, lru_cache
from itertools import combinations
from pathlib import Path

import ftqc_extractors
from ftqc_audit.model import propagate
from ftqc_extractors.__main__ import digest
from ftqc_extractors.algebra import (
    css_dimension,
    echelon,
    hgp,
    parity,
    reduce_row,
    seed,
)
from ftqc_extractors.audit_bridge import Bridge

from .contracts import BIT, SCOPE, Violation, contract_result, theorem
from .supports import Family, fields


class Instrument:
    def __init__(self, variant):
        if variant not in ("edge-extractor", "edge-hook"):
            raise ValueError("unknown one-edge HGP schedule")
        source = (
            Path(ftqc_extractors.__file__).parent / "examples" / (variant + ".json")
        )
        self.spec = json.loads(source.read_text())
        self.bridge = Bridge(self.spec)
        self.schedule = self.bridge.ft  # Constructor checks the ideal Heisenberg flows.
        n, hx, hz = hgp(*seed(self.spec["a"]), *seed(self.spec["b"]))
        self.n = n
        self.sites = len(self.schedule.ops)
        self.data = (1 << n) - 1
        self.measured = self.spec["measured_z"]
        # Logical (X1,Z1), (X2,Z2) in the fixture's coordinate order.
        self.logical = [(3, self.measured), (80, 48)]
        if n != 8 or css_dimension((n, hx, hz)) != 2:
            raise ValueError("unexpected HGP logical dimension")
        for i, (x, z) in enumerate(self.logical):
            if any(parity(x, c) for c in hz) or any(parity(z, c) for c in hx):
                raise ValueError("logical basis has syndrome")
            if any(parity(x, b) != (i == j) for j, (_, b) in enumerate(self.logical)):
                raise ValueError("incorrect logical commutation relations")
        bx, bz = echelon(hx), echelon(hz + [self.measured])
        self.rx = [reduce_row(x, bx) for x in range(1 << n)]
        self.rz = [reduce_row(z, bz) for z in range(1 << n)]
        self.sx = [any(parity(x, h) for h in hz) for x in range(1 << n)]
        self.sz = [any(parity(z, h) for h in hx) for z in range(1 << n)]
        self.detectors = [sum(1 << i for i in det) for det in self.spec["detectors"]]
        self.readout = sum(1 << i for i in self.spec["readout"])
        self.columns = {}
        for t, modes in enumerate(self.schedule.modes):
            for m in modes:
                x, z, r = self.bridge.transfer(t, m)
                self.columns[t, m] = x & self.data, z & self.data, r

    def rejected(self, state):
        x, z, r = state
        return self.sx[x] or self.sz[z] or any(parity(r, h) for h in self.detectors)

    def output(self, state):
        x, z, r = state
        return {
            "readout_flip": parity(r, self.readout),
            "survivor_X_flip": parity(z, 80),
            "survivor_Z_flip": parity(x, 48),
            "residual_mod_branch_stabilizers": [self.rx[x], self.rz[z]],
        }

    def leaf(self, instance):
        counts = Counter()
        for fault, state in self.columns.items():
            sig = self.schedule.signature(dict([fault]))
            rejected = self.rejected(state)
            out = self.output(state)
            harm = bool(
                out["readout_flip"] or any(out["residual_mod_branch_stabilizers"])
            )
            if rejected != bool(sig & self.schedule.mask) or (
                not rejected and harm != bool(sig >> len(self.schedule.det))
            ):
                raise ValueError("forward residual disagrees with ideal logical tests")
            if not rejected and harm:
                raise Violation(
                    dict(
                        instance=instance,
                        schedule_sha256=digest(self.spec),
                        op=self.schedule.ops[fault[0]],
                        fault=list(fault),
                        detectors=[],
                        **out,
                        audit=self.bridge.replay([list(fault)]),
                    )
                )
            counts["rejected" if rejected else "accepted_correct"] += 1
        coords = list(map(str, range(self.n)))
        memory = {
            "kind": "Q",
            "wires": coords,
            "code": "one-edge-HGP-" + digest([self.spec["a"], self.spec["b"]]),
            "basis": f"Pauli(X,Z)={self.logical}",
            "decoder": "ideal zero syndrome; fixed logical basis",
            "bad": [[q] for q in coords],
        }
        sites = [f"{instance}/{i}" for i in range(self.sites)]
        return theorem(
            [memory],
            [memory, BIT],
            "nondestructive-Z1",
            Family(sites, combinations(sites, 2)),
            None,
            {
                "schedule_sha256": digest(self.spec),
                "ideal": "checked signed detector/readout/commutant flows",
                "single_modes": dict(counts),
                "measured_z": self.measured,
                "surviving_xz": list(self.logical[1]),
                "terminal_syndrome": "ideal",
                "residual": "identity modulo memory stabilizers and measured Z on each accepted branch",
            },
        )


@lru_cache(maxsize=2)
def instrument(variant):
    return Instrument(variant)


def hgp_leaf(node):
    fields(node, "op", "instance", "variant")
    if not isinstance(node["instance"], str) or not node["instance"]:
        raise ValueError("expected a named HGP instance")
    return instrument(node["variant"]).leaf(node["instance"])


def xor_state(a, b):
    return tuple(x ^ y for x, y in zip(a, b, strict=True))


class RepeatedCircuit:
    """Independent flat Pauli simulation with classical controls rebased at the join."""

    def __init__(self, g):
        self.g = g
        self.cut = g.sites
        self.records = len(g.schedule.records)
        self.ops = list(g.schedule.ops)
        for op in g.schedule.ops:
            self.ops.append([*op[:2], op[2] + self.records] if op[0] == "Z_IF" else op)

    def run(self, faults):
        x = z = records = count = 0
        snapshots = []
        for t, (op, q, *tail) in enumerate(self.ops):
            b = 1 << q
            m = faults.get(t, 0)
            if op in ("RZ", "RX"):
                x &= ~b
                z &= ~b
                if m:
                    if op == "RZ":
                        x |= b
                    else:
                        z |= b
            elif op in ("MZ", "MX"):
                flip = bool((x if op == "MZ" else z) & b) ^ bool(m)
                records |= int(flip) << count
                count += 1
                x &= ~b
                z &= ~b
            elif op in ("CX", "Z_IF"):
                if op == "CX":
                    target = 1 << tail[0]
                    if x & b:
                        x ^= target
                    if z & target:
                        z ^= b
                    if m & 4:
                        x ^= target
                    if m & 8:
                        z ^= target
                elif records >> tail[0] & 1:
                    z ^= b
                if m & 1:
                    x ^= b
                if m & 2:
                    z ^= b
            else:
                raise ValueError("unsupported flat HGP operation")
            if (t + 1) % self.cut == 0:
                snapshots.append(
                    (
                        x & self.g.data,
                        z & self.g.data,
                        (records >> (count - self.records)) & ((1 << self.records) - 1),
                    )
                )
        return snapshots


def repeated_result(spec):
    fields(spec, "kind", "variant")
    g = instrument(spec["variant"])
    a = {"op": "hgp-z", "instance": "first", "variant": spec["variant"]}
    b = dict(a, instance="second")
    program = {
        "op": "sequence",
        "left": a,
        "right": {"op": "parallel", "left": b, "right": {"op": "bit-wire"}},
    }
    flat = RepeatedCircuit(g)
    try:
        composition = contract_result(
            {
                "kind": "contract",
                "scope": SCOPE,
                "program": program,
                "require_closed": True,
                "require_friendly": False,
            }
        )
    except Violation as failure:
        t, m = failure.evidence["fault"]
        s1, s2 = flat.run({t: m})
        final = g.output(s2)
        if (
            g.rejected(s1)
            or g.rejected(s2)
            or not (
                parity(s1[2], g.readout)
                or final["readout_flip"]
                or any(final["residual_mod_branch_stabilizers"])
            )
        ):
            raise ValueError("hook failed whole-circuit replay")
        return {
            "status": "counterexample",
            "scope": SCOPE,
            "witness": failure.evidence,
            "flat_check": {
                "faults": [[t, m]],
                "detectors": [],
                "first": g.output(s1),
                "second": final,
            },
        }
    labels = [f"{stage}/{t}" for stage in ("first", "second") for t in range(g.sites)]
    cover = {frozenset(t) for t in composition["bad_faults"]["minimal"]}
    if set(labels) != set(composition["bad_faults"]["sites"]) or any(
        len(t) != 2 for t in cover
    ):
        raise ValueError("unexpected composed location domain or cover")

    zero = (0, 0, 0)

    # Independent prediction from audit columns for each gadget and propagation
    # of the first output Pauli through the fault-free second gadget.
    @cache
    def across(x, z):
        ax, az, r = propagate(g.bridge.circuit, 0, x, z)
        return ax & g.data, az & g.data, r

    counts = Counter()
    witness = None
    first_sites = g.sites
    modes = list(g.columns)
    locations = [(stage * first_sites + i, m) for stage in range(2) for i, m in modes]
    for w in (0, 1, 2):
        for faults in combinations(locations, w):
            if w == 2 and faults[0][0] == faults[1][0]:
                continue
            local = [zero, zero]
            for i, m in faults:
                stage, t = divmod(i, first_sites)
                local[stage] = xor_state(local[stage], g.columns[t, m])
            s1, s2 = flat.run(dict(faults))
            if s1 != local[0] or s2 != xor_state(across(*local[0][:2]), local[1]):
                raise ValueError(f"flat and joined boundary states differ at {faults}")
            accepted = not (g.rejected(s1) or g.rejected(s2))
            bad = accepted and (
                parity(s1[2], g.readout)
                or parity(s2[2], g.readout)
                or g.rx[s2[0]]
                or g.rz[s2[1]]
            )
            if bad and frozenset(labels[i] for i, m in faults) not in cover:
                raise ValueError("harmful execution outside composed bad-support cover")
            if bad and witness is None:
                witness = {
                    "faults": [list(f) for f in faults],
                    "first": g.output(s1),
                    "second": g.output(s2),
                }
            counts[f"weight_{w}"] += 1
            if bad:
                counts["accepted_harmful"] += 1
    if witness is None:
        raise ValueError("missing upper witness at weight two")
    return {
        "status": "proved",
        "scope": SCOPE,
        "composition": composition,
        "flat_check": {
            "locations": len(flat.ops),
            "counts": dict(counts),
            "distance": len(witness["faults"]),
            "witness": witness,
            "boundaries": "ideal input and ideal syndrome test after each measurement; noiseless bit storage",
        },
    }
