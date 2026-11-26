"""Cross-check ftqc_extractors location faults with ftqc_audit's independent Clifford transport."""

import argparse
import json
import sys
from pathlib import Path

from ftqc_audit.analysis import min_undetected
from ftqc_audit.circuit import Circuit
from ftqc_audit.model import propagate

from .__main__ import digest, validate
from .algebra import bits
from .circuit import Schedule
from .t_faults import FaultCircuit

ARMS = {"gate": None, "X": 0, "Y": 1}
NAMES = {"RZ": "R", "MZ": "M"}


def schedule(spec, arm):
    validate(spec)
    if spec["kind"] == "circuit" and arm == "gate":
        return Schedule(spec), False
    if spec["kind"] in ("t_gate", "t_bounds") and spec["resource"] == "ideal-encoded-T":
        return FaultCircuit(spec, ARMS[arm]), True
    raise ValueError("requires a circuit or a fixed arm of a conditional T protocol")


class Bridge:
    def __init__(self, spec, arm="gate"):
        self.ft, two_blocks = schedule(spec, arm)
        inputs = self.ft.n * (2 if two_blocks else 1)
        # ftqc_audit's Circuit requires each wire to be introduced. These ideal
        # placeholders precede every fault. Only the input-independent Pauli
        # transfer function is used; their state is never used as the HGP input.
        ticks = [[("!R", q) for q in range(inputs)]]
        ticks.extend([(NAMES.get(g, g), *args)] for g, *args in self.ft.ops)
        self.circuit = Circuit(ticks, detectors=[], outputs=[])
        self.records = {
            tick - 1: j for j, (tick, _, _, _) in enumerate(self.circuit.meas)
        }
        self.tests = [
            (x, z, sum(1 << self.records[t] for t in bits(rec)))
            for x, z, rec in self.ft.flows
        ]
        self.locations = self.circuit.locations()

    def transfer(self, site, mode):
        ft = self.ft
        if (
            type(site) != int
            or not 0 <= site < len(self.locations)
            or type(mode) != int
            or mode not in ft.modes[site]
        ):
            raise ValueError("invalid ftqc_extractors fault location or mode")
        tick, _, op = self.locations[site]
        x = z = rec = 0
        if op.name in ("MX", "M"):
            rec = 1 << self.records[site]
        else:
            if op.name == "RX":
                mode = 2
            for i, q in enumerate(op.qs):
                if mode >> (2 * i) & 1:
                    x |= 1 << q
                if mode >> (2 * i + 1) & 1:
                    z |= 1 << q
        return propagate(self.circuit, tick, x, z, rec)

    def signature(self, faults):
        if len({i for i, m in faults}) != len(faults):
            raise ValueError("each faulty site must occur once")
        x = z = rec = 0
        for site, mode in faults:
            a, b, c = self.transfer(site, mode)
            x ^= a
            z ^= b
            rec ^= c
        return sum(
            (
                (
                    (px & z).bit_count()
                    + (pz & x).bit_count()
                    + (records & rec).bit_count()
                )
                % 2
            )
            << i
            for i, (px, pz, records) in enumerate(self.tests)
        )

    def replay(self, faults):
        sig = self.signature(faults)
        if sig != self.ft.signature(dict(faults)):
            raise ValueError("ftqc_audit and ftqc_extractors disagree on the witness")
        return {
            "faults": faults,
            "detectors": list(bits(sig & self.ft.mask)),
            "logical_tests": list(bits(sig >> len(self.ft.det))),
            "undetected_harm": bool(sig and not sig & self.ft.mask),
        }

    def compare(self, weight, cap):
        if type(weight) is not int or type(cap) is not int or weight < 1 or cap < 1:
            raise ValueError("positive search limits required")
        # ftqc_extractors uses reverse propagation. ftqc_audit simulates forward through
        # its physical circuit and we evaluate the shared endpoint predicate.
        columns = {}
        count = 0
        for site, group in enumerate(self.ft.events(reverse=True)):
            tick, _, op = self.locations[site]
            for mode, expected in group:
                observed = self.signature([[site, mode]])
                if observed != expected:
                    raise ValueError(
                        f"fault signatures disagree at {site}, mode {mode}"
                    )
                if observed:
                    key = (observed & self.ft.mask, observed >> len(self.ft.det))
                    columns.setdefault(key, (site, mode))
                count += 1
            # The other two ftqc_audit preparation Paulis should be either
            # harmless or equivalent to ftqc_extractors's one eigenstate flip.
            if op.name in ("R", "RX"):
                q = op.qs[0]
                x, z, rec = propagate(
                    self.circuit,
                    tick,
                    1 << q if op.name == "RX" else 0,
                    1 << q if op.name == "R" else 0,
                )
                if any(
                    (
                        (px & z).bit_count()
                        + (pz & x).bit_count()
                        + (r & rec).bit_count()
                    )
                    % 2
                    for px, pz, r in self.tests
                ):
                    raise ValueError(
                        f"preparation stabilizer has a nonzero signature at {site}"
                    )
        keys = sorted(columns)
        result = min_undetected(keys, any, weight, cap)
        selected = result.pop("witness")
        if selected is not None:
            sites = {}
            for index in selected:
                site, mode = columns[keys[index]]
                sites[site] = sites.get(site, 0) ^ mode
            faults = [[i, m] for i, m in sorted(sites.items()) if m]
            witness = self.replay(faults)
            if not witness["undetected_harm"] or len(faults) != result["upper"]:
                raise ValueError("distance witness violates the shared harm predicate")
            result["witness"] = witness
        return {
            "max_weight": weight,
            "cap": cap,
            "locations": len(self.locations),
            "checked_modes": count,
            "distinct_nonzero_signatures": len(keys),
            "distance": result,
        }


def cat_prefix(spec, arm, round_index):
    ft, two_blocks = schedule(spec, arm)
    if not two_blocks or not 0 <= round_index < len(ft.parts):
        raise ValueError("select a cat round of an ftqc_extractors T schedule")
    start, end, x, z, _, flags = ft.parts[round_index]
    size = (x | z).bit_count()
    # Singleton cat has no verification flags and consists of one preparation.
    stop = flags[-1] + 1 if flags else start + 1
    first = 2 * ft.n
    if not start < stop <= end:
        raise ValueError("invalid cat round boundary")
    ticks = []
    for g, *qs in ft.ops[start:stop]:
        if any(q < first for q in qs):
            raise ValueError("data interaction before end of cat preparation")
        ticks.append([(NAMES.get(g, g), *(q - first for q in qs))])
    c = Circuit(
        ticks,
        detectors=[[j] for j in range(len(flags))],
        outputs=list(range(size)),
        name=f"FTQC {arm} round {round_index} cat preparation",
    )
    native_measurements = [start + t for t, _, _, _ in c.meas]
    if native_measurements != flags:
        raise ValueError("cat flag record mismatch")
    return dict(
        c.to_json(),
        ftqc_origin={
            "spec_sha256": digest(spec),
            "arm": arm,
            "round": round_index,
            "first_location": start,
            "first_qubit": first,
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="cmd", required=True)
    for name in ("compare", "replay", "export-cat"):
        p = subs.add_parser(name)
        p.add_argument("spec", type=Path)
        p.add_argument("--arm", choices=ARMS, default="gate")
        p.add_argument("-o", type=Path)
        if name == "compare":
            p.add_argument("--max-weight", type=int, default=2)
            p.add_argument("--cap", type=int, default=3000000)
        elif name == "export-cat":
            p.add_argument("--round", type=int, default=0)
        else:
            p.add_argument(
                "--fault", action="append", required=True, metavar="LOCATION:MODE"
            )
    args = parser.parse_args()
    try:
        s = json.loads(args.spec.read_text())
        if args.cmd == "export-cat":
            out = cat_prefix(s, args.arm, args.round)
        else:
            bridge = Bridge(s, args.arm)
            if args.cmd == "compare":
                out = bridge.compare(args.max_weight, args.cap)
            else:
                out = bridge.replay([list(map(int, v.split(":"))) for v in args.fault])
            out.update(spec_sha256=digest(s), arm=args.arm)
        txt = json.dumps(out, sort_keys=True, indent=2) + "\n"
        if args.o:
            args.o.write_text(txt)
        else:
            print(txt, end="")
        return 0
    except (IndexError, KeyError, TypeError, ValueError, OSError) as exc:
        print(f"cross-check refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
