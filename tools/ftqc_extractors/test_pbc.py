"""Independent two-qubit amplitude check of all adaptive records of the bridge fixture."""

import cmath
import json
import math
import unittest
from pathlib import Path

from .__main__ import receipt, verify
from .fixed_graph import pauli
from .algebra import parity

EXAMPLE = Path(__file__).parent / "examples/bridged-pbc.json"


def apply(v, x, z, sign):
    phase = (-1) ** sign * 1j ** (x & z).bit_count()
    return [phase * (-1) ** parity(i ^ x, z) * v[i ^ x] for i in range(len(v))]


def logical_action(setting, b):
    # Restrict the actual physical operator to encoded space by its commutators
    # and signed multiplication (these fixtures have positive CSS logicals).
    w = setting["logical"]
    x = z = 0
    phase = 0
    for pos, idx in enumerate(setting["blocks"]):
        p = pauli("+" + w[1 + 4 * pos : 5 + 4 * pos], 4)
        for j, (xs, zs) in enumerate(b["basis"]):
            a, c = pauli(xs, 4), pauli(zs, 4)
            xb = parity(p.x, c.z) ^ parity(p.z, c.x)
            zb = parity(p.x, a.z) ^ parity(p.z, a.x)
            x |= xb << (2 * idx + j)
            z |= zb << (2 * idx + j)
            phase += xb * zb
        # Physical ordered product X*Z differs from canonical physical Y.
        phase -= (p.x & p.z).bit_count()
    return x, z, (int(w[0] == "-") + phase // 2) % 2


class PBCTests(unittest.TestCase):
    def test_signed_branches(self):
        s = json.loads(EXAMPLE.read_text())
        # H,S before and after the rotations exercise signs at the output.
        s["gates"].extend([["H", [0, 0]], ["S", [1, 0]], ["H", [1, 0]]])
        r = receipt(s)
        verify(s, r)
        out = r["result"]
        self.assertEqual(len(out["branches"]), 16)
        root = math.sqrt(2)
        t = cmath.exp(1j * math.pi / 4)
        initial = [complex(1, 2), complex(-2, 1), complex(3, -1), complex(2, 0)]

        def had(v, site):
            return [
                (v[i ^ (1 << site)] + (-1 if i & (1 << site) else 1) * v[i]) / root
                for i in range(4)
            ]

        ideal = had(initial, 0)
        ideal = [ideal[i ^ (2 if i & 1 else 0)] for i in range(4)]
        ideal = had([a * (t if i & 1 else 1) for i, a in enumerate(ideal)], 0)
        ideal = had([a * (1j if i & 2 else 1) for i, a in enumerate(ideal)], 1)
        for branch in out["branches"]:
            state = initial[:]
            measurements = [
                e for e in branch["program"] if e.get("role") == "rotation parity b"
            ]
            for idx, (e, raw) in enumerate(zip(measurements, branch["outcomes"])):
                x, z, sgn = logical_action(
                    out["settings"][e["measurement"]], s["blocks"][0]
                )

                # Remove the resource Z and contract the measured resource in X.
                def keep_data(mask):
                    return (mask & 1) | ((mask >> 1) & 2)

                a = apply(state, keep_data(x), keep_data(z), sgn)
                amp = (-1j if idx == 0 else t) * (-1 if raw >> 1 else 1)
                plus, minus = (amp, 1) if raw & 1 else (1, amp)
                state = [
                    ((plus + minus) * v + (plus - minus) * p) / 2
                    for v, p in zip(state, a)
                ]
            for measurement in (
                e for e in branch["program"] if e.get("role", "").startswith("output")
            ):
                idx = int(measurement["role"].split("(")[1].split(",")[0])
                x, z, sgn = logical_action(
                    out["settings"][measurement["measurement"]], s["blocks"][0]
                )
                a = apply(
                    state, (x & 1) | ((x >> 1) & 2), (z & 1) | ((z >> 1) & 2), sgn
                )
                actual = sum(p.conjugate() * v for p, v in zip(state, a)).real
                expected = sum(
                    abs(v) ** 2 * (-1 if i & (1 << idx) else 1)
                    for i, v in enumerate(ideal)
                )
                self.assertAlmostEqual(actual, expected)

    def test_missing_or_shared_ports(self):
        s = json.loads(EXAMPLE.read_text())
        with self.assertRaises(ValueError):
            receipt(dict(s, bridges=[]))
        s["bridges"][0]["matching"] = [[0, 0], [0, 1]]
        with self.assertRaises(ValueError):
            receipt(s)

    def test_frame_availability_and_disjoint_rotations(self):
        s = json.loads(EXAMPLE.read_text())
        s["gates"] = [["T", [0, 0]]]
        branch = receipt(s)["result"]["branches"][0]["program"]
        cycles = {op["role"]: op["cycle"] for op in branch if "role" in op}
        self.assertGreater(cycles["output Z(0,0)"], cycles["resource X m"])
        s.update(
            gates=[["T", [0, 0]], ["T", [1, 0]]], magic_supplies=[0, 1], bridges=[]
        )
        branches = receipt(s)["result"]["branches"]
        self.assertEqual(len(branches), 16)
        for b in branches:
            self.assertEqual(
                [
                    op["cycle"]
                    for op in b["program"]
                    if op.get("role") == "rotation parity b"
                ],
                [0, 0],
            )
            self.assertEqual(
                [
                    op["cycle"]
                    for op in b["program"]
                    if op.get("role", "").startswith("output")
                ],
                [2, 2],
            )

    def test_partial_ports_need_not_cover_first_basis_z(self):
        from .block_map import Hardware

        s = json.loads(EXAMPLE.read_text())
        s.update(blocks=s["blocks"][:1], bridges=[])
        s["blocks"][0]["ports"] = [0, 1, -1, -1]
        hw = Hardware(s, 2000000)
        self.assertEqual(hw.select((0,), "+ZZII")["support"], [0, 1])
        with self.assertRaises(ValueError):
            hw.select((0,), "+ZIZI")
