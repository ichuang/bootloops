"""Signed stabilizer inputs for code-action replay, including non-CSS codes."""

import json
from dataclasses import dataclass

from . import basis, codes


@dataclass
class StabilizerCode:
    name: str
    n: int
    signed_checks: tuple
    slots: list

    def commutes(self):
        return all(
            basis.commute_pauli(a, b)
            for i, a in enumerate(self.signed_checks)
            for b in self.signed_checks[:i]
        )

    def sha256(self):
        return basis._canonical_sha256(
            {"n": self.n, "signed_checks": self.signed_checks}
        )


def load(spec):
    if spec.endswith(".json"):
        with open(spec) as stream:
            data = json.load(stream)
        if isinstance(data, dict) and data.get("schema") == "ftqc_codes-stabilizer-1":
            n = data.get("n")
            if type(n) is not int or n < 1 or not isinstance(data.get("checks"), list):
                raise codes.SpecError("expected n > 0 and signed stabilizer checks")
            if not isinstance(data.get("slots"), list):
                raise codes.SpecError("expected ordered logical slots")
            checks = tuple(basis.parse_pauli(p, n) for p in data["checks"])
            return StabilizerCode(
                data.get("name", "stabilizer"), n, checks, data["slots"]
            )
    return codes.load(spec)


def default_basis(code):
    if isinstance(code, StabilizerCode):
        return basis.build_basis_receipt(code, code.slots)
    return basis.derive_paired_basis(code)
