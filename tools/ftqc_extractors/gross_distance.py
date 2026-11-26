"""Generate and verify static BB144 memory and unprimed-Z merged distance receipts."""

from ftqc_codes.codes import CSSCode

from . import css_cover
from .algebra import bits, echelon, parity, reduce_row
from .gross import LEFT_Z, RIGHT_Z, cell, gross_z, logical_action, shifted_z


def witness_mask(sites, n, checks, stabilizers):
    if (
        not isinstance(sites, list)
        or not 0 < len(sites) <= 12
        or any(type(q) is not int or not 0 <= q < n for q in sites)
        or sorted(set(sites)) != sites
    ):
        raise ValueError("invalid gross-code logical support")
    word = sum(1 << q for q in sites)
    if any(parity(word, h) for h in checks) or not reduce_row(
        word, echelon(stabilizers)
    ):
        raise ValueError("upper witness is not a nontrivial logical")
    return word


def gross_certificate(spec, budget, saved=None):
    memory, merged, _, _ = gross_z()
    codes = {"memory": (memory.n, memory.hx, memory.hz), "merged": merged}
    # Translation by y gives another Z logical. ZX duality (invert cells
    # and swap the halves) supplies an X logical away from the ports.
    x_word = sum(
        1 << (half * 72 + cell(-i, -j))
        for half, terms in enumerate((RIGHT_Z, LEFT_Z))
        for i, j in terms
    )
    initial = {"X": x_word, "Z": shifted_z(0, 1)}
    records = {}
    for name, (n, hx, hz) in codes.items():
        axes = {}
        for axis, checks, stabs in (("X", hz, hx), ("Z", hx, hz)):
            if saved is None:
                word = initial[axis]
                while True:
                    sites = list(bits(word))
                    witness_mask(sites, n, checks, stabs)
                    shorter, nodes, proof = css_cover.cover(
                        n, checks, stabs, len(sites) - 1, budget
                    )
                    budget -= nodes
                    if not shorter:
                        break
                    word = shorter
            else:
                entry = saved["codes"][name]["axes"][axis]
                sites, nodes, proof = entry["witness"], entry["nodes"], entry["cover"]
                witness_mask(sites, n, checks, stabs)
                css_cover.verify_cover(
                    n, checks, stabs, len(sites) - 1, nodes, proof, budget
                )
                budget -= nodes
            axes[axis] = dict(
                distance=len(sites), witness=sites, nodes=nodes, cover=proof
            )
        records[name] = dict(
            matrix_sha256=CSSCode(name, n, hx, hz).sha256(),
            axes=axes,
            distance=min(a["distance"] for a in axes.values()),
        )
    return dict(
        status="proved"
        if records["merged"]["distance"] >= spec["target"]
        else "counterexample",
        model="static-css",
        proof_format="odd-check-cover-v1",
        action=logical_action(),
        codes=records,
    )
