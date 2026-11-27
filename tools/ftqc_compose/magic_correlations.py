"""Finite joint-mask propagation for ideal CCZ twirls and independent batches."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, product

from ftqc_magic.delivery import read_delivery
from .supports import fields


def ccz_twirl(masses, blocks):
    """Recheck each GL(3,2) permutation, its Clifford compensation and Z action."""
    if not 1 <= blocks <= 2:
        raise ValueError("joint twirl supports one or two CCZ states")
    counts = [Counter() for _ in range(8)]
    for a, b, c in permutations(range(1, 8), 3):
        if a ^ b == c:
            continue
        perm = [
            (a if x & 1 else 0) ^ (b if x & 2 else 0) ^ (c if x & 4 else 0)
            for x in range(8)
        ]
        # CCZ phase discrepancy is a Boolean polynomial of degree <=2
        # exactly when its coefficient on xyz (parity of the table) vanishes.
        phases = [(x == 7) ^ (perm[x] == 7) for x in range(8)]
        if sum(phases) % 2:
            raise ValueError("non-Clifford twirl correction")
        for mask in range(8):
            image = next(
                z
                for z in range(8)
                if all(
                    (x & mask).bit_count() % 2 == (perm[x] & z).bit_count() % 2
                    for x in range(8)
                )
            )
            counts[mask][image] += 1
    result = defaultdict(Fraction)
    for key, prob in masses.items():
        masks = [counts[(int(key) >> (3 * i)) & 7] for i in range(blocks)]
        for terms in product(*(m.items() for m in masks)):
            out = 0
            mass = Fraction(prob)
            for i, (z, count) in enumerate(terms):
                out |= z << (3 * i)
                mass *= Fraction(count, 168)
            result[out] += mass
    if sum(result.values()) != sum(map(Fraction, masses.values())):
        raise ValueError("joint twirl lost mass")
    return {str(z): str(p) for z, p in sorted(result.items())}


def interleaved(spec):
    """One named T member from each of several independently operated factories."""
    fields(spec, "kind", "instances", "independent_noise")
    if spec["independent_noise"] is not True or not 1 <= len(spec["instances"]) <= 12:
        raise ValueError(
            "requires independent source noise and 1 to 12 distinct instances"
        )
    ids = set()
    data = set()
    proofs = []
    law = {0: Fraction(1)}
    acceptance = Fraction(1)
    for i, row in enumerate(spec["instances"]):
        fields(row, "id", "batch", "member", "factory", "receipt", "data")
        if (
            any(
                not isinstance(row[k], str) or not row[k]
                for k in ("id", "batch", "data")
            )
            or row["id"] in ids
            or row["data"] in data
        ):
            raise ValueError(
                "interleaving needs distinct factory instances and data ports"
            )
        ids.add(row["id"])
        data.add(row["data"])
        src = read_delivery(row["factory"], row["receipt"])
        member = row["member"]
        if (
            src["state"] != "T"
            or src["law"] != "exact-accepted-Z"
            or type(member) is not int
            or not 0 <= member < len(src["groups"])
        ):
            raise ValueError("expected one verified member from an exact T batch")
        (wire,) = src["groups"][member]
        p = sum(
            (
                Fraction(v)
                for m, v in src["conditional_z_mass"].items()
                if int(m) >> wire & 1
            ),
            Fraction(0),
        )
        law = {
            z | (b << i): v * f
            for z, v in law.items()
            for b, f in enumerate((1 - p, p))
        }
        acceptance *= Fraction(src["acceptance_lower"])
        proofs.append(
            dict(
                instance=row["id"],
                batch=row["batch"],
                member=member,
                data=row["data"],
                source=src,
            )
        )
    from .magic_batch import injection

    return dict(
        status="proved",
        instances=proofs,
        instrument=injection("T"),
        all_accepted=str(acceptance),
        conditional_data_Z={str(k): str(v) for k, v in law.items()},
        premise="independent factory instances including acceptance, independent of the joint data/reference; ideal consumers",
    )
