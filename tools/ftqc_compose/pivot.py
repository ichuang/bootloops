"""Exact signed branches for Xu et al. Fig. 4b in positive-Y convention."""

from itertools import product
from ftqc_codes.gf2 import dot


def branch_table():
    rows = []
    for a, b, c, d in product(range(2), repeat=4):
        y = a ^ c
        # Our Y readout uses |+i> as the positive eigenstate. Convert its
        # physical record to the signed correction convention of the pivot.
        frame = b ^ d ^ (y & (1 ^ c))
        for error in (0, 1):
            phases = []
            for parity in (0, 1):
                pivot = parity ^ c
                raw = pivot ^ a
                phases.append(
                    (
                        raw
                        + 4 * error * raw
                        + 4 * b * raw
                        + 4 * d * pivot
                        - 2 * y * pivot
                        + 4 * frame * parity
                    )
                    % 8
                )
            if (phases[1] - phases[0]) % 8 != 1 + 4 * error:
                raise ValueError("pivot branch has incorrect relative phase")
        rows.append(
            dict(
                outcomes=[a, b, c, d],
                basis="Y" if y else "X",
                frame=frame,
                probability="1/16",
            )
        )
    return rows


def ideal_table(mapping):
    import hashlib

    rows = mapping["rotations"]
    axes = [x for x, z in mapping["basis"]]
    masks = []
    for r in rows:
        mask = sum(dot(x, r["physical_z"]) << i for i, x in enumerate(axes))
        mask ^= 32
        for slot in r["mask_zero_slots"]:
            mask ^= 1 << slot
        if mask != r["parity"]:
            raise ValueError("physical rotation acts on incorrect slots")
        masks.append(mask)
    # Contract the 32 computational columns of the five-plus input.
    correction = mapping["clifford_corrections"]
    if any(v for k, v in correction.items() if k != "output_s_corrections"):
        raise ValueError("unexpected RM Clifford")
    s = correction["output_s_corrections"].get(
        0, correction["output_s_corrections"].get("0")
    )
    phases = [(sum(dot(x, p) for p in masks) + 2 * s * (x & 1)) % 8 for x in range(32)]
    if phases != [x & 1 for x in range(32)]:
        raise ValueError("incorrect output T or nontrivial check phase")
    # A resource Z adds pi on its actual physical parity; compare the
    # resulting four X check signs and output Z to the native oracle.
    from ftqc_magic.catalog import get_catalog_factory

    fac = get_catalog_factory("bk_15_1_3_t")
    signatures = bytearray(32768)
    oracle = bytearray(32768)
    for m in range(1, 32768):
        bit = m & -m
        i = bit.bit_length() - 1
        signatures[m] = signatures[m ^ bit] ^ masks[i]
        oracle[m] = oracle[m ^ bit] ^ (
            fac.column_label(i) | (fac.column_syndrome(i) << 1)
        )
    if signatures != oracle:
        raise ValueError("compiled instrument differs on input masks")
    return dict(
        branches=branch_table(),
        input_masks=len(signatures),
        truth_sha256=hashlib.sha256(signatures).hexdigest(),
        clean_output_phases=phases[:2],
        accept_counts=[signatures.count(0), signatures.count(1)],
    )
