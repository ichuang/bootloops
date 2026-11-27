"""Postselected RM15 driven by accepted, potentially correlated injection batches."""

from fractions import Fraction
from itertools import product
from .catalog import get_catalog_factory
from .certify import canonical_json_sha256, compute_package_source_sha256


def rm15_program():
    """Unrecycled five-wire PPR program in the positive odd-parity T convention."""
    fac = get_catalog_factory("bk_15_1_3_t")
    algebra = fac.verify_algebraic_structure()
    return dict(
        schema="ftqc-magic-ppr-1",
        factory=fac.name,
        prepare_plus=list(range(fac.k + fac.r)),
        rotations=[
            dict(
                column=i,
                parity=fac.column_mask(i),
                phase_on_odd_mod8=1,
                exponent_i_pi_P_over_8=-1,
            )
            for i in range(fac.n)
        ],
        corrections={k: v for k, v in algebra.items() if k.endswith("_corrections")},
        postselect_plus_X=list(range(fac.k, fac.k + fac.r)),
        output_slot=0,
    )


def integer(value, limit):
    return type(value) is int and 0 <= value < limit


def evaluate(batches, independent_attempts):
    """Input laws must have been verified by their native producer before this call."""
    fac = get_catalog_factory("bk_15_1_3_t")
    if not batches or (len(batches) > 1 and independent_attempts is not True):
        raise ValueError(
            "distinct source attempts require an independent-noise premise"
        )
    ids = set()
    allocated = set()
    sources = []
    distribution = {0: Fraction(1)}
    accepted = Fraction(1)
    for src, ports in batches:
        if (
            src["schema"] != "ftqc-injection-batch-1"
            or src["evidence"] != "exact-T-Z"
            or src["target"] != "T"
            or src["modulus"] != 8
            or src["phases"] != [0, 1]
            or src["attempt"] in ids
        ):
            raise ValueError("expected a distinct exact T source")
        ids.add(src["attempt"])
        members = set()
        mapping = []
        for member, column in ports:
            if (
                not integer(member, len(src["slots"]))
                or src["slots"][member] != member
                or not integer(column, 15)
                or member in members
                or column in allocated
            ):
                raise ValueError("reused or invalid member/rotation column")
            members.add(member)
            allocated.add(column)
            mapping.append((member, column))
        if not mapping:
            raise ValueError("unused attempt")
        contributions = {}
        for key, value in src["accepted_mass"].items():
            mask = int(key)
            p = Fraction(value)
            if p < 0 or not integer(mask, 1 << len(src["slots"])):
                raise ValueError("invalid source law")
            signature = 0
            for member, column in mapping:
                if mask >> member & 1:
                    signature ^= fac.column_syndrome(column) | (
                        fac.column_label(column) << 4
                    )
            contributions[signature] = contributions.get(signature, 0) + p
        a = sum(contributions.values(), Fraction(0))
        if a != Fraction(src["acceptance"]) or not 0 <= a <= 1:
            raise ValueError("source normalization")
        accepted *= a
        nxt = {}
        for (i, p), (j, q) in product(distribution.items(), contributions.items()):
            nxt[i ^ j] = nxt.get(i ^ j, 0) + p * q
        distribution = nxt
        sources.append(
            dict(
                attempt=src["attempt"],
                block=src["block"],
                source_sha256=canonical_json_sha256(src),
                mapping=ports,
                unused=sorted(set(src["slots"]) - members),
            )
        )
    if allocated != set(range(15)):
        raise ValueError("15 rotation columns required")
    clean = distribution.get(0, Fraction(0))
    wrong = distribution.get(16, Fraction(0))
    final = clean + wrong
    return dict(
        schema="ftqc-magic-insitu-rm15-1",
        source_sha256=compute_package_source_sha256(),
        catalog=fac.name,
        state="T",
        phases=[0, 1],
        modulus=8,
        inputs=sources,
        source_accept=str(accepted),
        final_accept=str(final),
        accepted_wrong=str(wrong),
        factory_accept_given_delivery=str(final / accepted) if accepted else None,
        conditional_output_Z=str(wrong / final) if final else None,
        accepted_mass={"0": str(clean), "1": str(wrong)},
        independent_attempts=independent_attempts,
        scope="one synchronized group of attempts; postselection; ideal factory Cliffords",
    )


def output_delivery(spec, proof):
    from ftqc_compose.__main__ import receipt as replay

    rebuilt = replay(spec)
    if canonical_json_sha256(rebuilt) != canonical_json_sha256(proof):
        raise ValueError("in-situ output receipt differs from replay")
    f = rebuilt["result"]["factory"]
    if f["conditional_output_Z"] is None:
        raise ValueError("no accepted output")
    delta = Fraction(f["conditional_output_Z"])
    return dict(
        schema="ftqc-magic-delivery-1",
        producer=dict(
            spec_sha256=canonical_json_sha256(spec),
            receipt_sha256=canonical_json_sha256(proof),
            source_sha256=f["source_sha256"],
        ),
        source_attempts=[row["attempt"] for row in f["inputs"]],
        native_output_count=1,
        states_per_batch=1,
        resource_wires=1,
        wire_order="little-endian",
        correlation="INDEPENDENT_SINGLE",
        acceptance_rule="all source heralds and four positive RM checks",
        conditional_infidelity_upper=str(delta),
        accepted_bad_mass_upper=f["accepted_wrong"],
        distance=None,
        leading_count=None,
        accepted_fault_witness=None,
        active_qubits=None,
        expected_wire_ticks_per_batch_upper=None,
        expected_wire_ticks_per_state_upper=None,
        expected_wire_ticks_per_output_upper=None,
        producer_clifford_corrections={"output_s_corrections": {"0": 1}},
        state="T",
        groups=[[0]],
        phases=[0, 1],
        modulus=8,
        frame="canonical-positive",
        law="exact-accepted-Z",
        conditional_z_mass={"0": str(1 - delta), "1": str(delta)},
        accepted_z_mass=f["accepted_mass"],
        acceptance_exact=True,
        acceptance_lower=f["final_accept"],
        conditional_trace_distance_upper=str(delta),
        premises=[f["scope"], "ideal logical input transfers and raw output interface"],
    )
