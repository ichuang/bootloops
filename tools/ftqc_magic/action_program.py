"""Transport the RM15 Pauli-product rotations into a verified logical frame."""

from ftqc_codes import action_domain, automorphisms, basis
from .insitu import rm15_program


def transport(spec):
    if (
        set(spec) != {"kind", "code", "action", "slots"}
        or spec["kind"] != "rm15-action"
    ):
        raise ValueError("expected an rm15-action specification")
    code = action_domain.load(spec["code"])
    action = automorphisms.verify(code, spec["action"])
    paired = action["basis"]
    program = rm15_program()
    slots = spec["slots"]
    if not isinstance(slots, list) or len(slots) != len(program["prepare_plus"]):
        raise ValueError("expected one ordered destination slot per factory wire")
    sites = [basis.slot_index(paired, name) for name in slots]
    if len(set(sites)) != len(sites):
        raise ValueError("factory slots must be distinct")

    def image(axis, mask):
        word = sum(1 << site for j, site in enumerate(sites) if mask >> j & 1)
        logical = (word, 0, 0) if axis == "X" else (0, word, 0)
        return basis.pauli_from_xz(
            paired["k"], *automorphisms.logical_image(action, logical)
        )

    out = program["output_slot"]
    return dict(
        status="proved",
        scope="conjugated RM15 observables and signed rotation axes",
        code_sha256=action["code_sha256"],
        basis_sha256=action["basis_sha256"],
        prepare_plus=[image("X", 1 << j) for j in program["prepare_plus"]],
        rotations=[
            dict(row, signed_axis=image("Z", row["parity"]))
            for row in program["rotations"]
        ],
        postselect_plus=[image("X", 1 << j) for j in program["postselect_plus_X"]],
        output_x=image("X", 1 << out),
        output_z=image("Z", 1 << out),
    )
