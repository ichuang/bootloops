"""Map RM15 pivot parities to the translated Gross mono-layer Z instrument."""

from functools import lru_cache, reduce
from operator import xor
from ftqc_codes.gf2 import Echelon, nullspace, dot
from ftqc_gadgets.surgery import Instrument
from ftqc_magic.insitu import rm15_program
from .gross import shifted_z, gross_z, port_frame, cell
from .gross_cycle import cycle, X0, Z0, NWIRES
from .algebra import bits
from .block_map import word
from .block_reading import row_solver


def combine(rows, mask):
    return reduce(xor, (rows[i] for i in bits(mask)), 0)


def translate(v, dx, dy):
    return sum(
        1 << (72 * (q // 72) + cell(q % 72 // 6 + dx, q % 6 + dy)) for q in bits(v)
    )


@lru_cache(maxsize=1)
def native():
    """Single-round Gross merge, public Z parity, split and physical Z frame."""
    m, _, ports, incident = gross_z()
    n = m.n
    edge = list(range(n, X0))
    ticks = [
        [["RX", q] for q in edge],
        *[[list(op) for op in tick] for tick in cycle()],
        [["MX", q] for q in edge],
    ]
    paths, cycles = port_frame(ports, n, len(edge))
    width = NWIRES - X0
    for q, path in zip(bits(shifted_z()), paths, strict=True):
        for j in bits(path):
            ticks.append([["Z_IF", q, width + j]])
    spec = dict(
        n=n,
        checks=[word(n, x, 0) for x in m.hx] + [word(n, 0, z) for z in m.hz],
        targets=[word(n, 0, shifted_z())],
        ticks=ticks,
        readouts=[dict(records=list(range(Z0 - X0 + len(m.hz), width)), sign=0)],
        detectors=[[i] for i in range(Z0 - X0 + len(m.hz))]
        + [[width + i for i in bits(c)] for c in cycles],
        terminal_syndrome=False,
    )
    Instrument(spec)
    return spec


@lru_cache(maxsize=1)
def layout():
    m = gross_z()[0]
    translations = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1)]
    orbit_basis = [shifted_z(i, j) for i, j in translations]
    # A four-dimensional affine flat of the six-dimensional translation
    # orbit. On the five factory slots its points are Zout times a check
    # parity; the pivot Z makes each required parity directly measurable.
    zs = [combine(orbit_basis, c) for c in (2, 7, 9, 19, 35, 1)]
    group = Echelon(m.hz)
    if not all(group.add(v) for v in zs):
        raise ValueError("dependent factory layout")
    zs.extend(v for v in nullspace(m.hx, 144) if group.add(v))
    # Use one independent spare in the +Z frame to dress the first check.
    # Unmasking is legal only while this logical is initialized to zero.
    zs[1] ^= zs[6]
    xspace = nullspace(m.hz, 144)
    dual = row_solver([sum(dot(x, z) << j for j, z in enumerate(zs)) for x in xspace])
    xs = [combine(xspace, dual(1 << i)) for i in range(12)]
    if len(zs) != 12 or any(
        dot(x, z) != (i == j) for i, x in enumerate(xs) for j, z in enumerate(zs)
    ):
        raise ValueError("invalid paired basis")
    return m, xs, zs


def mapping():
    from .__main__ import digest

    m, xs, zs = layout()
    program = rm15_program()
    native_spec = native()
    stabs = Echelon(m.hz)
    orbit = {}
    for dx in range(12):
        for dy in range(6):
            physical = shifted_z(dx, dy)
            orbit.setdefault(stabs.reduce(physical), (dx, dy, physical))
    # Re-derive the entire reachable signed orbit in the factory's own basis.
    # Keep the established (dx,dy) representative and route for each class.
    from ftqc_codes import automorphisms as aut, basis
    from ftqc_codes.orbit_minweight import infer_default_automorphisms
    from .action_orbit import orbit as signed_orbit
    paired = basis.build_basis_receipt(m, [
        dict(name=f"q{j}", x=word(m.n, x, 0, 0), z=word(m.n, 0, z, 0))
        for j, (x, z) in enumerate(zip(xs, zs))
    ])
    actions = [aut.build(m, paired, [{"perm": p}]) for p in infer_default_automorphisms(m)]
    seed = (0, sum(dot(x, shifted_z()) << j for j, x in enumerate(xs)), 0)
    reached = signed_orbit(actions, seed, 72)
    if len(reached) != len(orbit):
        raise ValueError("signed logical orbit disagrees with physical translations")
    rotations = []
    for rotation in program["rotations"]:
        i, mask = rotation["column"], rotation["parity"]
        desired = combine(zs, mask) ^ zs[5]
        for extra in (0, 1 << 6):
            key = stabs.reduce(desired ^ combine(zs, extra))
            if key in orbit:
                break
        else:
            raise ValueError("unavailable Gross measurement")
        dx, dy, physical = orbit[key]
        if any(
            {translate(h, dx, dy) for h in checks} != set(checks)
            for checks in (m.hx, m.hz)
        ):
            raise ValueError("shift is not a BB automorphism")
        actual = sum(dot(x, physical) << j for j, x in enumerate(xs))
        image = seed
        for g, count in zip(actions, (dx, dy), strict=True):
            for _ in range(count):
                image = aut.logical_image(g, image)
        if image != (0, actual, 0) or image not in reached:
            raise ValueError("generator word disagrees with the installed Gross port")
        if actual != mask ^ (1 << 5) ^ extra:
            raise ValueError("native setting changes the signed port or spectators")
        rotations.append(
            dict(
                column=i,
                parity=mask,
                pivot=5,
                mask_zero_slots=list(bits(extra)),
                shift=[dx, dy],
                physical_z=physical,
                sign=0,
                exponent_i_pi_P_over_8=rotation["exponent_i_pi_P_over_8"],
                phase_on_odd_mod8=rotation["phase_on_odd_mod8"],
            )
        )
    return dict(
        schema="bb144-rm15-map-1",
        status="proved",
        memory="BB144",
        check_sha256=digest([m.hx, m.hz]),
        basis=[[x, z] for x, z in zip(xs, zs)],
        native_circuit_sha256=digest(native_spec),
        orbit_size=len(orbit),
        magic_program_sha256=digest(program),
        factory_plus_slots=program["prepare_plus"],
        pivot_slot=5,
        zero_slots=[6],
        spectators=list(range(7, 12)),
        output_slot=program["output_slot"],
        postselect_plus_X=program["postselect_plus_X"],
        clifford_corrections=program["corrections"],
        rotations=rotations,
        scope="one-round Gross Z Choi check transported by exact coordinate translations",
    )


def analyze(s):
    from .__main__ import digest

    m = mapping()
    if "candidate" in s and digest(s["candidate"]) != digest(m):
        raise ValueError("candidate differs from checked BB144 map")
    return m
