"""Addressable paired logical bases, signed Pauli embeddings, and basis conversions.

Schema `ftqc_codes-basis-1` binds ordered logical slot names, signed physical
X and Z Pauli representatives, coordinate ordering, `code_sha256`, and
`basis_sha256`.  Supports generic binary overlap elimination, external candidate
basis verification, exact mod-4 signed mixed-Pauli (Y) embeddings, verified
symplectic basis changes with stabilizer dressings, and consumer export for
`ftqc_gadgets.surgery` and `ftqc_extractors`.
"""
import copy
import hashlib
import json

from . import certify, codes, gf2

SCHEMA = "ftqc_codes-basis-1"


def _canonical_sha256(payload):
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def pauli_from_xz(n, x, z, sign=0):
    """Format `(x, z, sign)` as signed Pauli string `±P_0...P_{n-1}` (qubit 0 first)."""
    if sign not in (0, 1):
        raise codes.SpecError(f"Pauli sign bit must be 0 or 1, got {sign!r}")
    if x < 0 or z < 0 or (x >> n) or (z >> n):
        raise codes.SpecError(f"Pauli bitmasks out of range for n={n}")
    chars = ["IXZY"[((x >> q) & 1) | (((z >> q) & 1) << 1)] for q in range(n)]
    return ("-" if sign else "+") + "".join(chars)


def parse_pauli(word, n):
    """Parse `±P_0...P_{n-1}` into `(x_mask, z_mask, sign_bit)`."""
    if (
        not isinstance(word, str)
        or len(word) != n + 1
        or word[0] not in "+-"
        or set(word[1:]) - {"I", "X", "Y", "Z"}
    ):
        raise codes.SpecError(f"expected signed Pauli string of length {n + 1} (qubit 0 first)")
    x = 0
    z = 0
    for q, ch in enumerate(word[1:]):
        if ch in ("X", "Y"):
            x |= 1 << q
        if ch in ("Z", "Y"):
            z |= 1 << q
    return x, z, (1 if word[0] == "-" else 0)


def commute_pauli(p1, p2):
    """Return True iff `(x1, z1, ...)` and `(x2, z2, ...)` commute."""
    return (gf2.dot(p1[0], p2[1]) ^ gf2.dot(p1[1], p2[0])) == 0


def mul_pauli(p1, p2):
    """Multiply two Pauli tuples `(x, z, phase_mod4_or_sign_bit)`.

    Here `p = (x, z, s)` represents `(-1)^s P(x, z)` when `p` is a 3-tuple; for
    chained mod-4 accumulation see `mul_pauli_mod4`.
    Returns `(x_out, z_out, phase_mod4)` where `P1 * P2 = i^{phase_mod4} P(x_out, z_out, 0)`.
    """
    x1, z1, s1 = p1
    x2, z2, s2 = p2
    xo = x1 ^ x2
    zo = z1 ^ z2
    phase = (
        2 * (s1 + s2 + gf2.dot(z1, x2))
        + (x1 & z1).bit_count()
        + (x2 & z2).bit_count()
        - (xo & zo).bit_count()
    ) % 4
    return xo, zo, phase


def _mul_chain_signed(ops, sign_init=0, i_power_init=0):
    """Multiply a sequence of signed Hermitian Paulis `(x, z, s)` with initial `(-1)^sign_init * i^i_power_init`.

    Requires the resulting operator to be Hermitian (`final_phase % 2 == 0`),
    returning `(x, z, final_sign_bit)`.
    """
    px, pz = 0, 0
    phase = (2 * sign_init + i_power_init) % 4
    for x, z, s in ops:
        nx = px ^ x
        nz = pz ^ z
        phase = (
            phase
            + 2 * (s + gf2.dot(pz, x))
            + (px & pz).bit_count()
            + (x & z).bit_count()
            - (nx & nz).bit_count()
        ) % 4
        px, pz = nx, nz
    if phase % 2 != 0:
        raise certify.Refusal("product of signed Paulis is anti-Hermitian (odd power of i)")
    return px, pz, phase // 2


class _SignedStabilizerGroup:
    """Exact signed stabilizer group over GF(2)^{2n} detecting `-I` inconsistencies."""

    def __init__(self, n, gens=()):
        self.n = n
        self.rows = {}  # highest symplectic bit -> (x, z, s)
        for g in gens:
            self.add(g)

    def reduce(self, g):
        x, z, s = g
        while x or z:
            piv = (x | (z << self.n)).bit_length() - 1
            if piv not in self.rows:
                break
            rx, rz, rs = self.rows[piv]
            x, z, ph = mul_pauli((x, z, s), (rx, rz, rs))
            if ph % 2 != 0:
                raise certify.Refusal("non-commuting element reduced in stabilizer group")
            s = ph // 2
        return x, z, s

    def add(self, g):
        x, z, s = self.reduce(g)
        if x or z:
            piv = (x | (z << self.n)).bit_length() - 1
            self.rows[piv] = (x, z, s)
            return True
        if s:
            raise certify.Refusal("inconsistent signed stabilizer group (-I generated)")
        return False

    def contains_exact(self, g):
        x, z, s = self.reduce(g)
        return x == 0 and z == 0 and s == 0


def _code_signed_checks(code):
    if not code.commutes():
        raise certify.Refusal("H_X and H_Z do not commute")
    if hasattr(code, "signed_checks"):
        return list(code.signed_checks)
    return [(r, 0, 0) for r in code.hx] + [(0, r, 0) for r in code.hz]


def pair_css_representatives(code, raw_lx, raw_lz):
    """Pair X representatives in `ker(H_Z)/row(H_X)` and Z in `ker(H_X)/row(H_Z)` so `L_X L_Z^T = I_k`."""
    ex, ez = gf2.Echelon(code.hx), gf2.Echelon(code.hz)
    k = code.n - ex.rank - ez.rank
    if len(raw_lx) != k or len(raw_lz) != k:
        raise certify.Refusal(f"expected k={k} X and Z representatives, got ({len(raw_lx)}, {len(raw_lz)})")
    lx = list(raw_lx)
    lz = list(raw_lz)
    for i in range(k):
        piv = next((j for j in range(i, k) if gf2.dot(lx[i], lz[j]) == 1), None)
        if piv is None:
            raise certify.Refusal(f"degenerate X/Z overlap matrix at row {i}")
        if piv != i:
            lz[i], lz[piv] = lz[piv], lz[i]
        for j in range(k):
            if j != i and gf2.dot(lx[i], lz[j]) == 1:
                lz[j] ^= lz[i]
    return lx, lz


def derive_paired_basis(code, *, slot_names=None, coordinate_order="0..n-1", provenance=None):
    """Derive a canonical verified `ftqc_codes-basis-1` receipt by binary overlap elimination."""
    if not code.commutes():
        raise certify.Refusal("H_X and H_Z do not commute")
    ex = gf2.Echelon(code.hx)
    ez = gf2.Echelon(code.hz)
    k = code.n - ex.rank - ez.rank
    raw_lx = []
    ex_work = gf2.Echelon(code.hx)
    for v in gf2.nullspace(code.hz, code.n):
        if ex_work.add(v):
            raw_lx.append(v)
    raw_lz = []
    ez_work = gf2.Echelon(code.hz)
    for v in gf2.nullspace(code.hx, code.n):
        if ez_work.add(v):
            raw_lz.append(v)
    lx, lz = pair_css_representatives(code, raw_lx, raw_lz)
    names = list(slot_names) if slot_names is not None else [f"q{i}" for i in range(k)]
    slots = [
        {
            "name": names[i],
            "x": pauli_from_xz(code.n, lx[i], 0, 0),
            "z": pauli_from_xz(code.n, 0, lz[i], 0),
        }
        for i in range(k)
    ]
    return build_basis_receipt(
        code,
        slots,
        source="binary_elimination",
        coordinate_order=coordinate_order,
        provenance=provenance,
    )


def build_basis_receipt(
    code,
    slots_or_candidate=None,
    *,
    source="supplied",
    coordinate_order="0..n-1",
    provenance=None,
):
    """Validate `slots_or_candidate` against `code` and emit an `ftqc_codes-basis-1` receipt."""
    if not code.commutes():
        raise certify.Refusal("H_X and H_Z do not commute")
    if slots_or_candidate is None:
        cand = code.meta.get("candidate_logicals")
        if cand is None:
            raise codes.SpecError("no slots provided and code carries no candidate_logicals")
        slots_or_candidate = cand

    if isinstance(slots_or_candidate, dict) and "lx_supports" in slots_or_candidate:
        names = slots_or_candidate["slot_names"]
        sx = slots_or_candidate["signs_x"]
        sz = slots_or_candidate["signs_z"]
        lxs = [gf2.from_supp(s) for s in slots_or_candidate["lx_supports"]]
        lzs = [gf2.from_supp(s) for s in slots_or_candidate["lz_supports"]]
        slots = [
            {
                "name": names[i],
                "x": pauli_from_xz(code.n, lxs[i], 0, sx[i]),
                "z": pauli_from_xz(code.n, 0, lzs[i], sz[i]),
            }
            for i in range(len(names))
        ]
    elif isinstance(slots_or_candidate, list):
        slots = copy.deepcopy(slots_or_candidate)
    else:
        raise codes.SpecError("slots_or_candidate must be a slot list or candidate_logicals dict")

    rec = {
        "schema": SCHEMA,
        "code_name": code.name,
        "code_sha256": code.sha256(),
        "n": code.n,
        "k": len(slots),
        "coordinate_order": str(coordinate_order),
        "source": str(source),
        "provenance": copy.deepcopy(provenance or {}),
        "slots": slots,
    }
    rec["basis_sha256"] = _canonical_sha256(
        {
            "code_sha256": rec["code_sha256"],
            "n": rec["n"],
            "coordinate_order": rec["coordinate_order"],
            "slots": rec["slots"],
        }
    )
    verify_basis_receipt(code, rec)
    return rec


def verify_basis_receipt(code, basis_receipt):
    """Verify membership, rank, independence mod stabilizers, unique slot names, and symplectic pairing."""
    if isinstance(code, dict) and "code" in code:
        code = codes.CSSCode.from_json(code["code"])
    if not isinstance(basis_receipt, dict) or basis_receipt.get("schema") != SCHEMA:
        raise certify.Refusal(f"expected basis schema {SCHEMA!r}")
    if not code.commutes():
        raise certify.Refusal("code checks do not commute")
    if basis_receipt.get("code_sha256") != code.sha256():
        raise certify.Refusal(
            f"basis code_sha256 {basis_receipt.get('code_sha256')} != code {code.sha256()}"
        )
    if basis_receipt.get("n") != code.n:
        raise certify.Refusal("basis n does not match code.n")

    checks = _code_signed_checks(code)
    stab_symp = gf2.Echelon([x | (z << code.n) for x, z, _ in checks])
    k = code.n - stab_symp.rank
    slots = basis_receipt.get("slots")
    if not isinstance(slots, list) or basis_receipt.get("k") != k or len(slots) != k:
        raise certify.Refusal(f"basis has {len(slots) if isinstance(slots, list) else None} slots; code has k={k}")

    want_bsha = _canonical_sha256(
        {
            "code_sha256": basis_receipt["code_sha256"],
            "n": basis_receipt["n"],
            "coordinate_order": basis_receipt["coordinate_order"],
            "slots": slots,
        }
    )
    if basis_receipt.get("basis_sha256") != want_bsha:
        raise certify.Refusal("basis_sha256 does not match canonical slot payload")

    # Also verify checks themselves generate a consistent signed group
    _SignedStabilizerGroup(code.n, checks)

    seen_names = set()
    parsed_x = []
    parsed_z = []
    for i, item in enumerate(slots):
        if not isinstance(item, dict) or set(item) != {"name", "x", "z"}:
            raise codes.SpecError(f"slot {i} must have exact keys 'name', 'x', 'z'")
        nm = item["name"]
        if not isinstance(nm, str) or not nm:
            raise codes.SpecError(f"slot {i} name must be a non-empty string")
        if nm in seen_names:
            raise certify.Refusal(f"duplicated logical slot name {nm!r}")
        seen_names.add(nm)
        try:
            px = parse_pauli(item["x"], code.n)
            pz = parse_pauli(item["z"], code.n)
        except codes.SpecError as exc:
            raise certify.Refusal(f"slot {nm!r} invalid Pauli: {exc}") from None

        for axis_lbl, op in (("X", px), ("Z", pz)):
            if any(not commute_pauli(op, chk) for chk in checks):
                raise certify.Refusal(f"slot {nm!r} {axis_lbl} anticommutes with a code stabilizer")
        parsed_x.append(px)
        parsed_z.append(pz)

    # Symplectic pairing: X_i commutes with X_j, Z_i commutes with Z_j, X_i anticommutes iff j == i
    for i in range(k):
        for j in range(k):
            if not commute_pauli(parsed_x[i], parsed_x[j]):
                raise certify.Refusal(f"X[{i}] and X[{j}] anticommute")
            if not commute_pauli(parsed_z[i], parsed_z[j]):
                raise certify.Refusal(f"Z[{i}] and Z[{j}] anticommute")
            want_comm = i != j
            if commute_pauli(parsed_x[i], parsed_z[j]) != want_comm:
                raise certify.Refusal(
                    f"pairing failure between X[{i}] ({slots[i]['name']}) and Z[{j}] ({slots[j]['name']})"
                )

    # Independence modulo stabilizers over GF(2)^{2n}
    ext = gf2.Echelon([x | (z << code.n) for x, z, _ in checks])
    for i in range(k):
        if not ext.add(parsed_x[i][0] | (parsed_x[i][1] << code.n)):
            raise certify.Refusal(f"slot {slots[i]['name']!r} X is dependent modulo stabilizers")
        if not ext.add(parsed_z[i][0] | (parsed_z[i][1] << code.n)):
            raise certify.Refusal(f"slot {slots[i]['name']!r} Z is dependent modulo stabilizers")
    if ext.rank != stab_symp.rank + 2 * k:
        raise certify.Refusal("logical generators lack full rank 2k modulo stabilizers")

    return [
        f"basis_sha256 {want_bsha[:16]} verified on code {code.sha256()[:16]}",
        f"{k} named symplectic logical slots verified",
    ]


def slot_index(basis_receipt, name_or_idx):
    """Resolve a slot name or integer index into `0..k-1`."""
    slots = basis_receipt["slots"]
    if isinstance(name_or_idx, int) and not isinstance(name_or_idx, bool):
        if not (0 <= name_or_idx < len(slots)):
            raise certify.Refusal(f"slot index {name_or_idx} out of range 0..{len(slots) - 1}")
        return name_or_idx
    for i, s in enumerate(slots):
        if s["name"] == name_or_idx:
            return i
    raise certify.Refusal(f"unknown logical slot {name_or_idx!r}")


def embed_logical_pauli(basis_receipt, x_mask, z_mask, sign=0):
    """Embed a logical Pauli `(-1)^sign * prod_j L_j` into a physical signed Pauli string.

    For each slot `j` where both `x_mask` and `z_mask` have bit `j` set, uses
    `Y_j = +i X_j Z_j` with exact mod-4 phase tracking across all physical qubits.
    """
    n = basis_receipt["n"]
    slots = basis_receipt["slots"]
    k = len(slots)
    if x_mask < 0 or z_mask < 0 or (x_mask >> k) or (z_mask >> k) or sign not in (0, 1):
        raise codes.SpecError("logical Pauli mask or sign out of range")
    y_count = (x_mask & z_mask).bit_count()
    ops = []
    for j in gf2.supp(x_mask):
        ops.append(parse_pauli(slots[j]["x"], n))
    for j in gf2.supp(z_mask):
        ops.append(parse_pauli(slots[j]["z"], n))
    # Wait: in mul_j (i^{x_j & z_j} X_j^{x_j} Z_j^{z_j}), slots commute pairwise for i != j!
    # Ordering all X_j first and Z_j second has X_{j1} passing Z_{j2} only for j1 == j2 (which is already X_j before Z_j!).
    px, pz, ps = _mul_chain_signed(ops, sign_init=sign, i_power_init=y_count)
    return pauli_from_xz(n, px, pz, ps)


def embed_named_pauli(basis_receipt, terms, *, sign=0):
    """Embed `terms` like `[('q0', 'Z'), ('q1', 'Y')]` with overall `sign` (0 for +, 1 for -)."""
    seen = set()
    x_mask = 0
    z_mask = 0
    for slot_ref, axis in terms:
        idx = slot_index(basis_receipt, slot_ref)
        if idx in seen:
            raise certify.Refusal(f"duplicated slot {slot_ref!r} in logical target")
        seen.add(idx)
        if axis == "I":
            continue
        if axis not in ("X", "Y", "Z"):
            raise codes.SpecError(f"invalid Pauli axis {axis!r}")
        if axis in ("X", "Y"):
            x_mask |= 1 << idx
        if axis in ("Z", "Y"):
            z_mask |= 1 << idx
    return embed_logical_pauli(basis_receipt, x_mask, z_mask, sign=sign)


def verify_logical_target(code, basis_receipt, terms, claimed_physical_pauli, *, sign=0):
    """Verify `claimed_physical_pauli` equals the embedded logical operator modulo `+1` stabilizers.

    Catches wrong signs on mixed-Pauli (`Y`) targets, swapped slots, or stale bases.
    """
    verify_basis_receipt(code, basis_receipt)
    want_str = embed_named_pauli(basis_receipt, terms, sign=sign)
    want = parse_pauli(want_str, code.n)
    got = parse_pauli(claimed_physical_pauli, code.n)
    checks = _code_signed_checks(code)
    if any(not commute_pauli(got, c) for c in checks):
        raise certify.Refusal("claimed physical target anticommutes with a stabilizer")
    # Multiply want * got: since both are Hermitian and commute with stabilizers, if they represent
    # the same signed logical observable then want * got must lie in the +1 stabilizer group!
    dx, dz, ph = mul_pauli(want, got)
    if ph % 2 != 0:
        raise certify.Refusal("claimed physical target anticommutes with requested logical target")
    diff = (dx, dz, ph // 2)
    sgrp = _SignedStabilizerGroup(code.n, checks)
    rx, rz, rs = sgrp.reduce(diff)
    if rx or rz:
        raise certify.Refusal("claimed physical target implements a different logical operator")
    if rs != 0:
        raise certify.Refusal("claimed physical target has the wrong sign (-1 vs +1)")
    return want_str


def dress_basis(
    code,
    basis_receipt,
    *,
    permute_slots=None,
    new_slot_names=None,
    dress_x_hx=None,
    dress_x_hz=None,
    dress_z_hx=None,
    dress_z_hz=None,
    flip_signs_x=None,
    flip_signs_z=None,
):
    """Construct a dressed / relabeled `ftqc_codes-basis-1` receipt and an explicit transform proof."""
    verify_basis_receipt(code, basis_receipt)
    k = basis_receipt["k"]
    n = code.n
    perm = list(permute_slots) if permute_slots is not None else list(range(k))
    if len(perm) != k or set(perm) != set(range(k)):
        raise certify.Refusal("permute_slots must be a bijection on 0..k-1")

    dx_hx = dress_x_hx or [[] for _ in range(k)]
    dx_hz = dress_x_hz or [[] for _ in range(k)]
    dz_hx = dress_z_hx or [[] for _ in range(k)]
    dz_hz = dress_z_hz or [[] for _ in range(k)]
    fx = flip_signs_x or [0] * k
    fz = flip_signs_z or [0] * k

    old_slots = basis_receipt["slots"]
    names = list(new_slot_names) if new_slot_names is not None else [old_slots[i]["name"] for i in range(k)]

    new_slots = []
    symplectic_x = []
    symplectic_z = []
    for i in range(k):
        src = perm[i]
        base_x = parse_pauli(old_slots[src]["x"], n)
        ops_x = [base_x] + [(code.hx[r], 0, 0) for r in dx_hx[i]] + [(0, code.hz[r], 0) for r in dx_hz[i]]
        nx, nz, ns = _mul_chain_signed(ops_x, sign_init=fx[i], i_power_init=0)

        base_z = parse_pauli(old_slots[src]["z"], n)
        ops_z = [base_z] + [(code.hx[r], 0, 0) for r in dz_hx[i]] + [(0, code.hz[r], 0) for r in dz_hz[i]]
        zx, zz, zs = _mul_chain_signed(ops_z, sign_init=fz[i], i_power_init=0)

        new_slots.append(
            {
                "name": names[i],
                "x": pauli_from_xz(n, nx, nz, ns),
                "z": pauli_from_xz(n, zx, zz, zs),
            }
        )
        symplectic_x.append({"x_mask": 1 << src, "z_mask": 0, "sign": int(fx[i])})
        symplectic_z.append({"x_mask": 0, "z_mask": 1 << src, "sign": int(fz[i])})

    new_rec = build_basis_receipt(
        code,
        new_slots,
        source="dressed_transform",
        coordinate_order=basis_receipt["coordinate_order"],
        provenance={"from_basis_sha256": basis_receipt["basis_sha256"]},
    )
    transform = {
        "code_sha256": code.sha256(),
        "old_basis_sha256": basis_receipt["basis_sha256"],
        "new_basis_sha256": new_rec["basis_sha256"],
        "logical_x": symplectic_x,
        "logical_z": symplectic_z,
        "dressings": {
            "x_hx": [list(v) for v in dx_hx],
            "x_hz": [list(v) for v in dx_hz],
            "z_hx": [list(v) for v in dz_hx],
            "z_hz": [list(v) for v in dz_hz],
        },
    }
    verify_basis_transform(code, basis_receipt, new_rec, transform)
    return new_rec, transform


def verify_basis_transform(code, old_basis, new_basis, transform):
    """Verify an explicit symplectic logical transformation + stabilizer dressing against physical rows."""
    verify_basis_receipt(code, old_basis)
    verify_basis_receipt(code, new_basis)
    if old_basis["code_sha256"] != new_basis["code_sha256"] or old_basis["code_sha256"] != code.sha256():
        raise certify.Refusal("basis conversion requires matching code_sha256")
    if (
        transform.get("code_sha256") != code.sha256()
        or transform.get("old_basis_sha256") != old_basis["basis_sha256"]
        or transform.get("new_basis_sha256") != new_basis["basis_sha256"]
    ):
        raise certify.Refusal("transform digest binding mismatch")

    k = old_basis["k"]
    n = code.n
    lx_specs = transform.get("logical_x")
    lz_specs = transform.get("logical_z")
    dress = transform.get("dressings") or {}
    if not isinstance(lx_specs, list) or not isinstance(lz_specs, list) or len(lx_specs) != k or len(lz_specs) != k:
        raise certify.Refusal("transform requires k logical_x and k logical_z entries")

    def symp_comm(a, b):
        return (gf2.dot(a["x_mask"], b["z_mask"]) ^ gf2.dot(a["z_mask"], b["x_mask"])) == 0

    for i in range(k):
        for j in range(k):
            if not symp_comm(lx_specs[i], lx_specs[j]) or not symp_comm(lz_specs[i], lz_specs[j]):
                raise certify.Refusal("declared logical transformation violates X/X or Z/Z commutation")
            if symp_comm(lx_specs[i], lz_specs[j]) != (i != j):
                raise certify.Refusal("declared logical transformation is not symplectic")

    for i in range(k):
        for axis, spec_list, hx_key, hz_key, slot_key in (
            ("X", lx_specs, "x_hx", "x_hz", "x"),
            ("Z", lz_specs, "z_hx", "z_hz", "z"),
        ):
            sp = spec_list[i]
            emb_str = embed_logical_pauli(old_basis, sp["x_mask"], sp["z_mask"], sign=sp.get("sign", 0))
            base_op = parse_pauli(emb_str, n)
            hx_ids = dress.get(hx_key, [[] for _ in range(k)])[i]
            hz_ids = dress.get(hz_key, [[] for _ in range(k)])[i]
            if any(not (0 <= r < len(code.hx)) for r in hx_ids) or any(
                not (0 <= r < len(code.hz)) for r in hz_ids
            ):
                raise certify.Refusal("dressing stabilizer row index out of range")
            ops = [base_op] + [(code.hx[r], 0, 0) for r in hx_ids] + [(0, code.hz[r], 0) for r in hz_ids]
            px, pz, ps = _mul_chain_signed(ops)
            rebuilt = pauli_from_xz(n, px, pz, ps)
            if rebuilt != new_basis["slots"][i][slot_key]:
                raise certify.Refusal(
                    f"slot {i} ({new_basis['slots'][i]['name']}) {axis} mismatch: "
                    f"transform gives {rebuilt}, new_basis has {new_basis['slots'][i][slot_key]}"
                )
    return True


def export_surgery_basis(basis_receipt, *, measured_slots=(), surviving_slots=None):
    """Export a consumer-ready basis dictionary bound to `code_sha256` and `basis_sha256`.

    `result["surgery_basis"]` has the exact `{"input": [...], "surviving": [...]}`
    shape expected by `ftqc_gadgets.surgery.Instrument`.
    """
    slots = basis_receipt["slots"]
    names = [s["name"] for s in slots]
    meas_names = []
    for m in measured_slots:
        nm = names[slot_index(basis_receipt, m)]
        if nm in meas_names:
            raise certify.Refusal(f"duplicated measured slot {nm!r}")
        meas_names.append(nm)

    if surviving_slots is None:
        surv = [nm for nm in names if nm not in meas_names]
    else:
        surv = []
        for s in surviving_slots:
            nm = names[slot_index(basis_receipt, s)]
            if nm in surv:
                raise certify.Refusal(f"duplicated surviving slot {nm!r}")
            if nm in meas_names:
                raise certify.Refusal(f"slot {nm!r} listed as both measured and surviving")
            surv.append(nm)

    return {
        "code_sha256": basis_receipt["code_sha256"],
        "basis_sha256": basis_receipt["basis_sha256"],
        "measured": meas_names,
        "surgery_basis": {
            "input": copy.deepcopy(slots),
            "surviving": surv,
        },
    }
