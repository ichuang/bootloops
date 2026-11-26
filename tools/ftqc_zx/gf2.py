"""Exact bit-packed F_2 and symplectic Pauli linear algebra (pure stdlib, zero floats).

Conventions for N-element Pauli vectors (edges or circuit wires):
- Element i uses two bits in a Python int:
    bit (2*i)     : z_i in {0, 1} (green / Z component)
    bit (2*i + 1) : x_i in {0, 1} (red   / X component)
- (z_i, x_i) = (0,0) -> I, (1,0) -> Z, (0,1) -> X, (1,1) -> Y.
"""

from __future__ import annotations

from functools import lru_cache


@lru_cache(maxsize=256)
def even_bit_mask(n_items: int) -> int:
    """Mask with bit 2*i set for 0 <= i < n_items."""
    mask = 0
    for i in range(n_items):
        mask |= 1 << (2 * i)
    return mask


def pack_pauli_char(ch: str) -> tuple[int, int]:
    """Return (z_bit, x_bit) for 'I', 'Z', 'X', or 'Y'."""
    if ch == "I":
        return (0, 0)
    if ch == "Z":
        return (1, 0)
    if ch == "X":
        return (0, 1)
    if ch == "Y":
        return (1, 1)
    raise ValueError(f"Unknown Pauli character: {ch!r}")


def unpack_pauli_char(z_bit: int, x_bit: int) -> str:
    """Return 'I', 'Z', 'X', or 'Y' from (z_bit, x_bit)."""
    if z_bit == 0 and x_bit == 0:
        return "I"
    if z_bit == 1 and x_bit == 0:
        return "Z"
    if z_bit == 0 and x_bit == 1:
        return "X"
    return "Y"


def pauli_string_to_symplectic(pauli_str: str) -> int:
    """Pack an N-character Pauli string ('IZXY...') into a 2N-bit int."""
    vec = 0
    for i, ch in enumerate(pauli_str):
        zb, xb = pack_pauli_char(ch)
        if zb:
            vec |= 1 << (2 * i)
        if xb:
            vec |= 1 << (2 * i + 1)
    return vec


def symplectic_to_pauli_string(vec: int, n_items: int) -> str:
    """Unpack a 2N-bit symplectic int into an N-character Pauli string."""
    chars: list[str] = []
    for i in range(n_items):
        zb = (vec >> (2 * i)) & 1
        xb = (vec >> (2 * i + 1)) & 1
        chars.append(unpack_pauli_char(zb, xb))
    return "".join(chars)


def swap_zx_bits(vec: int, n_items: int) -> int:
    """Swap z_i <-> x_i on each of the n_items 2-bit slots."""
    emask = even_bit_mask(n_items)
    return ((vec & emask) << 1) | ((vec >> 1) & emask)


def symplectic_inner_product(u: int, v: int, n_items: int) -> int:
    """Return 1 iff Pauli vectors u and v anticommute, else 0."""
    return (u & swap_zx_bits(v, n_items)).bit_count() & 1


def pauli_product_phase_mod4(u: int, v: int, n_items: int) -> int:
    """Return k in {0, 1, 2, 3} such that P(u) * P(v) = i^k * P(u ^ v) on n_items qubits.

    Uses standard single-qubit identities:
      X*Y = +i Z, Y*Z = +i X, Z*X = +i Y,
      Y*X = -i Z, Z*Y = -i X, X*Z = -i Y.
    """
    acc = 0
    for i in range(n_items):
        zu = (u >> (2 * i)) & 1
        xu = (u >> (2 * i + 1)) & 1
        zv = (v >> (2 * i)) & 1
        xv = (v >> (2 * i + 1)) & 1
        if (zu == 0 and xu == 0) or (zv == 0 and xv == 0) or (zu == zv and xu == xv):
            continue
        # Cyclic order: X=(0,1) -> Y=(1,1) -> Z=(1,0) -> X=(0,1) gives +i (+1 mod 4)
        if (
            (zu == 0 and xu == 1 and zv == 1 and xv == 1)
            or (zu == 1 and xu == 1 and zv == 1 and xv == 0)
            or (zu == 1 and xu == 0 and zv == 0 and xv == 1)
        ):
            acc += 1
        else:
            acc -= 1
    return acc & 3


def pauli_weight(vec: int, n_items: int) -> int:
    """Return number of items i where (z_i, x_i) != (0, 0)."""
    emask = even_bit_mask(n_items)
    active = (vec | (vec >> 1)) & emask
    return active.bit_count()


def rref_gf2(
    rows: list[int],
    ncols: int,
    track_transform: bool = False,
) -> tuple[list[int], list[int], list[int]]:
    """Compute canonical minimum-column RREF over F_2.

    Args:
      rows: list of m bit-packed rows (column j is bit j, 0 <= j < ncols).
      ncols: number of columns.
      track_transform: if True, maintain m-bit combination witnesses per row.

    Returns:
      (nonzero_rref_rows, pivot_cols, transform_rows)
      where if track_transform is True, transform_rows has length len(rows):
      the first r entries correspond to nonzero_rref_rows, and the remaining
      entries are null-combination witnesses (combinations of original rows = 0).
    """
    work = list(rows)
    m = len(work)
    transforms = [1 << i for i in range(m)] if track_transform else []
    pivot_cols: list[int] = []
    r = 0
    for col in range(ncols):
        bit = 1 << col
        piv = -1
        for i in range(r, m):
            if work[i] & bit:
                piv = i
                break
        if piv == -1:
            continue
        if piv != r:
            work[r], work[piv] = work[piv], work[r]
            if track_transform:
                transforms[r], transforms[piv] = transforms[piv], transforms[r]
        piv_row = work[r]
        piv_tr = transforms[r] if track_transform else 0
        for i in range(m):
            if i != r and (work[i] & bit):
                work[i] ^= piv_row
                if track_transform:
                    transforms[i] ^= piv_tr
        pivot_cols.append(col)
        r += 1
        if r == m:
            break
    return work[:r], pivot_cols, transforms


def nullspace_gf2(rows: list[int], ncols: int) -> list[int]:
    """Return a canonical basis for {x in F_2^ncols : A x = 0 mod 2}."""
    rref_rows, pivot_cols, _ = rref_gf2(rows, ncols, track_transform=False)
    piv_set = set(pivot_cols)
    basis: list[int] = []
    for fcol in range(ncols):
        if fcol in piv_set:
            continue
        vec = 1 << fcol
        fbit = 1 << fcol
        for idx, pcol in enumerate(pivot_cols):
            if rref_rows[idx] & fbit:
                vec |= 1 << pcol
        basis.append(vec)
    return basis


def solve_gf2(rows: list[int], rhs: list[int], ncols: int) -> tuple[bool, int, list[int]]:
    """Solve A x = rhs (mod 2).

    Returns:
      (is_consistent, x_particular, nullspace_basis).
    """
    if len(rows) != len(rhs):
        raise ValueError("rows and rhs must have identical length")
    rhs_col = ncols
    aug_rows = [row | ((b & 1) << rhs_col) for row, b in zip(rows, rhs)]
    rref_rows, pivot_cols, _ = rref_gf2(aug_rows, ncols + 1, track_transform=False)
    if pivot_cols and pivot_cols[-1] == rhs_col:
        return False, 0, []
    rhs_bit = 1 << rhs_col
    x_part = 0
    clean_rref: list[int] = []
    for row, pcol in zip(rref_rows, pivot_cols):
        if row & rhs_bit:
            x_part |= 1 << pcol
        clean_rref.append(row & ~rhs_bit)
    piv_set = set(pivot_cols)
    basis: list[int] = []
    for fcol in range(ncols):
        if fcol in piv_set:
            continue
        vec = 1 << fcol
        fbit = 1 << fcol
        for idx, pcol in enumerate(pivot_cols):
            if clean_rref[idx] & fbit:
                vec |= 1 << pcol
        basis.append(vec)
    return True, x_part, basis


def canonicalize_signed_stabilizers(
    generators: list[tuple[int, int]],
    n_qubits: int,
) -> tuple[bool, list[tuple[int, int]]]:
    """Reduce commuting signed Pauli generators [(symp_2n, sign_bit), ...] to canonical RREF.

    Here sign_bit in {0, 1} represents eigenvalue (-1)^sign_bit in {+1, -1}.
    Row multiplication uses exact multi-qubit Pauli phase tracking:
      ( (-1)^s1 P1 ) * ( (-1)^s2 P2 ) = (-1)^(s1 ^ s2 ^ (k // 2)) P(u1 ^ u2)
    where P1 * P2 = i^k P(u1 ^ u2). If two generators fail to commute (k odd),
    raises ValueError. If elimination produces -I (symp == 0 and sign == 1),
    returns (False, []) indicating the zero state/operator!

    Returns:
      (is_consistent_nonzero, canonical_rref_generators).
    """
    work = [(u, s & 1) for u, s in generators]
    m = len(work)
    ncols = 2 * n_qubits

    def mul_rows(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
        ua, sa = a
        ub, sb = b
        k = pauli_product_phase_mod4(ua, ub, n_qubits)
        if k & 1:
            raise ValueError("Stabilizer generators do not commute")
        new_sign = sa ^ sb ^ ((k >> 1) & 1)
        return (ua ^ ub, new_sign)

    r = 0
    for col in range(ncols):
        bit = 1 << col
        piv = -1
        for i in range(r, m):
            if work[i][0] & bit:
                piv = i
                break
        if piv == -1:
            continue
        if piv != r:
            work[r], work[piv] = work[piv], work[r]
        piv_item = work[r]
        for i in range(m):
            if i != r and (work[i][0] & bit):
                work[i] = mul_rows(work[i], piv_item)
        r += 1
        if r == m:
            break

    for i in range(r, m):
        if work[i][0] == 0 and work[i][1] == 1:
            return False, []
    return True, work[:r]
