"""Joint record and Pauli-frame recovery for signed surgery instruments."""

from collections import Counter, defaultdict
from fractions import Fraction


def single_fault_table(groups):
    """Most likely joint response at first order for uniform categorical site faults.

    Zero syndrome preserves the ideal operation. Unknown syndromes use the
    identity frame. Ties choose the lowest (public mask, residual mask).
    """
    masses = defaultdict(Counter)
    for group in groups:
        weight = Fraction(1, len(group))
        for _, (det, public, residual) in group:
            masses[det][public, residual] += weight
    table = {
        det: min(law, key=lambda effect: (-law[effect], effect))
        for det, law in masses.items()
        if det
    }
    table[0] = (0, 0)
    return table


def decoded_damage(table, signature):
    det, bit, residual = signature
    flip, frame = table.get(det, (0, 0))
    return bit ^ flip, residual ^ frame
