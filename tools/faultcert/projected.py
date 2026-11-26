"""Exclude low-weight detector-silent sums in a projection of fault columns."""

from itertools import combinations


def indices(v):
    while v:
        b = v & -v
        yield b.bit_length() - 1
        v -= b


class NodeLimit(Exception):
    pass


def exclude(columns, nd, limit, budget):
    """Return (proved, nodes, witness) with a relaxed XOR-column weight.

    The final two choices use an unrestricted syndrome table. A table hit
    gives a valid global sum after cancelling duplicate columns. A miss rules
    out every completion of this branch. A projected hit alone supplies no
    physical circuit witness.
    """
    cols = sorted(set(columns) - {0})
    mask = (1 << nd) - 1
    n = len(cols)
    checks = [0] * nd
    for i, col in enumerate(cols):
        for r in indices(col & mask):
            checks[r] |= 1 << i
    # For a fixed syndrome, one of two distinct payloads must differ from
    # any requested payload. Keeping two suffices for a nonzero-harm test.
    table = {}

    def add(value, ids):
        entries = table.setdefault(value & mask, [])
        if len(entries) < 2 and all(v != value for v, _ in entries):
            entries.append((value, ids))

    add(0, ())
    if limit >= 1:
        for i, c in enumerate(cols):
            add(c, (i,))
    if limit >= 2:
        for a, b in combinations(range(n), 2):
            add(cols[a] ^ cols[b], (a, b))
    nodes = 0

    def walk(value, free, left, prefix):
        nonlocal nodes
        nodes += 1
        if nodes > budget:
            raise NodeLimit
        odd = value & mask
        if not odd:
            if value:
                return prefix
            if prefix:
                return None
        if left <= 2:
            for v, tail in table.get(odd, ()):
                if v != value:
                    # All lookups take place at min(2,initial limit) remaining.
                    if len(tail) > left:
                        raise ValueError("invalid suffix budget")
                    return tuple(sorted(set(prefix) ^ set(tail)))
            return None
        choices = None
        used = packed = 0
        for row in indices(odd):
            opt = checks[row] & free
            if not opt:
                return None
            if choices is None or opt.bit_count() < choices.bit_count():
                choices = opt
            if not opt & used:
                used |= opt
                packed += 1
                if packed > left:
                    return None
        for i in indices(choices):
            free ^= 1 << i
            found = walk(value ^ cols[i], free, left - 1, prefix + (i,))
            if found is not None:
                return found
        return None

    try:
        if limit <= 2:
            hit = walk(0, 0, limit, ())
        else:
            hit = None
            for i, c in enumerate(cols):
                hit = walk(c, (1 << n) - (2 << i), limit - 1, (i,))
                if hit is not None:
                    break
        return hit is None, nodes, None if hit is None else [cols[i] for i in hit]
    except NodeLimit:
        return False, budget, None
