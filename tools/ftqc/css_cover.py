"""Finite binary CSS distance proofs by branching on an unsatisfied parity check.

The proof stream supplies a check to split on, or a packing of disjoint
checks to close a branch. The verifier checks these choices instead of running
the generator's search heuristic. See GROSS.md for the coverage argument.
"""

from base64 import b64decode, b64encode
import zlib

from qeccert.gf2 import Echelon

from .algebra import Limit, bits, echelon, reduce_row, transpose


def encode(stream):
    return b64encode(zlib.compress(stream, 9)).decode("ascii")


def decode(encoded, cap):
    if not isinstance(encoded, str) or len(encoded) > (cap + 1024) * 2:
        raise ValueError("invalid cover encoding")
    try:
        decoder = zlib.decompressobj()
        raw = decoder.decompress(b64decode(encoded, validate=True), cap + 1)
    except zlib.error as exc:
        raise ValueError("invalid compressed cover") from exc
    if len(raw) > cap or not decoder.eof or decoder.unused_data:
        raise ValueError("oversized, truncated or trailing cover")
    return raw


def cover(n, checks, stabilizers, limit, budget):
    """Exclude nontrivial zero-syndrome supports through limit or return a witness."""
    if not 0 <= limit < 126 or len(checks) >= 128:
        raise ValueError(
            "cover format requires fewer than 128 checks and bounded depth"
        )
    columns = transpose(checks, n)
    span = echelon(stabilizers)
    stream = bytearray()
    nodes = 0

    def visit(chosen, available, syndrome, left):
        nonlocal nodes
        nodes += 1
        if nodes > budget:
            raise Limit("CSS distance cover node budget")
        if not syndrome:
            return chosen if reduce_row(chosen, span) else 0
        if not left:
            return 0

        occupied = 0
        packing = []
        pivot = None
        fewest = n + 1
        for row in bits(syndrome):
            choices = checks[row] & available
            size = choices.bit_count()
            if size < fewest:
                pivot, fewest = row, size
                if not size:
                    break
            if not choices & occupied:
                packing.append(row)
                occupied |= choices
                if len(packing) > left:
                    stream.append(128 + len(packing))
                    stream.extend(packing)
                    return 0
        stream.append(pivot)
        for site in bits(checks[pivot] & available):
            bit = 1 << site
            available ^= bit
            witness = visit(chosen | bit, available, syndrome ^ columns[site], left - 1)
            if witness:
                return witness
        return 0

    for root in range(n) if limit else ():
        bit = 1 << root
        witness = visit(bit, (1 << n) - 2 * bit, columns[root], limit - 1)
        if witness:
            return witness, nodes, None
    return 0, nodes, encode(stream)


def verify_cover(n, checks, stabilizers, limit, declared_nodes, encoded, budget):
    """Replay the split/packing tree using the declared rows and GF(2) membership."""
    if type(declared_nodes) is not int or not 0 <= declared_nodes <= budget:
        raise ValueError("invalid cover node count")
    raw = decode(encoded, declared_nodes * (limit + 2))
    at = nodes = 0
    group = Echelon(stabilizers)
    columns = [
        sum(1 << i for i, h in enumerate(checks) if h & (1 << q)) for q in range(n)
    ]

    def take():
        nonlocal at
        if at >= len(raw):
            raise ValueError("incomplete cover tree")
        value = raw[at]
        at += 1
        return value

    def check(chosen, free, odd, remaining):
        nonlocal nodes
        nodes += 1
        if nodes > declared_nodes:
            raise ValueError("cover exceeds declared node count")
        if odd == 0:
            if chosen not in group:
                raise ValueError("cover contains a surviving logical")
            return
        if remaining == 0:
            return
        token = take()
        if token >= 128:
            count = token - 128
            if count != remaining + 1:
                raise ValueError("invalid packing size")
            used_sites = used_checks = 0
            for _ in range(count):
                i = take()
                if not (odd >> i & 1) or used_checks >> i & 1:
                    raise ValueError("packing requires distinct unsatisfied checks")
                sites = free & checks[i]
                if used_sites & sites:
                    raise ValueError("packing supports overlap")
                used_sites |= sites
                used_checks |= 1 << i
            return
        if not odd >> token & 1:
            raise ValueError("branch requires an unsatisfied check")
        candidates = checks[token] & free
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            free -= bit
            check(
                chosen | bit, free, odd ^ columns[bit.bit_length() - 1], remaining - 1
            )

    if limit:
        for q, syndrome in enumerate(columns):
            check(1 << q, ((1 << n) - 1) ^ ((2 << q) - 1), syndrome, limit - 1)
    if at != len(raw) or nodes != declared_nodes:
        raise ValueError("unused cover bytes or inconsistent node count")
