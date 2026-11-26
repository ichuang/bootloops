"""Enumerate projected supports and check time-resolved detector choices."""

from collections import Counter
from contextlib import contextmanager
from functools import reduce
from itertools import chain
from operator import or_, xor
from pathlib import Path
import subprocess
import tempfile


@contextmanager
def support_enumerator():
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / "supports"
        src = Path(__file__).with_name("projected_supports.cpp")
        cmd = ["c++", "-std=c++17", "-O3", "-march=native", str(src), "-o", str(exe)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode:
            raise ValueError(f"could not compile support enumerator: {p.stderr[-500:]}")
        yield exe


def enumerate_supports(exe, cols, nd, width, budget):
    if max(cols).bit_length() > 128:
        raise ValueError("projection exceeds two machine limbs")
    data = [f"{nd} {len(cols)}"] + [f"{v % 2**64} {v >> 64}" for v in cols]
    p = subprocess.run(
        [str(exe), str(width), str(budget)],
        input="\n".join(data) + "\n",
        capture_output=True,
        text=True,
    )
    if p.returncode == 3:
        return None, budget
    if p.returncode:
        raise ValueError("unsupported native support search")
    lines = p.stdout.splitlines()
    if not lines or not lines[-1].startswith("nodes "):
        raise ValueError("incomplete support enumeration")
    nodes = int(lines.pop().split()[1])
    patterns = set()
    for line in lines:
        ids = tuple(sorted(map(int, line.split())))
        if (
            not ids
            or len(ids) > width
            or len(set(ids)) != len(ids)
            or ids[0] < 0
            or ids[-1] >= len(cols)
        ):
            raise ValueError("invalid projected support")
        v = reduce(xor, (cols[i] for i in ids))
        if not v or v & ((1 << nd) - 1):
            raise ValueError("invalid projected logical")
        patterns.add(ids)
    return sorted(patterns), nodes


def reduce_vector(v, basis):
    while v and v.bit_length() in basis:
        v ^= basis[v.bit_length()]
    return v


def span(values):
    basis = {}
    for v in values:
        v = reduce_vector(v, basis)
        if v:
            basis[v.bit_length()] = v
    return basis


def xor_feasible(groups, target=0):
    """Can one value from each finite domain XOR to target? Exact branch search."""
    groups = list(groups)
    while groups:
        reps = [next(iter(g)) for g in groups]
        rhs = reduce(xor, reps, target)
        varies = [reduce(or_, (v ^ a for v in g), 0) for g, a in zip(groups, reps)]
        left = [0]
        for v in varies:
            left.append(left[-1] | v)
        if rhs & ~left[-1]:
            return False
        right = 0
        changed = False
        for j in reversed(range(len(groups))):
            fixed = ~(left[j] | right)
            right |= varies[j]
            if not varies[j] & fixed:
                continue
            keep = [v for v in groups[j] if not (rhs ^ reps[j] ^ v) & fixed]
            if not keep:
                return False
            if len(keep) != len(groups[j]):
                groups[j] = keep
                changed = True
        if not changed:
            break
    if not groups:
        return target == 0
    i = min(range(len(groups)), key=lambda j: len(groups[j]))
    group = groups.pop(i)
    return any(xor_feasible(groups, target ^ v) for v in group)


def lift_bound(exe, tags, nd, width, budget):
    """Exclude lifts up to width, including one possible zero-projection fault.

    A projected solution below width-1 could allow two or more extra faults.
    In that case only report the proven projected bound.
    """
    cols = sorted(set(tags) - {0})
    patterns, nodes = enumerate_supports(exe, cols, nd, width, budget)
    if patterns is None:
        return dict(outcome="budget", nodes=nodes, lower=1, spatial_lower=1)
    minimum = min(map(len, patterns), default=width + 1)
    stats = dict(
        nodes=nodes,
        spatial_lower=minimum,
        pattern_counts={str(k): v for k, v in Counter(map(len, patterns)).items()},
    )
    if minimum < width - 1 or minimum < 2:
        return dict(stats, outcome="spatial relaxation", lower=minimum)
    reps = {c: min(tags[c]) for c in cols}
    directions = {c: span(t ^ reps[c] for t in tags[c]).values() for c in cols}
    affine_passing = 0
    # Relax location exclusivity. A failed lift remains a valid lower bound.
    for ids in sorted(patterns, key=len):
        chosen = [cols[i] for i in ids]
        basis = span(chain.from_iterable(directions[c] for c in chosen))
        rhs = reduce(xor, (reps[c] for c in chosen))
        optional = tags[0] if len(ids) < width else [0]
        tails = [v for v in optional if not reduce_vector(v ^ rhs, basis)]
        if not tails:
            continue
        affine_passing += 1
        if xor_feasible([*(tags[c] for c in chosen), tails]):
            return dict(
                stats,
                outcome="possible lift",
                lower=len(chosen),
                affine_passing=affine_passing,
            )
    return dict(
        stats, outcome="excluded", lower=width + 1, affine_passing=affine_passing
    )
