"""Choose X/Z bases off the logical supports using the Appendix G integer objectives."""

import argparse
from functools import reduce
import json
from operator import or_
from pathlib import Path
from .algebra import bits
from .injection import analyze, logicals, mask, subsystem_error


def optimize(s, objective):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix

    if (
        objective not in ("fixed", "protected")
        or s["method"] != "selected"
        or analyze(s)["status"] != "proved"
    ):
        raise ValueError(
            "a verified selected source and a named objective are required"
        )
    n, hx, hz, pairs = logicals(s)
    p = s["preparation"]
    constrained = reduce(or_, (x | z for x, z, _, _ in pairs), 0)
    excluded = mask(p["carriers"] + sum(p["bells"], []))
    stable = constrained | excluded
    free = list(bits(((1 << n) - 1) & ~stable))
    free_index = {q: i for i, q in enumerate(free)}
    forced_x = [q for q in p["plus"] if stable >> q & 1]
    forced_z = [q for q in p["zero"] if stable >> q & 1]
    cols = len(free)
    constraints = []
    flags = []
    cover = {}
    for axis, rows, wrong, targets in (
        ("X", hx, forced_z, forced_x + [a for a, b in p["bells"]]),
        ("Z", hz, forced_x, forced_z + [b for a, b in p["bells"]]),
    ):
        for row in rows:
            if row & mask(p["carriers"]) or subsystem_error(row, wrong, p["bells"]):
                continue
            f = cols
            cols += 1
            flags.append(f)
            for q in bits(row):
                if q in free_index:
                    coef = -1 if axis == "X" else 1
                    constraints.append(({f: 1, free_index[q]: coef}, int(axis == "Z")))
                if q in targets and constrained >> q & 1:
                    cover.setdefault(q, []).append(f)
    protected = []
    for checks in cover.values():
        v = cols
        protected.append(v)
        cols += 1
        constraints.append(({v: 1, **{f: -1 for f in checks}}, 0))
    target = flags if objective == "fixed" else protected
    if not cols:
        return dict(s, postselect=analyze(s)["fixed_checks"])
    loss = np.zeros(cols)
    loss[target] = -1
    ii, jj, vv = [], [], []
    for i, (coeff, _) in enumerate(constraints):
        for j, v in coeff.items():
            ii.append(i)
            jj.append(j)
            vv.append(v)
    mat = coo_matrix((vv, (ii, jj)), shape=(len(constraints), cols)).tocsc()
    res = milp(
        loss,
        integrality=np.ones(cols),
        bounds=Bounds(0, 1),
        constraints=LinearConstraint(mat, -np.inf, [rhs for _, rhs in constraints]),
        options={"mip_rel_gap": 0},
    )
    if not res.success:
        raise ValueError(f"integer optimizer: {res.message}")
    plus = [q for j, q in enumerate(free) if res.x[j] > 0.5]
    prep = dict(
        p,
        plus=sorted(forced_x + plus),
        zero=sorted(forced_z + [q for q in free if q not in plus]),
    )
    candidate = {**s, "preparation": prep}
    candidate.pop("postselect", None)
    evidence = analyze(candidate)
    if evidence["status"] != "proved":
        raise ValueError("optimizer proposal failed source verification")
    candidate["postselect"] = evidence["fixed_checks"]
    return candidate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--objective", choices=["protected", "fixed"], required=True)
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args()
    proposed = optimize(json.loads(args.source.read_text()), args.objective)
    args.o.write_text(json.dumps(proposed, indent=2) + "\n")


if __name__ == "__main__":
    main()
