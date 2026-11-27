"""Source-bound finite surgery designs accepted by the Gadgets CLI."""

from .parsimonious import cellulate
from .surgery_geometry import adapter, distance, expansion
from .surgery_merge import merge
from .surgery_build import graph


def geometry(s):
    method = s["method"]
    if method == "linear":
        return s["ports"], s["aux"], s.get("faces"), {}
    sites = s["sites"]
    if len(set(sites)) != len(sites) or any(
        type(q) != int or not 0 <= q < s["n"] for q in sites
    ):
        raise ValueError("invalid or duplicate physical port")
    ports = [1 << q for q in sites]
    if method in ("enhanced", "bridge"):
        from .surgery_bridges import enhanced, bridge

        if method == "enhanced":
            g = enhanced(
                sites, s["hx"] if s["axis"] == "Z" else s["hz"], s.get("budget", 100000)
            )
        else:
            m = s["left_size"]
            g = bridge(m, len(sites) - m, s["left"], s["right"], s["links"])
        return ports, [(1 << a) | (1 << b) for a, b in g["edges"]], g.get("faces"), g
    if method == "adapter":
        if len(sites) % 2:
            raise ValueError("adapter needs equally many ports at each end")
        g = adapter(len(sites) // 2, s["left"], s["right"])
        return ports, [(1 << a) | (1 << b) for a, b in g["edges"]], g["faces"], g
    if method == "parsimonious":
        g = cellulate(len(sites), s["edges"])
        extended = [0] * g["vertices"]
        for q, v in zip(sites, g["ports"]):
            extended[v] = 1 << q
        return extended, [(1 << a) | (1 << b) for a, b in g["edges"]], g["faces"], g
    if method == "gauging":
        edges = s["edges"]
        if any(
            len(p) != 2
            or any(type(i) != int or not 0 <= i < len(sites) for i in p)
            or p[0] == p[1]
            for p in edges
        ):
            raise ValueError("invalid gauging edges")
        from ftqc_codes.gf2 import Echelon

        rows = [(1 << a) | (1 << b) for a, b in edges]
        if len(set(rows)) != len(rows) or Echelon(rows).rank != len(sites) - 1:
            raise ValueError("gauging requires a connected simple graph")
        return ports, rows, None, dict(edges=edges)
    raise ValueError("unknown finite surgery constructor")


def build(s):
    n = s["n"]
    if s["method"] == "branch-code":
        from .surgery_branch_code import cylinders

        return cylinders(n, s["hx"], s["hz"], s["targets"], s["rounds"])
    if s["method"] == "branch":
        from .surgery_parallel import branch

        return dict(
            spec=branch(n, s["checks"], s["targets"], s["rounds"]),
            geometry=dict(method="branch"),
        )
    if s["method"] == "homomorphic":
        g = graph(
            n, s["checks"], s["target"], s["edges"], "homomorphic", s.get("rounds", 1)
        )
        return dict(spec=g, geometry=dict(method="homomorphic"))
    p, a, faces, details = geometry(s)
    model = s.get("noise_locations", "physical")
    if model not in ("physical", "phenomenological"):
        raise ValueError("unknown surgery noise locations")
    cone, g = merge(
        n,
        s["hx"],
        s["hz"],
        p,
        a,
        s["rounds"],
        s["axis"],
        faces,
        model == "phenomenological",
    )
    if g["targets"] != s["targets"]:
        raise ValueError("induced subgroup basis differs from requested signed outputs")
    target = s.get("distance_target", 1)
    budget = s.get("budget", 100000)
    if type(target) != int or target < 1 or type(budget) != int or budget < 0:
        raise ValueError("invalid geometry budget")
    return dict(
        spec=g,
        geometry=dict(
            method=s["method"],
            construction=details,
            cone=cone,
            static=distance(cone["n"], cone["hx"], cone["hz"], target, budget),
            relative_expansion=expansion(p, a, target, budget),
        ),
    )
