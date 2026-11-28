"""Verify a code action on a hardware graph with data sites first and optional spectators."""
from ftqc_codes import action_domain, automorphisms
from ftqc_faults.action import check_data_action, single_fault_spread


def analyze(spec):
    if (set(spec) != {"kind", "code", "action", "width", "ticks", "couplers"}
            or spec["kind"] != "installed-action"):
        raise ValueError("expected an installed-action specification")
    code = action_domain.load(spec["code"])
    action = automorphisms.verify(code, spec["action"])
    width, ticks, edges = spec["width"], spec["ticks"], spec["couplers"]
    if type(width) is not int or width < code.n or not isinstance(ticks, list):
        raise ValueError("invalid installed width or ticks")
    if not isinstance(edges, list) or any(
        not isinstance(edge, list) or len(edge) != 2 or edge[0] == edge[1]
        or any(type(q) is not int or not 0 <= q < width for q in edge) for edge in edges
    ):
        raise ValueError("expected physical coupler pairs")
    available = {tuple(sorted(edge)) for edge in edges}
    check_data_action(width, code.n, ticks, action["steps"])
    used = {tuple(sorted(op[1:])) for tick in ticks for op in tick if len(op) == 3}
    if used - available:
        raise ValueError(f"uninstalled coupler {min(used - available)}")
    counts = {}
    idle = width * len(ticks)
    for tick in ticks:
        for gate, *sites in tick:
            counts[gate] = counts.get(gate, 0) + 1
            idle -= len(sites)
    return dict(status="proved", scope="signed unitary data action", code_sha256=code.sha256(),
                basis_sha256=action["basis_sha256"], ticks=len(ticks), width=width,
                n_data=code.n, installed_couplers=len(available), used_couplers=len(used),
                gate_counts=counts, idle_qubit_ticks=idle, faults=single_fault_spread(width, code.n, ticks))

