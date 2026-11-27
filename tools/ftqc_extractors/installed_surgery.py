"""Install checked measurement instruments on one explicit coupler inventory."""

from ftqc_gadgets.surgery import receipt, verify, pauli


def compile_schedule(s):
    n = s["n"]
    if type(n) != int or n < 1 or type(s["idles"]) != bool:
        raise ValueError("invalid hardware memory")
    links = set()
    for pair in s["couplers"]:
        if (
            len(pair) != 2
            or any(type(q) != int or q < 0 for q in pair)
            or pair[0] == pair[1]
        ):
            raise ValueError("invalid hardware coupler")
        links.add(tuple(sorted(pair)))
    ticks, detectors, readouts, targets = [], [], [], []
    nr = 0
    peak = n
    for batch in s["batches"]:
        if not batch:
            raise ValueError("empty parallel batch")
        used = set()
        parts = []
        for setting in batch:
            if set(setting) - {"gadget", "placement", "control"}:
                raise ValueError("unknown installation field")
            proof = setting["gadget"]
            verify(proof)
            g = proof["spec"]
            pos = setting["placement"]
            width = 1 + max(
                [g["n"] - 1]
                + [
                    q
                    for row in g["ticks"]
                    for op in row
                    for q in (op[1:2] if op[0].lstrip("!") == "Z_IF" else op[1:])
                ]
            )
            if (
                len(pos) != width
                or any(type(q) != int or q < 0 for q in pos)
                or len(set(pos)) != len(pos)
            ):
                raise ValueError("incomplete or aliased placement")
            if any((q < n) != (j < g["n"]) for j, q in enumerate(pos)):
                raise ValueError("placement changes the data/ancilla partition")
            touched = {
                pos[q]
                for row in g["ticks"]
                for op in row
                for q in (op[1:2] if op[0].lstrip("!") == "Z_IF" else op[1:])
            }

            def embed(word):
                out = ["I"] * n
                for j, c in enumerate(word[1:]):
                    out[pos[j]] = c
                return word[0] + "".join(out)

            words = [embed(word) for word in g["targets"]]
            control = setting.get("control")
            guard = []
            if control is not None:
                index, sites = control
                sites = [sites] if type(sites) == int else sites
                if (
                    type(index) != int
                    or not 0 <= index < len(readouts)
                    or not isinstance(sites, list)
                    or not sites
                    or len(set(sites)) != len(sites)
                    or any(type(q) != int or q not in pos[: g["n"]] for q in sites)
                ):
                    raise ValueError("invalid earlier bit or conditional Z site")
                ctl = readouts[index]
                guard = [[["Z_IF", q, r] for q in sites] for r in ctl["records"]]
                if ctl["sign"]:
                    guard += [[["S", q] for q in sites] for _ in range(2)]
                touched.update(sites)
            if used & touched:
                raise ValueError(
                    "parallel settings need disjoint physical gate supports"
                )
            used.update(touched)
            parts.append(
                dict(
                    g=g, pos=pos, words=words, guard=guard, control=control, records=[]
                )
            )
        live = set(range(n))
        lengths = [2 * len(p["guard"]) + len(p["g"]["ticks"]) for p in parts]
        for t in range(max(lengths)):
            row = []
            for p, length in zip(parts, lengths):
                if t >= length:
                    continue
                guard, pos = p["guard"], p["pos"]
                j = t - len(guard)
                local = p["g"]["ticks"]
                if j < 0 or j >= len(local):
                    row += guard[t if j < 0 else j - len(local)]
                    continue
                for name, *args in local[j]:
                    gate = name.lstrip("!")
                    sites = [pos[i] for i in (args[:1] if gate == "Z_IF" else args)]
                    if gate == "CX" and tuple(sorted(sites)) not in links:
                        raise ValueError("surgery requires an uninstalled coupler")
                    if gate == "Z_IF":
                        sites.append(p["records"][args[1]])
                    row.append([name, *sites])
                    if gate in ("M", "MZ", "MX"):
                        p["records"].append(nr)
                        nr += 1
            busy = {
                q
                for name, *args in row
                for q in (args[:1] if name.lstrip("!") == "Z_IF" else args)
            }
            live.update(op[1] for op in row if op[0].lstrip("!") in ("R", "RX"))
            peak = max(peak, len(live))
            if s["idles"]:
                row += [["I", q] for q in sorted(live - busy)]
            live.difference_update(
                op[1] for op in row if op[0].lstrip("!") in ("M", "MZ", "MX")
            )
            ticks.append(row)
        for p in parts:
            g, ids = p["g"], p["records"]
            detectors.extend([[ids[r] for r in d] for d in g["detectors"]])
            for word, read in zip(p["words"], g["readouts"]):
                bits = {ids[r] for r in read["records"]}
                sign = read["sign"]
                if p["control"]:
                    index, sites = p["control"]
                    sites = [sites] if type(sites) == int else sites
                    if (pauli(word, n).x & sum(1 << q for q in sites)).bit_count() % 2:
                        bits.symmetric_difference_update(readouts[index]["records"])
                        sign ^= readouts[index]["sign"]
                targets.append(word)
                readouts.append(dict(records=sorted(bits), sign=sign))
    if targets != s["targets"]:
        raise ValueError(
            "installed target order, placement or sign differs from request"
        )
    return dict(
        n=n,
        checks=s["checks"],
        targets=targets,
        ticks=ticks,
        detectors=detectors,
        readouts=readouts,
        terminal_syndrome=s["terminal_syndrome"],
    ), peak


def analyze(s, budget):
    spec, peak = compile_schedule(s)
    proof = receipt(spec, s["target"] - 1, max(1, budget))
    timing = {}
    if "transport" in s:
        from .surgery_layout import transport

        layout = s["transport"]
        timing["transport"] = transport(
            layout["initial"], layout["layers"], spec["ticks"], layout["gate_range"]
        )
    return dict(
        **timing,
        status=proof["faults"]["status"],
        gadget=proof,
        peak_live=peak,
        installed_couplers=len({tuple(sorted(p)) for p in s["couplers"]}),
        ticks=len(spec["ticks"]),
    )
