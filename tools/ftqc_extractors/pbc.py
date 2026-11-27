"""Compile Clifford+T to signed measurements on a fixed extractor block map."""

from itertools import product
from ftqc_faults.zx.frame_oracle import _apply_1q_to_stab_list, _apply_2q_to_stab_list
from .algebra import Limit
from .block_map import Hardware, physical


def backprop(axis, slot, history, k):
    vec = 1 << (2 * slot + (axis == "X"))
    state = [(vec, 0)]
    for gate, *qs in reversed(history):
        if gate == "CX":
            state = _apply_2q_to_stab_list(state, *qs, "CX")
        else:
            state = _apply_1q_to_stab_list(state, qs[0], "SDG" if gate == "S" else gate)
    v, sign = state[0]
    return (
        sum(((v >> (2 * i + 1)) & 1) << i for i in range(k)),
        sum(((v >> (2 * i)) & 1) << i for i in range(k)),
        sign,
    )


def analyze(s, budget):
    hw = Hardware(s, budget)
    blocks = hw.blocks
    supplies = s["magic_supplies"]
    if len(supplies) != len(blocks) or any(
        type(i) is not int or not 0 <= i < len(blocks) for i in supplies
    ):
        raise ValueError("one magic-cache block index per block")
    k = [len(b["basis"]) for b in blocks]
    workspace = [(i, j) for i, n in enumerate(k) for j in range(n - 1)]
    settings = {}
    gates = s["gates"]
    rotations = 0
    for gate, *qs in gates:
        if gate not in ("H", "S", "X", "Y", "Z", "CX", "T") or len(qs) != (
            2 if gate == "CX" else 1
        ):
            raise ValueError("invalid Clifford+T gate")
        if any(tuple(q) not in workspace for q in qs) or (
            len(qs) == 2 and qs[0] == qs[1]
        ):
            raise ValueError("gate does not act on distinct workspace qubits")
        rotations += int(gate == "T" or (gate == "CX" and qs[0][0] != qs[1][0]))
    requested = s.get("outcomes")
    if requested is None:
        if 2 * rotations >= budget.bit_length() or 4**rotations > budget:
            raise Limit(
                "adaptive branch enumeration; supply an explicit outcome transcript"
            )
        transcripts = product(range(4), repeat=rotations)
    else:
        if len(requested) != rotations or any(
            type(v) is not int or not 0 <= v < 4 for v in requested
        ):
            raise ValueError("one two-bit outcome per rotation")
        transcripts = [requested]
    branches = []
    for transcript in transcripts:
        history = [[] for _ in blocks]
        time = [0] * len(blocks)
        program = []

        def issue(terms, label):
            ids = tuple(sorted(terms))
            reps = [physical(blocks[i], *terms[i]) for i in ids]
            logical = ("-" if sum(w[0] == "-" for w in reps) % 2 else "+") + "".join(
                w[1:] for w in reps
            )
            cert = hw.select(ids, logical)
            name = ",".join(map(str, ids)) + "/" + logical
            if name not in settings:
                settings[name] = dict(
                    blocks=list(ids),
                    logical=logical,
                    checks=cert["merged_checks"],
                    readout_sign=cert["readout_sign"],
                    split_sites=cert["split_sites"],
                    split_X_paths=cert["split_X_paths"],
                )
            start = max(time[i] for i in ids)
            for i in ids:
                time[i] = start + 1
            program.append(dict(cycle=start, measurement=name, role=label))

        def axis(i, j, p):
            return backprop(p, j, history[i], k[i])

        index = 0
        for gate, *qs in gates:
            i, j = qs[0]
            if gate != "T" and (gate != "CX" or qs[1][0] == i):
                history[i].append([gate, *[q[1] for q in qs]])
                continue
            raw = transcript[index]
            b, m = raw & 1, raw >> 1
            index += 1
            terms = {i: axis(i, j, "Z")}
            cache = supplies[i] if gate == "T" else i
            scratch = k[cache] - 1
            if gate == "CX":
                u, v = qs[1]
                terms[u] = axis(u, v, "X")
                issue(
                    {cache: (1 << scratch, 1 << scratch, 0)},
                    "prepare -Y; flip with Z on a +Y result",
                )
            else:
                program.append(
                    dict(
                        cycle=time[cache],
                        resource="fresh independent encoded |T>",
                        block=cache,
                        slot=scratch,
                    )
                )
            x, z, sgn = terms.get(cache, (0, 0, 0))
            terms[cache] = (x, z ^ (1 << scratch), sgn)
            issue(terms, "rotation parity b")
            issue({cache: (1 << scratch, 0, 0)}, "resource X m")
            # The local frames depend on both bits, including the cache's X
            # readout. A disjoint rotation may still start in the same cycle.
            for block in terms:
                time[block] = max(time[block], time[cache])
            if gate == "T":
                if b:
                    history[i].append(["S", j])
                if m:
                    history[i].append(["Z", j])
            else:
                if b ^ m:
                    history[i].append(["Z", j])
                    history[u].append(["X", v])
                history[i].append(["S", j])
                history[u].extend([["H", v], ["S", v], ["H", v]])
        for i, j in workspace:
            issue({i: axis(i, j, "Z")}, f"output Z({i},{j})")
        branches.append(
            dict(outcomes=list(transcript), program=program, cycles=max(time))
        )
    inventory = []
    for ids, (data, checks) in hw.labels.items():
        inventory.append(
            dict(
                blocks=list(ids),
                data=data,
                checks=checks,
                couplers=next(
                    v["physical_check_couplers"]
                    for (key, _), v in hw.cache.items()
                    if key == ids
                ),
            )
        )
    result=dict(
        status="proved",
        model="conditional Pauli-based compilation on a fixed block map",
        branch_scope="all" if requested is None else "specified transcript",
        rotations=rotations,
        settings=settings,
        hardware=inventory,
        branches=branches,
    )

    if "surgery_installation" in s:
        from .pbc_surgery import bind
        result["surgery_installation"]=bind(result,s,budget)
    return result
