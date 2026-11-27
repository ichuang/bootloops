"""Serialize source attempts, a live-block SWAP and native PBC operations."""

from collections import Counter, defaultdict
from ftqc_gadgets.occupied_transfer import occupied_swap
from ftqc_magic.delivery import read_delivery


def tick(i):
    if type(i) is not int or i < 0:
        raise ValueError("expected a nonnegative integer tick or wire")
    return i


def bind(compiled, source):
    install = source["magic_delivery"]
    if set(install) != {
        "factory",
        "receipt",
        "attempts",
        "requests",
        "couplers",
        "latency",
    }:
        raise ValueError("invalid delivery inventory")
    evidence = read_delivery(install["factory"], install["receipt"])
    if (
        evidence["groups"] != [[0]]
        or evidence["state"] != "T"
        or not evidence["acceptance_exact"]
    ):
        raise ValueError("PBC needs one T with a replayed exact acceptance law")
    latency = tick(install["latency"])
    blocks = source["blocks"]
    offsets = [sum(b["n"] for b in blocks[:i]) for i in range(len(blocks))]
    memory = sum(b["n"] for b in blocks)
    menu = compiled.get("surgery_installation", {})
    installed = menu.get("selected", {})
    links = set()
    for a, b in install["couplers"]:
        if tick(a) == tick(b):
            raise ValueError("invalid transfer coupler")
        links.add(tuple(sorted((a, b))))
    if (
        max(
            Counter(
                q
                for link in links | {tuple(p) for p in menu.get("couplers", [])}
                for q in link
            ).values(),
            default=0,
        )
        > source["max_degree"]
    ):
        raise ValueError("combined installed coupler degree too high")
    attempts = {}
    base = []
    for a in install["attempts"]:
        if set(a) != {"id", "start", "end", "wires", "raw", "accepted"}:
            raise ValueError("invalid factory attempt")
        key = a["id"]
        wires = [tick(q) for q in a["wires"]]
        if (
            not isinstance(key, str)
            or not key
            or key in attempts
            or len(set(wires)) != len(wires)
            or len(wires) != evidence["active_qubits"]
            or min(wires) < memory
            or a["raw"] not in wires
            or type(a["accepted"]) is not bool
            or tick(a["start"]) >= tick(a["end"])
        ):
            raise ValueError("invalid factory identity, wires or interval")
        if not a["accepted"] and evidence["acceptance_lower"] == "1":
            raise ValueError("impossible rejected attempt")
        attempts[key] = a
        base.append((a["start"], a["end"], set(wires), "factory"))
    prepared = {}
    ledgers = []
    for branch in compiled["branches"]:
        spent = set()
        reserved = list(base)
        steps = []
        now = used = buffer = active_memory = 0
        physical_idle = 0
        for op in branch["program"]:
            if "resource" in op:
                if used >= len(install["requests"]):
                    raise ValueError("unbound PBC request")
                req = install["requests"][used]
                used += 1
                if set(req) != {"attempt", "block", "slot", "start", "member"} or (
                    req["block"],
                    req["slot"],
                    req["member"],
                ) != (op["block"], op["slot"], 0):
                    raise ValueError("wrong destination slot or batch member")
                key = req["attempt"]
                if key not in attempts or key in spent or not attempts[key]["accepted"]:
                    raise ValueError("unavailable or already consumed attempt")
                spent.add(key)
                a = attempts[key]
                start = tick(req["start"])
                if start < a["end"] + latency or start < now:
                    raise ValueError(
                        "early delivery or collision with prior measurement"
                    )
                block, slot = op["block"], op["slot"]
                if (block, slot) not in prepared:
                    prepared[block, slot] = occupied_swap(blocks[block], slot)
                transfer = prepared[block, slot]
                raw = a["raw"]
                reserved.append((a["end"], start, {raw}, "buffer"))
                buffer += start - a["end"]
                live_ops = []
                for _, c, t in transfer["operations"]:
                    c, t = (
                        raw if q == transfer["raw"] else offsets[block] + q
                        for q in (c, t)
                    )
                    if tuple(sorted((c, t))) not in links:
                        raise ValueError("missing installed transfer coupler")
                    at = start + len(live_ops)
                    reserved.append((at, at + 1, {c, t}, "swap"))
                    live_ops.append([at, "CX", c, t])
                    active_memory += 1
                now = start + len(live_ops)
                steps.append(
                    dict(
                        attempt=key,
                        batch_member=0,
                        raw=raw,
                        block=block,
                        slot=slot,
                        ready_tick=a["end"] + latency,
                        consume_tick=now,
                        operations=live_ops,
                    )
                )
            else:
                setting = compiled["settings"][op["measurement"]]
                duration = installed[op["measurement"]]["ticks"] if installed else 1
                sites = {
                    q
                    for i in setting["blocks"]
                    for q in range(offsets[i], offsets[i] + blocks[i]["n"])
                }
                # Installed circuits include the full memory in the native schedule.
                if installed:
                    part = source["surgery_installation"]["settings"][
                        installed[op["measurement"]]["setting"]
                    ]
                    pos = part["placement"]
                    sites = set(range(memory)) | set(pos)
                    native = part["gadget"]["spec"]["ticks"]
                    if len(native) != duration:
                        raise ValueError("native surgery duration mismatch")
                    for row in native:
                        busy = {
                            pos[q]
                            for name, *args in row
                            if name.lstrip("!") != "I"
                            for q in (args[:1] if name.lstrip("!") == "Z_IF" else args)
                        }
                        physical_idle += memory - sum(q < memory for q in busy)
                reserved.append((now, now + duration, sites, "measurement"))
                active_memory += duration * sum(q < memory for q in sites)
                steps.append(
                    dict(
                        measurement=op["measurement"],
                        role=op["role"],
                        start=now,
                        ticks=duration,
                    )
                )
                if installed:
                    op["surgery_start_tick"] = now
                now += duration + latency
        if used != len(install["requests"]) or spent != {
            k for k, a in attempts.items() if a["accepted"]
        }:
            raise ValueError("every accepted batch member needs exactly one request")
        events = defaultdict(list)
        for start, end, sites, role in reserved:
            if start < end:
                events[start].append((1, sites, role))
                events[end].append((-1, sites, role))
        occupied = set()
        peak = memory
        for at, changes in sorted(events.items()):
            for sign, sites, role in sorted(changes, key=lambda item: item[0]):
                if sign < 0:
                    occupied.difference_update(sites)
                elif sites & occupied:
                    raise ValueError(
                        f"wire or installed surgery collision at tick {at} ({role})"
                    )
                else:
                    occupied.update(sites)
            peak = max(peak, memory + sum(q >= memory for q in occupied))
        elapsed = max([now] + [a["end"] for a in attempts.values()])
        if installed:
            branch["surgery_ticks"] = now
        ledgers.append(
            dict(
                outcomes=branch["outcomes"],
                steps=steps,
                elapsed_ticks=elapsed,
                peak_reserved_qubits=peak,
                buffer_wire_ticks=buffer,
                unreserved_memory_wire_ticks=elapsed * memory - active_memory,
                memory_idle_wire_ticks=(
                    elapsed * memory - active_memory + physical_idle
                    if installed
                    else None
                ),
            )
        )
    return dict(
        schema="ftqc-extractors-magic-installation-1",
        source=evidence,
        attempts=install["attempts"],
        transfers=[dict(block=b, proof=v) for (b, _), v in prepared.items()],
        branches=ledgers,
        timing="serialized PBC; declared factory intervals; ideal settings take one model tick unless physically installed",
        factory_cost="producer expected wire-ticks/output includes independent retries; transcript occupancy is realized",
        fault_scope="ideal free destination slot; arbitrary spectators; whole-program noisy schedule open",
    )
