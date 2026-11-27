"""Install the signed pivot protocol on BB144 and named external CSS source blocks."""

import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

from ftqc_codes.gf2 import supp
from ftqc_extractors.__main__ import digest
from ftqc_extractors.bicycle import layout, mapping, native, combine
from ftqc_extractors.injection import logicals
from ftqc_gadgets.bicycle_ports import generators, route, meter, verify_meter
from .insitu import factory
from .pivot import ideal_table
from .supports import fields

STAR = 319
FIRST_SOURCE = 320


def circuit(ticks, readouts=None, detectors=None):
    return dict(ticks=ticks, readouts=readouts or [], detectors=detectors or [])


@lru_cache(maxsize=512)
def footprint(blob):
    template = json.loads(blob)
    links = set()
    count = Counter()
    for tick in template["ticks"]:
        used = set()
        for name, *args in tick:
            sites = args[:1] if name == "Z_IF" else args
            if not sites or len(set(sites)) != len(sites) or used.intersection(sites):
                raise ValueError("physical tick collision")
            used.update(sites)
            if name == "CX":
                links.add(tuple(sorted(sites)))
            count[name] += 1
    return sorted(links), dict(count)


class Timeline:
    def __init__(self, latency, templates):
        self.latency = latency
        self.templates = {}
        self.inventory = {}
        self.expanded = templates
        self.steps = []
        self.now = 0
        self.ready = {}
        self.links = set()
        self.link_owner = {}

    def issue(self, label, options, *, controls=(), record=None, phase="factory"):
        # Conditional arms occupy the same reserved window. Conditions are
        # lookup tables indexed in the stated little-endian record order.
        if not set(controls) <= self.ready.keys() or len(options) != 1 << len(controls):
            raise ValueError("undefined or incomplete adaptive control")
        start = max([self.now] + [self.ready[k] for k in controls])
        ids = []
        for template in options:
            key = digest(template)
            if key not in self.inventory:
                blob = json.dumps(template, sort_keys=True, separators=(",", ":"))
                links, counts = footprint(blob)
                self.links.update(links)
                for pair in links:
                    self.link_owner.setdefault(pair, label)
                self.inventory[key] = dict(ticks=len(template["ticks"]), counts=counts)
                if self.expanded:
                    self.templates[key] = template
            if record and len(template["readouts"]) != 1:
                raise ValueError("missing public readout")
            ids.append(key)
        self.now = start + max(self.inventory[k]["ticks"] for k in ids)
        self.steps.append(
            dict(
                label=label,
                start=start,
                end=self.now,
                phase=phase,
                controls=list(controls),
                arms=ids,
                record=record,
            )
        )
        if record:
            if record in self.ready:
                raise ValueError("reused record")
            self.ready[record] = self.now + self.latency

    def check(self, requested):
        if requested is None:
            return
        if digest(requested) != digest(self.steps):
            raise ValueError("installed timing, placement or feed-forward mismatch")


@lru_cache(maxsize=256)
def port(x, z, sign=0):
    sites = supp(x | z)
    xx = sum(((x >> q) & 1) << j for j, q in enumerate(sites))
    zz = sum(((z >> q) & 1) << j for j, q in enumerate(sites))
    verify_meter(len(sites), (), xx, zz, sign)
    return meter(x, z, STAR, sign)


def pauli(x, z):
    # Products of X then Z agree up to a branch-global phase.
    return circuit(
        [[[axis, q] for q in supp(mask)] for axis, mask in (("X", x), ("Z", z)) if mask]
    )


def prepare_sources(t, source):
    # The native proof checks carrier Clifford, projection, fixed checks and the
    # frame. Serialize preparations and preserve each block until its last use.

    lanes = []
    acceptance = []
    base = FIRST_SOURCE
    for j, row in enumerate(source["inputs"]):
        model = row["spec"]["model"]
        exp = row["receipt"]["result"]
        spec = model["injection"]
        proof = exp["isometry"]
        n, hx, hz, pairs = logicals(spec)
        prefix = f"source.{j}"

        def add(label, g, **kw):
            t.issue(prefix + "." + label, [g], phase="preparation", **kw)

        init = (
            [["T", base + q] for q in proof["carriers"]]
            + [["RX", base + q] for q in proof["plus"]]
            + [["R", base + q] for q in proof["zero"]]
        )
        for a, b in proof["bells"]:
            init.extend([["RX", base + a], ["R", base + b]])
        ticks = [init]
        if proof["bells"]:
            ticks.append([["CX", base + a, base + b] for a, b in proof["bells"]])
        ticks.extend(
            [[op, *[base + q for q in qs]]] for op, *qs in proof["carrier_gates"]
        )
        add("prepare", circuit(ticks))
        for axis, checks in (("X", hx), ("Z", hz)):
            syndromes = {}
            for i, h in enumerate(checks):
                rec = f"{prefix}.{axis}{i}"
                syndromes[i] = rec
                add(
                    axis + str(i),
                    port(
                        (h << base) if axis == "X" else 0,
                        (h << base) if axis == "Z" else 0,
                    ),
                    record=rec,
                )
            for i, mask in proof["syndrome_frame"]["Z" if axis == "X" else "X"]:
                corr = pauli(0, mask << base) if axis == "X" else pauli(mask << base, 0)
                t.issue(
                    prefix + ".syndrome-frame",
                    [circuit([]), corr],
                    controls=[syndromes[i]],
                    phase="preparation",
                )
        acceptance.append(
            dict(
                attempt=row["spec"]["attempt"],
                block=row["spec"]["block"],
                require_zero=[
                    f"{prefix}.{axis}{i}"
                    for axis in ("X", "Z")
                    for i in proof["postselect"][axis]
                ],
                reject_if_event=[e["id"] for e in model["events"] if e["reject"]],
            )
        )
        lanes.append((base, n, pairs))
        base += n
    return lanes, base, acceptance


def install(s, expanded=False):
    fields(
        s, "kind", "source", "latency", *(k for k in ("schedule", "couplers") if k in s)
    )
    if (
        s["source"]["kind"] != "insitu-rm15"
        or type(s["latency"]) is not int
        or s["latency"] < 0
    ):
        raise ValueError("expected RM15 inputs and a nonnegative latency")
    src = factory(s["source"])["factory"]
    if len({b["block"] for b in src["inputs"]}) != len(src["inputs"]):
        raise ValueError("input blocks must be disjoint from each other")
    if any(b["block"] == "BB144/factory" for b in src["inputs"]):
        raise ValueError("source preparation aliases live factory memory")
    m = mapping()
    truth = ideal_table(m)
    _, xs, zs = layout()
    t = Timeline(s["latency"], expanded)
    lanes, width, source_conditions = prepare_sources(t, s["source"])
    # Resolve the source verification and heralding before use of any member.
    pre_end = t.now + s["latency"]
    t.now = pre_end
    inputs = {}
    for j, entry in enumerate(src["inputs"]):
        for member, col in entry["mapping"]:
            base, _, pairs = lanes[j]
            x, z, sx, sz = pairs[member]
            inputs[col] = (x << base, z << base, sx, sz, j)
    empty = circuit([])
    for slot in range(7):
        if slot == 5:
            continue
        xinit = slot < 5
        tag = f"init.{slot}"
        t.issue(
            tag, [port(xs[slot] if xinit else 0, 0 if xinit else zs[slot])], record=tag
        )
        t.issue(
            tag + ".correct",
            [empty, pauli(0, zs[slot]) if xinit else pauli(xs[slot], 0)],
            controls=[tag],
        )
    base_native = native()
    merged = circuit(
        base_native["ticks"], base_native["readouts"], base_native["detectors"]
    )
    pivot_x, pivot_z = xs[5], zs[5]
    location = []
    for row in m["rotations"]:
        i = row["column"]
        prefix = f"r{i}."
        a, b, c, d = [prefix + q for q in "abcd"]
        rx, rz, sx, sz, lane = inputs[i]
        t.issue(prefix + "pivot-plus", [port(pivot_x, 0)], record=prefix + "init")
        t.issue(
            prefix + "pivot-correct",
            [empty, pauli(0, pivot_z)],
            controls=[prefix + "init"],
        )
        t.issue(prefix + "source-pivot-ZZ", [port(0, rz | pivot_z, sz)], record=a)
        dx, dy = row["shift"]
        for k in route(-dx, -dy):
            t.issue(prefix + "shift", [circuit(generators()[k])])
        t.issue(prefix + "Gross-Z", [merged], record=c)
        for k in route(dx, dy):
            t.issue(prefix + "return", [circuit(generators()[k])])
        raw = port(rx, 0, sx)
        t.issue(prefix + "source-X", [raw], record=b)
        index = next(j for j, tick in enumerate(raw["ticks"]) if ["M", STAR] in tick)
        location.append(
            dict(
                column=i,
                step=prefix + "source-X",
                tick=index,
                wire=STAR,
                before="M",
                error="X",
                record=b,
            )
        )
        cross = (pivot_x & pivot_z).bit_count()
        y = port(pivot_x, pivot_z, ((1 - cross) // 2) % 2)
        x = port(pivot_x, 0)
        t.issue(prefix + "pivot-read", [x, y, y, x], controls=[a, c], record=d)
        correction = pauli(0, combine(zs, row["parity"]))
        arms = [None] * 16
        for branch in truth["branches"]:
            key = sum(v << j for j, v in enumerate(branch["outcomes"]))
            arms[key] = correction if branch["frame"] else empty
        t.issue(prefix + "frame", arms, controls=[a, b, c, d])
    # S on the output, using coherent parity accumulation and uncomputation.
    parity = [["CX", q, STAR] for q in supp(zs[0])]
    t.issue(
        "output-S",
        [
            circuit(
                [
                    [["R", STAR]],
                    *[[g] for g in parity],
                    [["S", STAR]],
                    *[[g] for g in reversed(parity)],
                    [["M", STAR]],
                ]
            )
        ],
    )
    for j in range(1, 5):
        t.issue(f"check.{j}", [port(xs[j], 0)], record=f"check.{j}")
    t.now = max(t.ready.values())
    t.check(s.get("schedule"))
    links = sorted(t.links)
    if "couplers" in s:
        offered = {tuple(sorted(pair)) for pair in s["couplers"]}
        if set(links) - offered:
            pair = min(set(links) - offered)
            raise ValueError(
                f"uninstalled physical coupler {pair} at {t.link_owner[pair]}"
            )
    result = dict(
        status="proved",
        map_sha256=digest(m),
        factory=src,
        ideal=truth,
        source_blocks=[b["block"] for b in src["inputs"]],
        width=width,
        factory_width=FIRST_SOURCE,
        preparation_end=pre_end,
        source_acceptance=source_conditions,
        elapsed=t.now,
        installed_couplers=links,
        steps=t.steps,
        inventory=t.inventory,
        output=dict(
            block="BB144/factory", slot=0, x=xs[0], z=zs[0], frame="canonical-positive"
        ),
        postselect=[f"check.{j}" for j in range(1, 5)],
        fault_locations=location,
        scope="clean encoded memory; source T preparations and named carrier events; full signed pivot circuit, Tanner shifts, Gross Z and additional parity-meter couplers; padded conditional windows",
    )
    if expanded:
        result["templates"] = t.templates
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args()
    args.o.write_text(
        json.dumps(
            install(json.loads(args.study.read_text()), True),
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    )
