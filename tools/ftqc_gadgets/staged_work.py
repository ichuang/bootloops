"""Resource counts for exactly the native locations and ideal boundaries tested."""

from ftqc_codes.codes import steane


def workload(spec, result):
    op = spec["operation"]
    nloc = result["faults"]["noisy_locations"]
    if op == "pieceable-ccz":
        widths = [h.bit_count() for h in steane().hz] * 3
        # Per check: one data interaction or idle per site, and a cat-output
        # noise site and one X readout on each of its legs.
        extra = 6 * sum(21 + 2 * w for w in widths)
        ngates = sum(map(len, spec["pieces"]))
        data, anc = 21, max(widths)
        clean, maximum = ngates + 6 * len(widths), ngates + 12 * len(widths)
        clean_reactions, max_reactions = 6, 12
    elif op == "rm15-phase":
        data, anc, extra = 15, 0, 0
        clean = maximum = 1
        clean_reactions = max_reactions = 0
    else:
        data, anc, extra = 15, max(v.bit_count() for v in spec["checks"]), 0
        clean = maximum = 9
        clean_reactions = max_reactions = 1
    return dict(
        data=data,
        peak_cat_legs=anc,
        clean_ticks=clean,
        max_ticks=maximum,
        clean_reactions=clean_reactions,
        max_reactions=max_reactions,
        clean_locations=nloc,
        max_locations=nloc + extra,
        tick="native layer or cat check with supplied cat; terminal EC priced separately",
    )
