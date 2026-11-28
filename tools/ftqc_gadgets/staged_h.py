"""Transversal H on the third Steane block, including idles on both controls."""

from ftqc_codes.codes import steane
from ftqc_codes.gf2 import Echelon
from ftqc_decoders.staged import Recovery
from .staged import endpoint, Counterexample


def check(spec, rule=None):
    if set(spec) != {"kind", "operation", "blocks"} or len(spec["blocks"]) != 3:
        raise ValueError("target H needs three ordered Steane blocks")
    code = steane()
    # H maps each positive CSS check to its complementary check and swaps the
    # positive logical X,Z. Equality on the Pauli generators fixes the channel.
    if any(h not in Echelon(code.hz) for h in code.hx) or any(
        h not in Echelon(code.hx) for h in code.hz
    ):
        raise ValueError("transversal H does not preserve the code")
    hx = [h << (7 * b) for b in range(3) for h in code.hx]
    hz = [h << (7 * b) for b in range(3) for h in code.hz]
    decoder = Recovery(hx, hz, rule)
    clean_failure = decoder.observe([], 0, 0, "clean")
    if clean_failure:
        raise Counterexample("target-H clean recovery", **clean_failure)
    for q in range(21):
        for x, z in ((1, 0), (1, 1), (0, 1)):
            for timing in ("incoming", "internal"):
                a, b = (z, x) if timing == "incoming" and q >= 14 else (x, z)
                collision = decoder.observe([], a << q, b << q, [timing, q, x, z])
                if collision:
                    raise Counterexample("target-H recovery collision", **collision)
    ports = endpoint("steane", spec["blocks"])
    return dict(
        status="proved",
        inputs=ports,
        outputs=ports,
        action="H-third",
        ideal=dict(
            target=spec["blocks"][2], logical_x_image="+Z", logical_z_image="+X"
        ),
        faults=dict(
            order=1,
            noisy_locations=21,
            alternatives={"incoming": 63, "H": 21, "idle": 42},
            recovery=decoder.rows(),
            ideal_resources=[],
            terminal="ideal full syndrome and recovery",
        ),
    )
