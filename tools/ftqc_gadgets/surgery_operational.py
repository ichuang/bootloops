"""Operational evidence for a checked signed instrument; all laws name their boundary."""

from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from .surgery import Instrument
from .surgery_faults import columns, residual_profile


def conditional_law(groups, probability, cap, initial=(0, 0, 0)):
    """Convolve independent physical sites, one categorical alternative per site."""
    p = Fraction(probability)
    if not 0 <= p <= 1 or type(cap) != int or cap < 1:
        raise ValueError("invalid channel or state cap")
    zero = (0, 0, 0)
    law = {initial: Fraction(1)}
    for i, group in enumerate(groups):
        local = defaultdict(Fraction)
        local[zero] = 1 - p
        for _, effect in group:
            local[effect] += p / len(group)
        updated = defaultdict(Fraction)
        for a, pa in law.items():
            for b, pb in local.items():
                if pb:
                    updated[tuple(x ^ y for x, y in zip(a, b))] += pa * pb
                    if len(updated) > cap:
                        return dict(
                            status="inconclusive", completed_locations=i, cap=cap
                        )
        law = updated
    keep = sum(v for (det, _, _), v in law.items() if not det)
    return dict(
        status="exact",
        locations=len(groups),
        acceptance=str(keep),
        conditional=[
            dict(readout=b, residual=r, p=str(law[(d, b, r)] / keep))
            for d, b, r in sorted(law)
            if not d and keep
        ],
    )


def sampled(spec, probability, shots, seed, decoder="signature", prefix=None):
    import stim
    import numpy as np
    from .surgery_stim import compile_probe
    from ftqc_decoders.dem import DecodingProblem
    from ftqc_decoders.effective_distance import run_decoder

    inst = Instrument(spec)
    if shots < 1 or not 0 < Fraction(probability) < 1:
        raise ValueError("expected positive shots and nonzero physical rate")
    probe = compile_probe(inst.c, spec, inst.initial, float(Fraction(probability)))
    dem = probe.detector_error_model().flattened()
    if prefix:
        Path(str(prefix) + ".stim").write_text(str(probe) + "\n")
        Path(str(prefix) + ".dem").write_text(str(dem) + "\n")
    model = DecodingProblem.from_stim_dem_text(str(dem), merge_parallel=False)
    ds, obs = probe.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True, bit_packed=True
    )
    guesses = {}
    public_guesses = {}
    joint = defaultdict(Fraction)
    marginal = defaultdict(Fraction)
    # Last two Choi boundary parities per target are post-split data/readout
    # and untouched reference/readout. The second gives the public-bit error.
    offset = probe.num_observables - 2 * len(inst.targets)

    def public(mask):
        return sum(
            ((mask >> (offset + 2 * j + 1)) & 1) << j for j in range(len(inst.targets))
        )

    for d, o, p in zip(model.mech_detectors, model.mech_observables, model.priors):
        if d:
            joint[(d, o)] += p
            marginal[(d, public(o))] += p
    for d, o in joint:
        if d not in guesses or joint[d, o] > joint[d, guesses[d]]:
            guesses[d] = o
    for d, o in marginal:
        if d not in public_guesses or marginal[d, o] > marginal[d, public_guesses[d]]:
            public_guesses[d] = o
    cache = {}
    wrong = wrong_read = modular = accepted = bad_accepted = 0
    for row, obsrow in zip(ds, obs):
        d = int.from_bytes(row.tobytes(), "little")
        o = int.from_bytes(obsrow.tobytes(), "little")
        if decoder != "signature" and d not in cache:
            cache[d] = run_decoder(model, d, decoder).observable_mask
        pred = guesses.get(d, 0) if decoder == "signature" else cache[d]
        wrong += bool(o ^ pred)
        wrong_read += bool(public(o ^ pred))
        modular += bool(public(o) ^ public_guesses.get(d, 0))
        if d == 0:
            accepted += 1
            bad_accepted += bool(o)

    def rate(k):
        # Wilson two-sided 95% score interval; k=0 still has finite upper limit.
        z = 1.959963984540054
        f = k / shots
        h = z * (f * (1 - f) / shots + z * z / (4 * shots**2)) ** 0.5
        c = f + z * z / (2 * shots)
        return dict(
            errors=k,
            wilson95=[
                max(0, (c - h) / (1 + z * z / shots)),
                min(1, (c + h) / (1 + z * z / shots)),
            ],
        )

    return dict(
        stim=stim.__version__,
        numpy=np.__version__,
        noise="independent depolarizing gates/resets and measurement flips",
        decoder=decoder,
        probability=str(probability),
        shots=shots,
        seed=seed,
        ideal_terminal_syndrome=inst.boundary,
        detectors=probe.num_detectors,
        observables=probe.num_observables,
        joint=rate(wrong),
        public=rate(wrong_read),
        public_marginal_first_order=rate(modular),
        accepted=accepted,
        damaged_accepted=bad_accepted,
        circuit_sha256=inst.c.sha256(),
    )


def main():
    import argparse, json

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=["residual", "law", "sample"])
    ap.add_argument("spec", type=Path)
    ap.add_argument("--incoming", type=int, default=0)
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--cap", type=int, default=100000)
    ap.add_argument("--p", default="1/1000")
    ap.add_argument("--shots", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--decoder", default="signature")
    ap.add_argument("--export-prefix", type=Path)
    ap.add_argument("-o", type=Path)
    args = ap.parse_args()
    spec = json.loads(args.spec.read_text())
    if args.action == "sample":
        out = sampled(
            spec, args.p, args.shots, args.seed, args.decoder, args.export_prefix
        )
    else:
        inst = Instrument(spec)
        if args.action == "law":
            out = conditional_law(columns(inst)[0], args.p, args.cap)
        else:
            out = residual_profile(inst, args.incoming, args.order, args.cap)
        out["circuit_sha256"] = inst.c.sha256()
    body = json.dumps(out, indent=2) + "\n"
    if args.o:
        args.o.write_text(body)
    else:
        print(body, end="")


if __name__ == "__main__":
    main()
