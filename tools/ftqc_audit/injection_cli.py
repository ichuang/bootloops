"""Export accepted elementary faults, including decoder residuals when requested."""

import json
from contextlib import ExitStack
from functools import lru_cache


def command(args):
    from ftqc_extractors.injection_circuit import paths, run

    if args.decode and not args.faults:
        raise ValueError("--decode requires --faults")
    predict = None
    if args.decode:
        import stim
        import numpy as np
        from .injection_sampling import build_decoder

        _, path, _, _ = paths(args.case)
        circuit = stim.Circuit.from_file(path)
        decoder = build_decoder(circuit)

        @lru_cache(maxsize=65536)
        def predict(ds):
            syndrome = np.array(
                [(ds >> i) & 1 for i in range(circuit.num_detectors)], dtype=np.uint8
            )
            return sum(int(v) << i for i, v in enumerate(decoder.decode(syndrome)))

    with ExitStack() as ctx:
        stream = ctx.enter_context(open(args.faults, "w")) if args.faults else None

        def write(op, sites, mode, ds, ls, p):
            entry = dict(
                operation=op,
                sites=sites,
                mode=mode,
                decoder_syndrome=ds,
                observable_flip_mask=ls,
                probability=str(p),
            )
            if predict:
                entry["decoded_error_mask"] = ls ^ predict(ds)
            stream.write(json.dumps(entry) + "\n")

        result = run(dict(case=args.case), write if stream else None)
        if result.get("status") != "proved":
            raise ValueError("circuit audit unavailable")
        txt = json.dumps(result, sort_keys=True, indent=2) + "\n"
        if args.out:
            with open(args.out, "w") as dest:
                dest.write(txt)
        else:
            print(txt, end="")
    return 0
