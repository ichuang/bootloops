"""Circuit fault analysis, ZX equivalence, and the Steane cross-check."""

import argparse
import json
import sys
from pathlib import Path


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(prog="ftqc_faults", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    bridge = commands.add_parser(
        "steane", help="check both five-ancilla Steane extractors"
    )
    output = bridge.add_mutually_exclusive_group()
    output.add_argument("-o", type=Path, help="save the cross-check receipt")
    output.add_argument("--verify", type=Path, help="replay a saved receipt")
    commands.add_parser("audit", help="analyze a physical Clifford schedule")
    commands.add_parser("zx", help="synthesize or verify ZX fault equivalence")
    if argv[:1] == ["audit"]:
        from .audit.__main__ import main as audit_main

        return audit_main(argv[1:], prog="ftqc_faults audit")
    if argv[:1] == ["zx"]:
        from .zx.cli import main as zx_main

        return zx_main(argv[1:], prog="ftqc_faults zx")
    args = parser.parse_args(argv)
    from .steane import receipt

    try:
        result = receipt()
        if args.verify:
            saved = json.loads(args.verify.read_text())
            if json.dumps(saved, sort_keys=True) != json.dumps(result, sort_keys=True):
                raise ValueError("Steane receipt differs from replay")
            print("Verified both Steane extractor sectors")
        elif args.o:
            args.o.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        else:
            print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
