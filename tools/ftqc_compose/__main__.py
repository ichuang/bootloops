"""Generate and replay content-bound deterministic composition receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

from .supports import family_result


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def receipt(spec):
    if spec["kind"] == "family":
        result = family_result(spec)
    else:
        raise ValueError("unknown study kind")
    return dict(schema="ftqc_compose-1", input_sha256=digest(spec), result=result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check")
    check.add_argument("spec", type=Path)
    check.add_argument("-o", type=Path)
    verify = commands.add_parser("verify")
    verify.add_argument("spec", type=Path)
    verify.add_argument("receipt", type=Path)
    args = parser.parse_args()
    try:
        result = receipt(json.loads(args.spec.read_text()))
        if args.command == "verify":
            if canonical(result) != canonical(json.loads(args.receipt.read_text())):
                raise ValueError("receipt differs from recomputed evidence")
            print("verified", result["result"]["status"])
            return 0
        content = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.o:
            args.o.write_text(content)
        else:
            print(content, end="")
        return int(result["result"]["status"] == "counterexample")
    except (ValueError, TypeError, KeyError, AttributeError, OSError) as error:
        print(f"ftqc_compose: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
