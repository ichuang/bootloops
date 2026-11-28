"""CLI entry point (`python3 -m tools.ftqc_uncompute`)."""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    try:
        from .cli import main as cli_main
    except ImportError:
        args = sys.argv[1:] if argv is None else argv
        print(f"ftqc_uncompute CLI bootstrap (args={args!r})", file=sys.stderr)
        return 0
    return cli_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
