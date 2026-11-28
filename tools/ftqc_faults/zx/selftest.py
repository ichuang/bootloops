#!/usr/bin/env python3
"""Self-contained BootLoops acceptance battery entrypoint for `ftqc_faults.zx` (`selftest` class).

Runs all positive legs (L1-L5), planted negative/refusal controls (N1-N5), and all 4
scratch-copy mutation-sabotage controls (M1-M4) from an isolated temporary directory.

Zero float literals or float builtin references appear anywhere in this module.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile

PKG_DIR = Path(__file__).resolve().parent
TOOLS_PARENT = PKG_DIR.parents[1]
if str(TOOLS_PARENT) not in sys.path:
    sys.path.insert(0, str(TOOLS_PARENT))

from ftqc_faults.zx.battery import run_battery  # noqa: E402


def main() -> int:
    orig_cwd = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as scratch_cwd:
            os.chdir(scratch_cwd)
            report = run_battery(include_mutations=True)
    finally:
        os.chdir(orig_cwd)

    for leg in report["legs"]:
        tag = "PASS" if leg["passed"] else "FAIL"
        print(f"  [{tag}] {leg['id']} ({leg['elapsed_ms']} ms): {leg['summary']}")

    overall = "PASS" if report["passed"] else "FAIL"
    print(f"[{overall}] ftqc_faults.zx selftest ({report['passed_count']}/{report['total_count']} legs)")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
