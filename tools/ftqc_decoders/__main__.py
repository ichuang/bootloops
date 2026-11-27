"""Dual-mode CLI entrypoint supporting both `PYTHONPATH=tools python3 -m ftqc_decoders`
and `python3 -m tools.ftqc_decoders` from the BootLoops repository root.
"""

from __future__ import annotations

from pathlib import Path
import sys

# Ensure `tools/` is on `sys.path` regardless of how the module was invoked
_tools_dir = str(Path(__file__).resolve().parents[1])
if _tools_dir not in sys.path:
    sys.path.insert(0, _tools_dir)

from .cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
