#!/usr/bin/env python3
"""Direct script entrypoint (`python3 tools/faultweb/faultweb.py`)."""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from faultweb.cli import main
else:
    from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
