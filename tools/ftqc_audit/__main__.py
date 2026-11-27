"""Compatibility entry point for ftqc_faults.audit.__main__."""
import sys
from importlib import import_module

if __name__ == "__main__":
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    raise SystemExit(import_module("ftqc_faults.audit.__main__").main())

sys.modules[__name__] = import_module("ftqc_faults.audit.__main__")
