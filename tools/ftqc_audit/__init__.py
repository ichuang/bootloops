"""Compatibility entry point for ftqc_faults.audit.__init__."""
from pathlib import Path as _Path
import sys as _sys

if __name__.startswith("tools."):
    _tools = str(_Path(__file__).resolve().parent.parent)
    if _tools not in _sys.path:
        _sys.path.insert(0, _tools)

from ftqc_faults.audit import *  # noqa: F401,F403

from ftqc_faults import audit as _engine


def __getattr__(name):
    return getattr(_engine, name)
