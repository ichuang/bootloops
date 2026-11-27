"""Compatibility entry point for ftqc_faults.zx.rewrites."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("ftqc_faults.zx.rewrites")
