"""Compatibility entry point for ftqc_faults.audit.injection_records."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("ftqc_faults.audit.injection_records")
