"""Compatibility shim: `qeccert` was renamed to `ftqc_codes`."""
import importlib as _importlib
import sys as _sys

_pkg = _importlib.import_module("ftqc_codes")
for _m in (
    "gf2",
    "codes",
    "minweight",
    "certify",
    "cnf",
    "lrat",
    "sat",
    "bch",
    "cuts",
    "concat",
    "leangen",
):
    _sub = _importlib.import_module(f"ftqc_codes.{_m}")
    setattr(_pkg, _m, _sub)
    _sys.modules[f"{__name__}.{_m}"] = _sub
_sys.modules[__name__] = _pkg
