"""Optional runtime parity probe for upstream external `kickmix` (Apache-2.0, NOT bundled).

Zero Google code is included or translated here; if an external `kickmix` binary is
installed in the user's Python environment, this adapter cross-checks `.kmx` parse/print
parity at runtime and skips cleanly by name when unavailable or ABI-incompatible.
"""

from __future__ import annotations

from typing import Any

from ..ir_kickmix import KickMixCircuit


def external_kickmix_available() -> tuple[bool, str]:
    try:
        import kickmix  # type: ignore[import-not-found]

        _ = kickmix
        return True, "available"
    except Exception as exc:  # ImportError on missing or ABI-mismatched Cython .so
        return False, f"skipped ({type(exc).__name__}: {exc})"


def verify_parity_with_external_kickmix(circuit: KickMixCircuit) -> dict[str, Any]:
    ok, reason = external_kickmix_available()
    if not ok:
        return {"checked": False, "status": reason}

    import kickmix  # type: ignore[import-not-found]

    kmx_str = circuit.to_kmx(include_directives=False)
    ext_circ = kickmix.Circuit.from_kmx(kmx_str)
    roundtrip = KickMixCircuit.from_kmx(str(ext_circ))
    match = roundtrip.sha256() == circuit.sha256()
    return {
        "checked": True,
        "status": "parity_ok" if match else "mismatch",
        "sha256_match": match,
    }
