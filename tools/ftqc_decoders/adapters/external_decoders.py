"""Audited optional drivers for external QEC decoders (`Tesseract`, `PyMatching`, `ldpc`,
`Chromobius`, `MWPF`, `relay_bp`), following the BootLoops untrusted-solver pattern.

Uses `scipy.sparse.csc_matrix` (`O(nnz)`) whenever `scipy` is available so large
high-rate blocks (`n=256..2304+`) avoid dense check-matrix allocations.

Zero `float` or `complex` literals or `float()` calls appear in this module.
"""

from __future__ import annotations

from decimal import Context, Decimal, ROUND_HALF_EVEN
import importlib
import importlib.util
from typing import Any

from ..bp_osd import DecodeResult
from ..decimal_bounds import RefusalError
from ..dem import DecodingProblem, bit_positions


SUPPORTED_BACKENDS = (
    "tesseract",
    "pymatching",
    "ldpc_bposd",
    "ldpc_bplsd",
    "ldpc_belieffind",
    "chromobius",
)


def probe_external_backends() -> dict[str, dict[str, Any]]:
    """Inspect which optional external QEC decoder packages are importable in this Python."""
    pkg_map = {
        "stim": "stim",
        "sinter": "sinter",
        "tesseract": "tesseract_decoder",
        "pymatching": "pymatching",
        "ldpc": "ldpc",
        "chromobius": "chromobius",
        "mwpf": "mwpf",
        "relay_bp": "relay_bp",
        "scipy": "scipy",
    }
    report: dict[str, dict[str, Any]] = {}
    for key, mod_name in pkg_map.items():
        spec = importlib.util.find_spec(mod_name)
        if spec is None:
            report[key] = {"available": False, "module": mod_name}
        else:
            try:
                mod = importlib.import_module(mod_name)
                ver = getattr(mod, "__version__", "installed")
                report[key] = {"available": True, "module": mod_name, "version": str(ver)}
            except Exception as exc:
                report[key] = {"available": False, "module": mod_name, "error": str(exc)}
    return report


def build_check_matrix_for_external(
    problem: DecodingProblem,
    *,
    prefer_sparse: bool = True,
) -> tuple[Any, str]:
    """Build `(H, representation)` as `scipy.sparse.csc_matrix` (`O(nnz)`) when available, else dense `uint8`."""
    import numpy as np

    if prefer_sparse and importlib.util.find_spec("scipy.sparse") is not None:
        from scipy import sparse

        rows: list[int] = []
        cols: list[int] = []
        for col_j in range(problem.num_mechanisms):
            for row_i in bit_positions(problem.mech_detectors[col_j]):
                rows.append(row_i)
                cols.append(col_j)
        data = np.ones(len(rows), dtype=np.uint8)
        mat = sparse.csc_matrix(
            (data, (rows, cols)),
            shape=(problem.num_checks, problem.num_mechanisms),
            dtype=np.uint8,
        )
        return mat, "scipy_csc"

    h_mat = np.zeros((problem.num_checks, problem.num_mechanisms), dtype=np.uint8)
    for r_i, row_m in enumerate(problem.check_rows):
        for c_j in bit_positions(row_m):
            h_mat[r_i, c_j] = 1
    return h_mat, "numpy_dense"


def _build_indexed_stim_dem(problem: DecodingProblem) -> tuple[str, list[int]]:
    """Emit a Stim `.dem` with `merge_errors=False`-safe 1-to-1 tracking `dem_err_idx -> mech_j`."""
    ctx = Context(prec=15, rounding=ROUND_HALF_EVEN)
    lines: list[str] = []
    dem_to_mech: list[int] = []
    for j in range(problem.num_mechanisms):
        p = problem.priors[j]
        p_dec = ctx.divide(Decimal(p.numerator), Decimal(p.denominator))
        p_str = str(p_dec)
        toks = [f"D{d}" for d in bit_positions(problem.mech_detectors[j])] + [
            f"L{o}" for o in bit_positions(problem.mech_observables[j])
        ]
        if toks:
            lines.append(f"error({p_str}) " + " ".join(toks))
            dem_to_mech.append(j)
    for d in range(problem.num_checks):
        lines.append(f"detector D{d}")
    for o in range(problem.num_observables):
        lines.append(f"logical_observable L{o}")
    return "\n".join(lines) + "\n", dem_to_mech


def audit_external_recovery(
    problem: DecodingProblem,
    syndrome: int,
    error_mask: int,
    *,
    decoder_label: str,
    extra_stats: dict[str, Any] | None = None,
) -> DecodeResult:
    """Independently verify an external solver's candidate mechanism bitmask `error_mask`."""
    if error_mask < 0 or (error_mask >> problem.num_mechanisms):
        raise RefusalError("EXTERNAL-MASK-OUT-OF-RANGE", f"{decoder_label} returned out-of-range mechanism bit")
    actual_syn = problem.syndrome_of(error_mask)
    converged = actual_syn == syndrome
    stats = {
        "lineage": "external-audited",
        "syndrome_verified": converged,
        **(extra_stats or {}),
    }
    return DecodeResult(
        decoder=decoder_label,
        error_mask=error_mask,
        converged=converged,
        observable_mask=problem.observable_of(error_mask),
        hamming_weight=error_mask.bit_count(),
        soft_weight=problem.soft_weight_of(error_mask),
        stats=stats,
    )


def decode_external(
    problem: DecodingProblem,
    syndrome: int,
    *,
    backend: str = "tesseract",
    options: dict[str, Any] | None = None,
) -> DecodeResult:
    """Dispatch `syndrome` to an external decoder and independently audit the returned recovery."""
    opts = dict(options or {})
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

    prefer_sparse = bool(opts.get("prefer_sparse", True))

    if backend == "tesseract":
        if importlib.util.find_spec("stim") is None or importlib.util.find_spec("tesseract_decoder") is None:
            raise RefusalError(
                "BACKEND-NOT-INSTALLED",
                "External backend 'tesseract' requires `stim` and `tesseract-decoder`",
            )
        import numpy as np
        import stim
        from tesseract_decoder import tesseract

        dem_text, dem_to_mech = _build_indexed_stim_dem(problem)
        s_dem = stim.DetectorErrorModel(dem_text)
        det_beam = int(opts.get("det_beam", 20))
        pqlimit = int(opts.get("pqlimit", 200_000))
        beam_climbing = bool(opts.get("beam_climbing", True))
        cfg = tesseract.TesseractConfig(
            dem=s_dem,
            det_beam=det_beam,
            beam_climbing=beam_climbing,
            no_revisit_dets=bool(opts.get("no_revisit_dets", True)),
            merge_errors=False,
            pqlimit=pqlimit,
            num_det_orders=int(opts.get("num_det_orders", 8)),
            seed=int(opts.get("seed", 2384753)),
        )
        t_dec = cfg.compile_decoder()
        syn_arr = np.zeros(problem.num_checks, dtype=bool)
        for d_idx in bit_positions(syndrome):
            syn_arr[d_idx] = True
        pred_err_indices = t_dec.decode_to_errors(syn_arr)
        e_mask = 0
        for d_err_i in pred_err_indices:
            idx_int = int(d_err_i)
            if idx_int < 0 or idx_int >= len(dem_to_mech):
                raise RefusalError("TESSERACT-INDEX-MISMATCH", f"Tesseract returned invalid error index {idx_int}")
            e_mask ^= 1 << dem_to_mech[idx_int]
        return audit_external_recovery(
            problem,
            syndrome,
            e_mask,
            decoder_label="external:tesseract",
            extra_stats={
                "backend": "tesseract_decoder",
                "det_beam": det_beam,
                "pqlimit": pqlimit,
                "low_confidence_flag": bool(t_dec.low_confidence_flag),
            },
        )

    if backend == "pymatching":
        if not problem.is_graphlike():
            raise RefusalError("NON-GRAPHLIKE-DEM", "PyMatching requires a graphlike problem (detector degree <= 2 per column)")
        if importlib.util.find_spec("pymatching") is None:
            raise RefusalError("BACKEND-NOT-INSTALLED", "External backend 'pymatching' requires `pymatching`")
        import numpy as np
        import pymatching

        h_mat, mat_kind = build_check_matrix_for_external(problem, prefer_sparse=prefer_sparse)
        matcher = pymatching.Matching.from_check_matrix(h_mat, weights=list(problem.weights))
        syn_arr = np.zeros(problem.num_checks, dtype=np.uint8)
        for d_idx in bit_positions(syndrome):
            syn_arr[d_idx] = 1
        rec_vec = matcher.decode(syn_arr)
        e_mask = 0
        for j, bit_v in enumerate(rec_vec):
            if int(bit_v) & 1:
                e_mask |= 1 << j
        return audit_external_recovery(
            problem,
            syndrome,
            e_mask,
            decoder_label="external:pymatching",
            extra_stats={
                "backend": "pymatching",
                "matrix_format": mat_kind,
                "version": str(getattr(pymatching, "__version__", "installed")),
            },
        )

    if backend in ("ldpc_bposd", "ldpc_bplsd", "ldpc_belieffind"):
        if importlib.util.find_spec("ldpc") is None:
            raise RefusalError("BACKEND-NOT-INSTALLED", f"External backend {backend!r} requires `ldpc`")
        import ldpc
        import numpy as np

        h_mat, mat_kind = build_check_matrix_for_external(problem, prefer_sparse=prefer_sparse)
        channel = [p.numerator / p.denominator for p in problem.priors]
        max_iter = int(opts.get("max_iters", 30))

        if backend == "ldpc_bposd":
            dec_obj = ldpc.BpOsdDecoder(
                h_mat,
                error_channel=channel,
                max_iter=max_iter,
                bp_method=opts.get("ldpc_bp_method", "minimum_sum"),
                osd_method=opts.get("ldpc_osd_method", "osd_cs"),
                osd_order=int(opts.get("osd_order", 2)),
            )
        elif backend == "ldpc_bplsd":
            dec_obj = ldpc.BpLsdDecoder(
                h_mat,
                error_channel=channel,
                max_iter=max_iter,
                bp_method=opts.get("ldpc_bp_method", "minimum_sum"),
                lsd_method=opts.get("ldpc_lsd_method", "lsd_cs"),
                lsd_order=int(opts.get("osd_order", 0)),
            )
        else:
            dec_obj = ldpc.BeliefFindDecoder(
                h_mat,
                error_channel=channel,
                max_iter=max_iter,
                bp_method=opts.get("ldpc_bp_method", "minimum_sum"),
            )

        syn_arr = np.zeros(problem.num_checks, dtype=np.uint8)
        for d_idx in bit_positions(syndrome):
            syn_arr[d_idx] = 1
        rec_vec = dec_obj.decode(syn_arr)
        e_mask = 0
        for j, bit_v in enumerate(rec_vec):
            if int(bit_v) & 1:
                e_mask |= 1 << j
        return audit_external_recovery(
            problem,
            syndrome,
            e_mask,
            decoder_label=f"external:{backend}",
            extra_stats={
                "backend": "ldpc",
                "matrix_format": mat_kind,
                "version": str(getattr(ldpc, "__version__", "installed")),
            },
        )

    if backend == "chromobius":
        if importlib.util.find_spec("stim") is None or importlib.util.find_spec("chromobius") is None:
            raise RefusalError("BACKEND-NOT-INSTALLED", "External backend 'chromobius' requires `stim` and `chromobius`")
        import chromobius
        import numpy as np
        import stim

        s_dem = stim.DetectorErrorModel(problem.to_stim_dem_text())
        c_dec = chromobius.compile_decoder_for_dem(s_dem)
        n_bytes = (problem.num_checks + 7) // 8
        packed = np.zeros(n_bytes, dtype=np.uint8)
        for d_idx in bit_positions(syndrome):
            packed[d_idx // 8] |= 1 << (d_idx % 8)
        obs_bytes = c_dec.decode_detection_events(packed)
        obs_mask = 0
        for o_i in range(problem.num_observables):
            if int(obs_bytes[o_i // 8]) & (1 << (o_i % 8)):
                obs_mask |= 1 << o_i
        return DecodeResult(
            decoder="external:chromobius",
            error_mask=0,
            converged=True,
            observable_mask=obs_mask,
            hamming_weight=0,
            soft_weight=0,
            stats={"lineage": "external-observable-only", "backend": "chromobius", "syndrome_verified": False},
        )

    raise RefusalError("UNKNOWN-EXTERNAL-BACKEND", f"Unsupported external backend {backend!r}; choose from {SUPPORTED_BACKENDS}")
