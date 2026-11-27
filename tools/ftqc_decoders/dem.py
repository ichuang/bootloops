"""Canonical parity-check & Detector Error Model (DEM) IR for `ftqc_decoders`.

Provides a single exact, zero-float representation (`DecodingProblem`) unifying:
  1. Static CSS codes from `ftqc_codes` (`CSSCode`, family specs, `ftqc_codes-receipt-1`),
     carrying check rows `H`, logical observables `L`, and dual stabilizer generators `G`
     (used both for coset-equivalence verification and quantum split-stabilizer trapping
     set discovery: Poulin & Chung, QIC 8, 987, 2008; Raveendran & Vasic, Quantum 5, 562, 2021).
  2. Circuit-level spacetime fault models from `ftqc_audit.RecordModel` and
     `ftqc_extractors` (`gross_circuit`, `reading_audit`, HGP schedules, and
     `surgery_sequence` mapping cones with metachecks: Cowtan, He, Williamson & Yoder,
     arXiv:2510.14895; Campbell, QST 4, 025006, 2019).
  3. Universal Stim `.dem` (`DetectorErrorModel`) text files (Gidney, Quantum 5, 497, 2021),
     parsed and emitted in pure Python stdlib with exact `fractions.Fraction` probabilities
     (never `float`) and lossless `^` graphlike decomposition legs for interoperability
     with `Tesseract` (Beni et al., arXiv:2503.10988), `PyMatching` v2 (Higgott & Gidney,
     arXiv:2303.15933), `ldpc` (Roffe et al., arXiv:2005.07016; Hillmann et al.,
     arXiv:2406.18655), `Chromobius` (arXiv:2312.08813), and `MWPF` (arXiv:2508.04969).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any

from .decimal_bounds import (
    RefusalError,
    parse_exact_probability,
    scaled_log_odds_int,
    xor_prob,
)


def _ensure_tools_on_path() -> None:
    tools_dir = str(Path(__file__).resolve().parents[1])
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)


def bit_positions(mask: int) -> list[int]:
    """Return sorted 0-based bit indices set in nonnegative integer `mask`."""
    out: list[int] = []
    cur = mask
    while cur:
        lsb = cur & -cur
        out.append(lsb.bit_length() - 1)
        cur ^= lsb
    return out


def transpose_bitmasks(rows: tuple[int, ...] | list[int], ncols: int) -> tuple[int, ...]:
    """Transpose `m` bitmask rows over `0..ncols-1` into `ncols` column bitmasks over `0..m-1`."""
    cols = [0] * ncols
    for r_idx, row in enumerate(rows):
        cur = row
        r_bit = 1 << r_idx
        while cur:
            lsb = cur & -cur
            c_idx = lsb.bit_length() - 1
            if c_idx >= ncols:
                raise RefusalError("MATRIX-COLUMN-OUT-OF-RANGE", f"Bit {c_idx} >= ncols={ncols}")
            cols[c_idx] ^= r_bit
            cur ^= lsb
    return tuple(cols)


@dataclass(frozen=True)
class DecodingProblem:
    """Exact binary parity/detector decoding instance `(H, L, priors, G_stab, M_meta)`.

    Given a syndrome bitmask `syndrome` in `GF(2)^num_checks` (and optional noisy
    syndrome vector repaired via `metacheck_rows`), a valid recovery mask `e_hat`
    over `0..num_mechanisms-1` satisfies `syndrome_of(e_hat) == syndrome`, and
    succeeds logically on true fault `e_true` iff `observable_of(e_true ^ e_hat) == 0`.
    """

    name: str
    num_checks: int
    num_mechanisms: int
    num_observables: int
    check_rows: tuple[int, ...]
    mech_detectors: tuple[int, ...]
    observable_rows: tuple[int, ...]
    mech_observables: tuple[int, ...]
    priors: tuple[Fraction, ...]
    weights: tuple[int, ...]
    stabilizer_generators: tuple[int, ...] = ()
    metacheck_rows: tuple[int, ...] = ()
    decomposed_legs: tuple[tuple[tuple[int, int], ...], ...] | None = None
    witness_map: dict[int, tuple[int, str]] = field(default_factory=dict)
    source_meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.num_checks < 0 or self.num_mechanisms <= 0 or self.num_observables < 0:
            raise RefusalError(
                "INVALID-PROBLEM-DIMENSIONS",
                f"Invalid dimensions checks={self.num_checks}, mechs={self.num_mechanisms}, obs={self.num_observables}",
            )
        if len(self.check_rows) != self.num_checks:
            raise RefusalError("CHECK-ROW-MISMATCH", "len(check_rows) != num_checks")
        if len(self.mech_detectors) != self.num_mechanisms:
            raise RefusalError("MECH-DET-MISMATCH", "len(mech_detectors) != num_mechanisms")
        if len(self.observable_rows) != self.num_observables:
            raise RefusalError("OBS-ROW-MISMATCH", "len(observable_rows) != num_observables")
        if len(self.mech_observables) != self.num_mechanisms:
            raise RefusalError("MECH-OBS-MISMATCH", "len(mech_observables) != num_mechanisms")
        if len(self.priors) != self.num_mechanisms or len(self.weights) != self.num_mechanisms:
            raise RefusalError("PRIOR-LEN-MISMATCH", "priors/weights length must equal num_mechanisms")

        for g in self.stabilizer_generators:
            if g < 0 or (g >> self.num_mechanisms):
                raise RefusalError("STAB-OUT-OF-RANGE", "Stabilizer generator support exceeds num_mechanisms")
            if self.syndrome_of(g) != 0 or self.observable_of(g) != 0:
                raise RefusalError(
                    "INVALID-STABILIZER-GENERATOR",
                    "Declared stabilizer generator either triggers a check or flips a logical observable",
                )

        for m_row in self.metacheck_rows:
            if m_row < 0 or (m_row >> self.num_checks):
                raise RefusalError("METACHECK-OUT-OF-RANGE", "Metacheck row support exceeds num_checks")
            for det_col in self.mech_detectors:
                if (m_row & det_col).bit_count() & 1:
                    raise RefusalError(
                        "INVALID-METACHECK",
                        "Metacheck row does not annihilate column space of H (M * H != 0 mod 2)",
                    )

    @classmethod
    def from_matrices(
        cls,
        check_rows: list[int] | tuple[int, ...],
        observable_rows: list[int] | tuple[int, ...],
        num_mechanisms: int,
        *,
        name: str = "parity_problem",
        priors: list[str | int | Fraction] | tuple[str | int | Fraction, ...] | None = None,
        default_prior: str | Fraction = "1/1000",
        scale: int = 1024,
        stabilizer_generators: list[int] | tuple[int, ...] = (),
        metacheck_rows: list[int] | tuple[int, ...] = (),
        decomposed_legs: tuple[tuple[tuple[int, int], ...], ...] | None = None,
        witness_map: dict[int, tuple[int, str]] | None = None,
        source_meta: dict[str, Any] | None = None,
    ) -> DecodingProblem:
        """Construct a `DecodingProblem` from row bitmasks `H` and `L`."""
        c_rows = tuple(int(r) for r in check_rows)
        o_rows = tuple(int(r) for r in observable_rows)
        mech_det = transpose_bitmasks(c_rows, num_mechanisms)
        mech_obs = transpose_bitmasks(o_rows, num_mechanisms)
        if priors is None:
            p0 = parse_exact_probability(default_prior, allow_zero=False, allow_one=False, max_half=True)
            p_tup = (p0,) * num_mechanisms
        else:
            p_tup = tuple(
                parse_exact_probability(p, allow_zero=False, allow_one=False, max_half=True)
                for p in priors
            )
        w_tup = tuple(scaled_log_odds_int(p, scale=scale) for p in p_tup)
        return cls(
            name=name,
            num_checks=len(c_rows),
            num_mechanisms=num_mechanisms,
            num_observables=len(o_rows),
            check_rows=c_rows,
            mech_detectors=mech_det,
            observable_rows=o_rows,
            mech_observables=mech_obs,
            priors=p_tup,
            weights=w_tup,
            stabilizer_generators=tuple(int(g) for g in stabilizer_generators),
            metacheck_rows=tuple(int(m) for m in metacheck_rows),
            decomposed_legs=decomposed_legs,
            witness_map=dict(witness_map or {}),
            source_meta=dict(source_meta or {}),
        )

    def syndrome_of(self, error_mask: int) -> int:
        """Compute `H * e (mod 2)` as a check bitmask over `0..num_checks-1`."""
        syn = 0
        cur = error_mask
        while cur:
            lsb = cur & -cur
            j = lsb.bit_length() - 1
            if j >= self.num_mechanisms:
                raise RefusalError("ERROR-MASK-OUT-OF-RANGE", f"Mechanism index {j} >= {self.num_mechanisms}")
            syn ^= self.mech_detectors[j]
            cur ^= lsb
        return syn

    def observable_of(self, error_mask: int) -> int:
        """Compute `L * e (mod 2)` as an observable bitmask over `0..num_observables-1`."""
        obs = 0
        cur = error_mask
        while cur:
            lsb = cur & -cur
            j = lsb.bit_length() - 1
            if j >= self.num_mechanisms:
                raise RefusalError("ERROR-MASK-OUT-OF-RANGE", f"Mechanism index {j} >= {self.num_mechanisms}")
            obs ^= self.mech_observables[j]
            cur ^= lsb
        return obs

    def metacheck_syndrome_of(self, noisy_syndrome: int) -> int:
        """Compute `M * s_noisy (mod 2)` over `0..len(metacheck_rows)-1`."""
        ms = 0
        for i, m_row in enumerate(self.metacheck_rows):
            if (m_row & noisy_syndrome).bit_count() & 1:
                ms |= 1 << i
        return ms

    def soft_weight_of(self, error_mask: int) -> int:
        """Exact integer log-odds cost `sum_{j in supp(error_mask)} weights[j]`."""
        total = 0
        cur = error_mask
        while cur:
            lsb = cur & -cur
            total += self.weights[lsb.bit_length() - 1]
            cur ^= lsb
        return total

    def is_graphlike(self) -> bool:
        """True iff every mechanism flips at most 2 detectors directly."""
        return all(d.bit_count() <= 2 for d in self.mech_detectors)

    def has_graphlike_decomposition(self) -> bool:
        """True iff `is_graphlike()` or `decomposed_legs` decomposes all mechanisms into <=2-detector components."""
        if self.is_graphlike():
            return True
        if self.decomposed_legs is None:
            return False
        return all(
            all(det_m.bit_count() <= 2 for det_m, _ in legs)
            for legs in self.decomposed_legs
        )

    def max_check_degree(self) -> int:
        return max((r.bit_count() for r in self.check_rows), default=0)

    def max_column_degree(self) -> int:
        return max((c.bit_count() for c in self.mech_detectors), default=0)

    def max_interaction_degree(self) -> int:
        """Maximum syndrome-graph neighbor degree among fault mechanisms (for Peierls cluster bounds)."""
        check_members = [bit_positions(r) for r in self.check_rows]
        best = 0
        for j in range(self.num_mechanisms):
            nbrs: set[int] = set()
            for c_idx in bit_positions(self.mech_detectors[j]):
                nbrs.update(check_members[c_idx])
            nbrs.discard(j)
            if len(nbrs) > best:
                best = len(nbrs)
        return best

    def canonical_sha256(self) -> str:
        payload = {
            "num_checks": self.num_checks,
            "num_mechanisms": self.num_mechanisms,
            "num_observables": self.num_observables,
            "check_rows": [str(r) for r in self.check_rows],
            "observable_rows": [str(r) for r in self.observable_rows],
            "priors": [f"{p.numerator}/{p.denominator}" for p in self.priors],
            "stabilizer_generators": [str(g) for g in self.stabilizer_generators],
            "metacheck_rows": [str(m) for m in self.metacheck_rows],
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def to_stim_dem_text(self, *, decimal_digits: int = 15) -> str:
        """Serialize to Stim `DetectorErrorModel` (`.dem`) text format."""
        from decimal import Context, Decimal, ROUND_HALF_EVEN

        ctx = Context(prec=decimal_digits, rounding=ROUND_HALF_EVEN)
        lines: list[str] = [f"# ftqc_decoders DEM: {self.name} (sha256={self.canonical_sha256()[:16]})"]
        for j in range(self.num_mechanisms):
            p = self.priors[j]
            p_dec = ctx.divide(Decimal(p.numerator), Decimal(p.denominator))
            p_str = format(p_dec, "f").rstrip("0").rstrip(".") if "E" not in str(p_dec) and "e" not in str(p_dec) else str(p_dec)
            if not p_str or p_str == "0":
                p_str = str(p_dec)

            if self.decomposed_legs is not None and len(self.decomposed_legs[j]) > 1:
                comp_strs: list[str] = []
                for det_m, obs_m in self.decomposed_legs[j]:
                    toks = [f"D{d}" for d in bit_positions(det_m)] + [f"L{o}" for o in bit_positions(obs_m)]
                    if toks:
                        comp_strs.append(" ".join(toks))
                target_str = " ^ ".join(comp_strs)
            else:
                toks = [f"D{d}" for d in bit_positions(self.mech_detectors[j])] + [
                    f"L{o}" for o in bit_positions(self.mech_observables[j])
                ]
                target_str = " ".join(toks)

            if target_str:
                lines.append(f"error({p_str}) {target_str}")
        for d in range(self.num_checks):
            lines.append(f"detector D{d}")
        for o in range(self.num_observables):
            lines.append(f"logical_observable L{o}")
        return "\n".join(lines) + "\n"

    @classmethod
    def from_stim_dem_text(
        cls,
        text: str,
        *,
        name: str = "stim_dem",
        merge_parallel: bool = True,
        scale: int = 1024,
    ) -> DecodingProblem:
        """Parse a Stim `.dem` string with zero floats."""
        raw_lines = [ln.split("#", 1)[0].strip() for ln in text.splitlines()]
        tokens_lines = [ln for ln in raw_lines if ln]

        def parse_block(idx: int, det_offset: int) -> tuple[int, int, list[tuple[Fraction, int, int, tuple[tuple[int, int], ...]]], int, int]:
            mechs: list[tuple[Fraction, int, int, tuple[tuple[int, int], ...]]] = []
            max_det = -1
            max_obs = -1
            cur_offset = det_offset
            i = idx
            while i < len(tokens_lines):
                ln = tokens_lines[i]
                if ln == "}":
                    return i + 1, cur_offset, mechs, max_det, max_obs
                rep_m = re.match(r"^repeat\s+(\d+)\s*\{\s*$", ln)
                if rep_m:
                    reps = int(rep_m.group(1))
                    end_i = i + 1
                    for _ in range(reps):
                        end_i, cur_offset, sub_m, sub_d, sub_o = parse_block(i + 1, cur_offset)
                        mechs.extend(sub_m)
                        max_det = max(max_det, sub_d)
                        max_obs = max(max_obs, sub_o)
                    i = end_i
                    continue

                err_m = re.match(r"^error\s*\(\s*([^)]+)\s*\)\s*(.*)$", ln)
                if err_m:
                    p = parse_exact_probability(err_m.group(1), allow_zero=True, allow_one=False, max_half=True)
                    rest = err_m.group(2).strip()
                    if p > 0 and rest:
                        raw_comps = rest.split("^")
                        legs: list[tuple[int, int]] = []
                        tot_det = 0
                        tot_obs = 0
                        for comp in raw_comps:
                            c_det = 0
                            c_obs = 0
                            for tok in comp.split():
                                if tok.startswith("D"):
                                    d_id = int(tok[1:]) + cur_offset
                                    c_det ^= 1 << d_id
                                    tot_det ^= 1 << d_id
                                    max_det = max(max_det, d_id)
                                elif tok.startswith("L"):
                                    o_id = int(tok[1:])
                                    c_obs ^= 1 << o_id
                                    tot_obs ^= 1 << o_id
                                    max_obs = max(max_obs, o_id)
                                else:
                                    raise RefusalError("INVALID-DEM-TARGET", f"Unsupported DEM target token {tok!r}")
                            if c_det or c_obs:
                                legs.append((c_det, c_obs))
                        if tot_det or tot_obs:
                            mechs.append((p, tot_det, tot_obs, tuple(legs)))
                    i += 1
                    continue

                if ln.startswith("shift_detectors"):
                    after = re.sub(r"^shift_detectors(\s*\([^)]*\))?\s*", "", ln).strip()
                    if after:
                        cur_offset += int(after.split()[0])
                    i += 1
                    continue

                if ln.startswith("detector"):
                    after = re.sub(r"^detector(\s*\([^)]*\))?\s*", "", ln).strip()
                    for tok in after.split():
                        if tok.startswith("D"):
                            max_det = max(max_det, int(tok[1:]) + cur_offset)
                    i += 1
                    continue

                if ln.startswith("logical_observable"):
                    after = ln[len("logical_observable") :].strip()
                    for tok in after.split():
                        if tok.startswith("L"):
                            max_obs = max(max_obs, int(tok[1:]))
                    i += 1
                    continue

                raise RefusalError("UNSUPPORTED-DEM-INSTRUCTION", f"Cannot parse DEM line: {ln!r}")

            return i, cur_offset, mechs, max_det, max_obs

        _, _, raw_mechs, max_det, max_obs = parse_block(0, 0)
        if not raw_mechs:
            raise RefusalError("EMPTY-DEM", "DEM contains no nontrivial error mechanisms")

        num_checks = max_det + 1 if max_det >= 0 else 0
        num_obs = max_obs + 1 if max_obs >= 0 else 0

        if merge_parallel:
            merged: dict[tuple[int, int], tuple[Fraction, tuple[tuple[int, int], ...]]] = {}
            for p, d_m, o_m, legs in raw_mechs:
                key = (d_m, o_m)
                if key in merged:
                    prev_p, prev_legs = merged[key]
                    best_legs = prev_legs if all(dm.bit_count() <= 2 for dm, _ in prev_legs) else legs
                    merged[key] = (xor_prob(prev_p, p), best_legs)
                else:
                    merged[key] = (p, legs)
            ordered = [(d_m, o_m, p, legs) for (d_m, o_m), (p, legs) in merged.items()]
        else:
            ordered = [(d_m, o_m, p, legs) for p, d_m, o_m, legs in raw_mechs]

        num_mechs = len(ordered)
        mech_det = tuple(item[0] for item in ordered)
        mech_obs = tuple(item[1] for item in ordered)
        priors_list = []
        for item in ordered:
            pv = item[2]
            if pv > Fraction(1, 2):
                pv = Fraction(1, 2)
            priors_list.append(pv)
        priors_t = tuple(priors_list)
        weights_t = tuple(scaled_log_odds_int(p, scale=scale) for p in priors_t)
        check_rows = transpose_bitmasks(mech_det, num_checks) if num_checks > 0 else ()
        obs_rows = transpose_bitmasks(mech_obs, num_obs) if num_obs > 0 else ()
        decomp = tuple(item[3] for item in ordered)
        return cls(
            name=name,
            num_checks=num_checks,
            num_mechanisms=num_mechs,
            num_observables=num_obs,
            check_rows=check_rows,
            mech_detectors=mech_det,
            observable_rows=obs_rows,
            mech_observables=mech_obs,
            priors=priors_t,
            weights=weights_t,
            decomposed_legs=decomp,
            source_meta={"format": "stim_dem", "merged_parallel": merge_parallel},
        )

    @classmethod
    def from_css_code(
        cls,
        code_or_spec: Any,
        *,
        sector: str = "X",
        prior: str | Fraction = "1/1000",
        scale: int = 1024,
    ) -> DecodingProblem:
        """Build a single-sector CSS code `DecodingProblem` from an `ftqc_codes.CSSCode`,
        family name, or receipt dict, preserving `tanner` local neighborhoods when present.
        """
        _ensure_tools_on_path()
        from ftqc_codes import codes as fc_codes, gf2 as fc_gf2

        if isinstance(code_or_spec, str):
            code = fc_codes.load(code_or_spec)
        elif isinstance(code_or_spec, dict):
            if "code" in code_or_spec and isinstance(code_or_spec["code"], dict):
                code = fc_codes.CSSCode.from_json(code_or_spec["code"])
            else:
                code = fc_codes.CSSCode.from_json(code_or_spec)
        else:
            code = code_or_spec

        if not code.commutes():
            raise RefusalError("NON-COMMUTING-CSS", f"CSS checks of {code.name!r} do not commute (H_X * H_Z^T != 0)")

        n = int(code.n)
        if sector == "X":
            checks = list(code.hz)
            stabs = list(code.hx)
        elif sector == "Z":
            checks = list(code.hx)
            stabs = list(code.hz)
        else:
            raise RefusalError("INVALID-CSS-SECTOR", f"sector must be 'X' or 'Z', got {sector!r}")

        ech_checks = fc_gf2.Echelon(checks)
        logical_duals = [v for v in fc_gf2.nullspace(stabs, n) if ech_checks.add(v)]
        smeta: dict[str, Any] = {
            "kind": "css_code",
            "code_name": code.name,
            "code_sha256": code.sha256(),
            "sector": sector,
            "n": n,
            "k": len(logical_duals),
        }
        if isinstance(getattr(code, "meta", None), dict) and "tanner" in code.meta:
            smeta["tanner"] = code.meta["tanner"]

        return cls.from_matrices(
            checks,
            logical_duals,
            n,
            name=f"{code.name}:{sector}",
            default_prior=prior,
            scale=scale,
            stabilizer_generators=stabs,
            source_meta=smeta,
        )

    @classmethod
    def from_record_model(
        cls,
        record_model: Any,
        *,
        name: str | None = None,
        location_prior: str | Fraction = "1/1000",
        scale: int = 1024,
    ) -> DecodingProblem:
        """Convert an `ftqc_audit.RecordModel` into an exact spacetime `DecodingProblem`."""
        p_loc = parse_exact_probability(location_prior, allow_zero=False, allow_one=False, max_half=True)
        nd = int(record_model.nd)
        det_mask_all = (1 << nd) - 1
        n_obs = len(record_model.circuit.observables)
        obs_mask_all = (1 << n_obs) - 1

        merged: dict[tuple[int, int], Fraction] = {}
        first_witness: dict[tuple[int, int], tuple[int, str]] = {}

        for loc_idx, group in enumerate(record_model.groups):
            if not group:
                continue
            mode_p = p_loc / len(group)
            for label in sorted(group):
                sig = int(group[label])
                d_m = sig & det_mask_all
                o_m = (sig >> nd) & obs_mask_all
                if d_m == 0 and o_m == 0:
                    continue
                key = (d_m, o_m)
                if key in merged:
                    merged[key] = xor_prob(merged[key], mode_p)
                else:
                    merged[key] = mode_p
                    first_witness[key] = (loc_idx, label)

        if not merged:
            raise RefusalError("TRIVIAL-RECORD-MODEL", "RecordModel has zero nontrivial detector/observable mechanisms")

        ordered_keys = sorted(merged.keys())
        num_mechs = len(ordered_keys)
        mech_det = tuple(k[0] for k in ordered_keys)
        mech_obs = tuple(k[1] for k in ordered_keys)
        priors_t = tuple(min(Fraction(1, 2), merged[k]) for k in ordered_keys)
        weights_t = tuple(scaled_log_odds_int(p, scale=scale) for p in priors_t)
        check_rows = transpose_bitmasks(mech_det, nd) if nd > 0 else ()
        obs_rows = transpose_bitmasks(mech_obs, n_obs) if n_obs > 0 else ()
        w_map = {j: first_witness[k] for j, k in enumerate(ordered_keys)}

        circ_name = getattr(record_model.circuit, "name", None) or "record_model"
        return cls(
            name=name or circ_name,
            num_checks=nd,
            num_mechanisms=num_mechs,
            num_observables=n_obs,
            check_rows=check_rows,
            mech_detectors=mech_det,
            observable_rows=obs_rows,
            mech_observables=mech_obs,
            priors=priors_t,
            weights=weights_t,
            witness_map=w_map,
            source_meta={
                "kind": "ftqc_audit_record_model",
                "circuit_sha256": record_model.circuit.sha256(),
                "locations": len(record_model.groups),
                "deduplicated_mechanisms": num_mechs,
            },
        )


def repetition_problem(distance: int, *, prior: str | Fraction = "1/100") -> DecodingProblem:
    """1D repetition code `[d, 1, d]` with open boundaries (graphlike, column weight <= 2)."""
    if distance < 3:
        raise RefusalError("INVALID-DISTANCE", "Repetition code distance must be >= 3")
    checks = [(1 << i) | (1 << (i + 1)) for i in range(distance - 1)]
    obs = [(1 << distance) - 1]
    return DecodingProblem.from_matrices(
        checks, obs, distance, name=f"repetition_d{distance}", default_prior=prior
    )


def toric_x_problem(length: int, *, prior: str | Fraction = "1/100") -> DecodingProblem:
    """2D `[[2 L^2, 2, L]]` Kitaev toric code X-error decoding instance."""
    if length < 2:
        raise RefusalError("INVALID-LATTICE-SIZE", "Toric lattice length L must be >= 2")
    L = length

    def h_edge(x: int, y: int) -> int:
        return (y % L) * L + (x % L)

    def v_edge(x: int, y: int) -> int:
        return L * L + (y % L) * L + (x % L)

    hz_rows: list[int] = []
    for y in range(L):
        for x in range(L):
            row = (1 << h_edge(x, y)) | (1 << h_edge(x - 1, y)) | (1 << v_edge(x, y)) | (1 << v_edge(x, y - 1))
            hz_rows.append(row)

    hx_rows: list[int] = []
    for y in range(L):
        for x in range(L):
            plaq = (1 << h_edge(x, y)) | (1 << h_edge(x, y + 1)) | (1 << v_edge(x, y)) | (1 << v_edge(x + 1, y))
            hx_rows.append(plaq)

    lz1 = 0
    lz2 = 0
    for k in range(L):
        lz1 |= 1 << h_edge(0, k)
        lz2 |= 1 << v_edge(k, 0)

    return DecodingProblem.from_matrices(
        hz_rows,
        [lz1, lz2],
        2 * L * L,
        name=f"toric_{L}x{L}:X",
        default_prior=prior,
        stabilizer_generators=hx_rows,
        source_meta={"kind": "toric", "L": L, "distance": L},
    )


def bb72_problem(*, sector: str = "X", prior: str | Fraction = "1/1000") -> DecodingProblem:
    """Bravyi et al. (Nature 627, 778, 2024, arXiv:2308.07915) `[[72, 12, 6]]` bivariate bicycle code."""
    _ensure_tools_on_path()
    from ftqc_codes import codes as fc_codes

    code = fc_codes.bb(6, 6, ["x3", "y1", "y2"], ["y3", "x1", "x2"], name="BB72")
    return DecodingProblem.from_css_code(code, sector=sector, prior=prior)


def bb144_gross_problem(*, sector: str = "X", prior: str | Fraction = "1/1000") -> DecodingProblem:
    """Bravyi et al. (Nature 627, 778, 2024, arXiv:2308.07915) `[[144, 12, 12]]` gross bivariate bicycle code."""
    _ensure_tools_on_path()
    from ftqc_codes import codes as fc_codes

    code = fc_codes.bb(12, 6, ["x3", "y1", "y2"], ["y3", "x1", "x2"], name="BB144")
    return DecodingProblem.from_css_code(code, sector=sector, prior=prior)
