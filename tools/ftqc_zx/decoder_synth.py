"""Flag-Conditioned Syndrome Lookup & Decoder Compiler (Appendix C.2).

Automates Section C.2 (pages 57-59) of Rodatz, Poor, and Kissinger (arXiv:2506.17181v4):
  1) Computes the punctured syndrome parity-check matrix `H_prime` mapping the surviving
     ancilla measurements `w` onto the code's exact syndrome bits when `flag == 0`.
  2) Propagates all single-location faults in `circuit` to extract:
       - which faults trigger each flag measurement,
       - what residual physical data error `e_data` each flagged fault leaves on the
         output code block, and
       - its canonical stabilizer-coset representative `(min_weight_rep, logical_flip)`.
  3) Synthesizes a certified two-stage lookup table `(flag_mask, next_round_syndrome) -> recovery`
     and verifies zero ambiguous collisions across all weight-1 faults (Proposition 2.9)!

Zero float literals appear anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .diagram import CliffordCircuit, RefusalError, SINGLE_MEAS
from .equivalence import CSSCodeSpec
from .gf2 import solve_gf2


def _expand_rowspace(rows: tuple[int, ...], n_bits: int) -> list[int]:
    space = [0]
    for r in rows:
        space = space + [x ^ r for x in space]
    # Deduplicate preserving order
    seen: set[int] = set()
    out: list[int] = []
    for x in space:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _min_coset_rep(err_mask: int, stab_space: list[int], n_bits: int) -> tuple[int, int]:
    """Return (min_hamming_weight, canonical_min_rep_mask) in `err_mask ^ stab_space`."""
    best_wt = n_bits + 1
    best_rep = err_mask
    for s in stab_space:
        cand = err_mask ^ s
        wt = cand.bit_count()
        if wt < best_wt or (wt == best_wt and cand < best_rep):
            best_wt = wt
            best_rep = cand
    return best_wt, best_rep


def _syndrome_bits(err_mask: int, check_rows: tuple[int, ...]) -> int:
    syn = 0
    for i, r in enumerate(check_rows):
        if (err_mask & r).bit_count() & 1:
            syn |= 1 << i
    return syn


def _mask_to_qubit_list_1indexed(mask: int, n: int) -> list[int]:
    return [q + 1 for q in range(n) if (mask >> q) & 1]


def _mask_to_pauli_label(mask: int, pauli_char: str, n: int) -> str:
    if mask == 0:
        return "I"
    return "".join(f"{pauli_char}_{q + 1}" for q in range(n) if (mask >> q) & 1)


@dataclass(frozen=True)
class FlagConditionedDecoder:
    """Certified flag-conditioned syndrome decoder table (Appendix C.2)."""

    circuit_name: str
    code_name: str
    n_data: int
    extracted_sector: str                 # 'X' (extracts H_X, flips Z on data; ancilla X faults propagate to data X)
    propagated_error_type: str            # 'X' if extracted_sector == 'X', else 'Z'
    flag_meas_labels: tuple[str, ...]
    syndrome_meas_labels: tuple[str, ...]
    h_prime_rows: tuple[list[int], ...]   # r_checks x n_syndrome_meas binary matrix mapping `w` -> syndrome!
    unflagged_max_coset_weight: int       # must be <= 1 for t=1 fault tolerance!
    flagged_footprints_1indexed: tuple[list[int], ...]
    lookup_unflagged: dict[int, int]      # next_syndrome -> recovery_mask on n data qubits
    lookup_flagged: dict[tuple[int, int], int]  # (flag_mask, next_syndrome) -> recovery_mask
    unambiguous_weight1: bool
    collisions: tuple[dict[str, Any], ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "circuit_name": self.circuit_name,
            "code_name": self.code_name,
            "n_data": self.n_data,
            "extracted_sector": self.extracted_sector,
            "propagated_error_type": self.propagated_error_type,
            "flag_meas_labels": list(self.flag_meas_labels),
            "syndrome_meas_labels": list(self.syndrome_meas_labels),
            "h_prime_rows": [list(r) for r in self.h_prime_rows],
            "unflagged_max_coset_weight": self.unflagged_max_coset_weight,
            "flagged_footprints_1indexed": [list(fp) for fp in self.flagged_footprints_1indexed],
            "unambiguous_weight1": self.unambiguous_weight1,
            "lookup_unflagged": {
                str(syn): _mask_to_pauli_label(rec, self.propagated_error_type, self.n_data)
                for syn, rec in sorted(self.lookup_unflagged.items())
            },
            "lookup_flagged": {
                f"flag={fm},syn={syn}": _mask_to_pauli_label(
                    rec, self.propagated_error_type, self.n_data
                )
                for (fm, syn), rec in sorted(self.lookup_flagged.items())
            },
            "collisions": list(self.collisions),
        }


def synthesize_flag_decoder(
    circuit: CliffordCircuit,
    code: CSSCodeSpec,
    extracted_sector: str = "X",
    flag_labels: tuple[str, ...] | None = None,
) -> FlagConditionedDecoder:
    """Synthesize and certify the Appendix C.2 flag-conditioned syndrome decoder for `circuit`.

    1. Propagates incoming single-qubit data errors through the noiseless circuit to solve
       for `H_prime` such that `H_prime * w == H_code * e_in` across all single-qubit data errors!
    2. Propagates every single-location physical fault inside `circuit` forward to `(meas_flips, data_out)`,
       separating unflagged faults (`flag_mask == 0`) from flagged faults (`flag_mask != 0`),
       and builds + certifies the exact `(flag_mask, next_syndrome) -> recovery` lookup table!
    """
    circuit.validate()
    if tuple(circuit.input_qubits) != tuple(range(code.n)) or tuple(circuit.output_qubits) != tuple(range(code.n)):
        raise RefusalError(
            "DECODER-DATA-QUBIT-MISMATCH",
            f"circuit inputs/outputs must be 0..{code.n - 1}",
        )

    meas_ops: list[tuple[int, str, str]] = []
    for idx, op in enumerate(circuit.ops):
        if op.kind in SINGLE_MEAS:
            lbl = op.label or f"{op.kind}_q{op.qubits[0]}_{idx}"
            meas_ops.append((idx, op.kind, lbl))

    if flag_labels is None:
        inferred_flags = tuple(lbl for _, _, lbl in meas_ops if "flag" in lbl.lower())
    else:
        inferred_flags = flag_labels
    flag_set = set(inferred_flags)
    flag_list = [lbl for _, _, lbl in meas_ops if lbl in flag_set]
    syn_list = [lbl for _, _, lbl in meas_ops if lbl not in flag_set]

    # Forward concrete Pauli propagator `(x_mask, z_mask)` over `circuit.n_qubits`:
    def run_pauli_through_ops(
        init_x: int,
        init_z: int,
        inject_step: int | None = None,
        inject_x: int = 0,
        inject_z: int = 0,
        inject_before_meas: bool = False,
    ) -> tuple[int, int, dict[str, int]]:
        xm, zm = init_x, init_z
        outcomes: dict[str, int] = {}
        for idx, op in enumerate(circuit.ops):
            k = op.kind
            if inject_step == idx and inject_before_meas:
                xm ^= inject_x
                zm ^= inject_z
            if k in ("PREP_0", "PREP_1", "PREP_PLUS", "PREP_MINUS"):
                q = op.qubits[0]
                xm &= ~(1 << q)
                zm &= ~(1 << q)
            elif k == "H":
                q = op.qubits[0]
                xb = (xm >> q) & 1
                zb = (zm >> q) & 1
                xm = (xm & ~(1 << q)) | (zb << q)
                zm = (zm & ~(1 << q)) | (xb << q)
            elif k in ("S", "SDG"):
                q = op.qubits[0]
                if (xm >> q) & 1:
                    zm ^= 1 << q
            elif k == "CX":
                c, t = op.qubits
                if (xm >> c) & 1:
                    xm ^= 1 << t
                if (zm >> t) & 1:
                    zm ^= 1 << c
            elif k == "CZ":
                q0, q1 = op.qubits
                if (xm >> q0) & 1:
                    zm ^= 1 << q1
                if (xm >> q1) & 1:
                    zm ^= 1 << q0
            elif k == "SWAP":
                q0, q1 = op.qubits
                x0, x1 = (xm >> q0) & 1, (xm >> q1) & 1
                z0, z1 = (zm >> q0) & 1, (zm >> q1) & 1
                xm = (xm & ~((1 << q0) | (1 << q1))) | (x1 << q0) | (x0 << q1)
                zm = (zm & ~((1 << q0) | (1 << q1))) | (z1 << q0) | (z0 << q1)
            elif k in SINGLE_MEAS:
                q = op.qubits[0]
                lbl = op.label or f"{k}_q{q}_{idx}"
                if k == "MEAS_Z":
                    outcomes[lbl] = (xm >> q) & 1
                elif k == "MEAS_X":
                    outcomes[lbl] = (zm >> q) & 1
                else:
                    outcomes[lbl] = ((xm >> q) ^ (zm >> q)) & 1
                xm &= ~(1 << q)
                zm &= ~(1 << q)
            if inject_step == idx and not inject_before_meas:
                xm ^= inject_x
                zm ^= inject_z
        return xm, zm, outcomes

    # 1. Solve for punctured check matrix `H_prime` (r_checks x len(syn_list))!
    # When `extracted_sector == 'X'`, incoming Z errors on data qubits 0..n-1 produce `w` on `syn_list`,
    # and we want `H_prime * w == H_X * e_in`!
    check_rows_extract = code.hx_rows if extracted_sector == "X" else code.hz_rows
    n_syn_m = len(syn_list)
    w_samples: list[int] = []
    target_syns: list[int] = []
    for dq in range(code.n):
        ix = (1 << dq) if extracted_sector == "Z" else 0
        iz = (1 << dq) if extracted_sector == "X" else 0
        _, _, outs = run_pauli_through_ops(ix, iz)
        w_vec = 0
        for m_i, slbl in enumerate(syn_list):
            if outs[slbl]:
                w_vec |= 1 << m_i
        w_samples.append(w_vec)
        target_syns.append(_syndrome_bits(1 << dq, check_rows_extract))

    h_prime: list[list[int]] = []
    for r_idx in range(len(check_rows_extract)):
        rhs = [(ts >> r_idx) & 1 for ts in target_syns]
        ok, hp_row_mask, _ = solve_gf2(w_samples, rhs, n_syn_m)
        if not ok:
            raise RefusalError(
                "SYNDROME-EXTRACTION-INCOMPLETE",
                f"surviving measurements {syn_list} do not span syndrome check row {r_idx}",
            )
        h_prime.append([(hp_row_mask >> j) & 1 for j in range(n_syn_m)])

    # 2. Propagate all single-location faults in the bad propagation sector (`X` when `extracted_sector == 'X'`,
    #    `Z` when `extracted_sector == 'Z'`, and under full Pauli faults check all single-op faults)!
    prop_type = "X" if extracted_sector == "X" else "Z"
    stab_same = _expand_rowspace(code.hx_rows if prop_type == "X" else code.hz_rows, code.n)
    next_check_rows = code.hz_rows if prop_type == "X" else code.hx_rows
    data_mask_all = (1 << code.n) - 1

    # Standard unflagged lookup table: all weight-0 and weight-1 errors on `0..code.n - 1`
    lookup_unflagged: dict[int, int] = {0: 0}
    for dq in range(code.n):
        err = 1 << dq
        syn = _syndrome_bits(err, next_check_rows)
        lookup_unflagged[syn] = err

    unflagged_max_wt = 0
    flagged_footprints_seen: set[int] = set()
    # Track `(flag_mask, next_syndrome) -> set of stabilizer-equivalent canonical reps`
    flagged_candidates: dict[tuple[int, int], dict[int, int]] = {}

    for idx, op in enumerate(circuit.ops):
        # Enumerate CSS fault injections on `op` in `prop_type` sector (before meas for MEAS, after op otherwise):
        before_m = op.kind in SINGLE_MEAS
        q_list = op.qubits
        # Under CSS noise on `op`, any non-empty subset of `op.qubits` in `prop_type` sector is weight 1!
        for sub in range(1, 1 << len(q_list)):
            fmask = 0
            for j_i, q in enumerate(q_list):
                if (sub >> j_i) & 1:
                    fmask |= 1 << q
            inj_x = fmask if prop_type == "X" else 0
            inj_z = fmask if prop_type == "Z" else 0
            out_x, out_z, outs = run_pauli_through_ops(
                0, 0, inject_step=idx, inject_x=inj_x, inject_z=inj_z, inject_before_meas=before_m
            )
            data_err = (out_x if prop_type == "X" else out_z) & data_mask_all
            flag_mask = 0
            for f_i, flbl in enumerate(flag_list):
                if outs[flbl]:
                    flag_mask |= 1 << f_i
            c_wt, c_rep = _min_coset_rep(data_err, stab_same, code.n)
            if flag_mask == 0:
                if c_wt > unflagged_max_wt:
                    unflagged_max_wt = c_wt
            else:
                flagged_footprints_seen.add(data_err)
                nsyn = _syndrome_bits(data_err, next_check_rows)
                bucket = flagged_candidates.setdefault((flag_mask, nsyn), {})
                # Keep raw `data_err` mapped by canonical coset `c_rep` so any two stabilizer-equivalent
                # errors (`e1 ^ e2 in Stab`) merge into ONE coset entry!
                bucket.setdefault(c_rep, data_err)

    collisions: list[dict[str, Any]] = []
    lookup_flagged: dict[tuple[int, int], int] = {}

    # Seed `lookup_flagged` with standard weight-1 errors (in case flag flipped due to a pure measurement error!)
    for f_m in range(1, 1 << len(flag_list)):
        for syn, std_rec in lookup_unflagged.items():
            if (f_m, syn) not in flagged_candidates:
                lookup_flagged[(f_m, syn)] = std_rec

    for (f_m, nsyn), coset_map in sorted(flagged_candidates.items()):
        if len(coset_map) > 1:
            collisions.append(
                {
                    "flag_mask": f_m,
                    "next_syndrome": nsyn,
                    "distinct_coset_reps": [_mask_to_qubit_list_1indexed(cr, code.n) for cr in sorted(coset_map)],
                }
            )
        # Pick minimum-weight representative for recovery
        chosen_crep = sorted(coset_map.keys(), key=lambda x: (x.bit_count(), x))[0]
        lookup_flagged[(f_m, nsyn)] = chosen_crep

    unambig = (unflagged_max_wt <= (code.d - 1) // 2) and (len(collisions) == 0)
    sorted_fps = tuple(
        sorted(
            (_mask_to_qubit_list_1indexed(fp, code.n) for fp in flagged_footprints_seen),
            key=lambda lst: (len(lst), lst),
        )
    )

    return FlagConditionedDecoder(
        circuit_name=circuit.name,
        code_name=code.name,
        n_data=code.n,
        extracted_sector=extracted_sector,
        propagated_error_type=prop_type,
        flag_meas_labels=tuple(flag_list),
        syndrome_meas_labels=tuple(syn_list),
        h_prime_rows=tuple(h_prime),
        unflagged_max_coset_weight=unflagged_max_wt,
        flagged_footprints_1indexed=sorted_fps,
        lookup_unflagged=lookup_unflagged,
        lookup_flagged=lookup_flagged,
        unambiguous_weight1=unambig,
        collisions=tuple(collisions),
    )
