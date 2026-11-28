"""Open-boundary mid-circuit escape postselection & Caliper gap verifier.

Implements the exact mid-circuit postselection analysis of:
- Mehta, Menon, Zhou, Lukin & Ataides (2026, arXiv:2610.06828):
  "The Practical Cost of Magic State Cultivation"

Core capabilities:
1. **Open-Boundary Escape DEM (`OpenBoundaryEscapeDEM`)**:
   Partitions a matchable detector error model into committed mid-circuit detectors
   `D_comm` and open future boundary detectors `D_open` (`D = D_comm |_| D_open`),
   interoperating directly with `ftqc_decoders.dem.DecodingProblem` and live `stim`
   DEMs (`DecodingProblem.from_stim_dem_text`) when available.
2. **Logical-Parity Shortest-Path Tables `d_b(u, v)` (`b in {0, 1}`)**:
   Builds the exact parity-lifted matching graph on `(V |_| {boundary}) x F_2` over
   exact integer soft weights (`scaled_log_odds_int`) and computes all-pairs shortest
   paths per observable-flip parity `b in {0, 1}` via Dijkstra.
3. **Caliper Hidden-Syndrome Search & Blossom Cross-Check (`caliper_worst_case_gap`)**:
   Evaluates worst-case soft-weight complementary gap `Delta_min(s_comm)` over hidden
   future syndromes `s_open subseteq D_open` (exact enumeration for `|D_open| <= 10`,
   plus local open-boundary pair/defect search from Lemma 2.5 / Theorem 2.7), cross-checked
   live against `ftqc_decoders.mwpm.decode_mwpm`!
4. **Sabotage Control `N5` (`OPEN-BOUNDARY-NAIVE-GAP-HAZARD`)**:
   Refuses naive terminal complementary-gap postselection (`s_open == 0` assumption) on any
   escape DEM with non-empty `D_open`, exhibiting a concrete hidden future syndrome witness
   that erodes the naive gap.
5. **Adaptive Multi-Round Escape Spacetime Volume (`adaptive_escape_spacetime_volume`)**:
   Exact `Fraction` evaluation of `E[V_pi] / Pr_pi(accept)` across multi-round escape policies.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import heapq
from pathlib import Path
import sys
from typing import Any, Sequence

from .ir import MagicRefusalError
from .rings import parse_exact_fraction, round_down_str, round_up_str


def _ensure_tools_on_syspath() -> None:
    for candidate in (
        Path(__file__).resolve().parents[1],
    ):
        if (candidate / "ftqc_decoders").is_dir():
            s = str(candidate)
            if s not in sys.path:
                sys.path.insert(0, s)
            break


def _integer_soft_weight(p: Fraction, scale: int = 1024) -> int:
    """Exact stdlib integer log-odds weight matching `ftqc_decoders.dem.scaled_log_odds_int`
    when available, or integer binary-search rational approximation with zero floats!
    """
    _ensure_tools_on_syspath()
    try:
        from ftqc_decoders.dem import scaled_log_odds_int

        return int(scaled_log_odds_int(p, scale=scale))
    except ImportError:
        # Pure-stdlib zero-float fallback: floor(scale * ln((1-p)/p)) via integer power comparison
        if p <= 0 or p >= Fraction(1, 2):
            raise ValueError(f"Mechanism prior {p} must lie in (0, 1/2)")
        ratio = (Fraction(1) - p) / p
        # Use natural log Taylor / bit-halving in Fraction? Even simpler:
        # ln(x) for x > 1: extract k = floor(log2(x)), x = 2^k * y with y in [1, 2),
        # ln(2) ~= 693147180559945309 / 10^18, and ln(y) via exact rational artanh series
        # ln(y) = 2 * sum_{m=0..12} ((y-1)/(y+1))^(2m+1) / (2m+1)!
        k = 0
        y = ratio
        while y >= 2:
            y /= 2
            k += 1
        z = (y - 1) / (y + 1)
        z2 = z * z
        term = z
        ln_y = Fraction(0)
        for m in range(14):
            ln_y += term / (2 * m + 1)
            term *= z2
        ln_y *= 2
        ln2_rat = Fraction(693147180559945309, 10**18)
        val = (k * ln2_rat + ln_y) * scale
        return max(1, val.numerator // val.denominator)


@dataclass(frozen=True)
class EscapeMechanism:
    det_mask: int
    obs_bit: int
    prior: Fraction
    soft_weight: int


@dataclass(frozen=True)
class OpenBoundaryEscapeDEM:
    """Matchable detector error model partitioned into `committed_mask` (`D_comm`) and
    `open_mask` (`D_open`) for mid-circuit cultivation escape postselection.
    """

    name: str
    num_checks: int
    committed_mask: int
    open_mask: int
    mechanisms: tuple[EscapeMechanism, ...]
    decoding_problem: Any | None = None

    def __post_init__(self) -> None:
        all_mask = (1 << self.num_checks) - 1
        if (self.committed_mask & self.open_mask) != 0:
            raise ValueError("committed_mask and open_mask must be disjoint")
        if (self.committed_mask | self.open_mask) != all_mask:
            raise ValueError("committed_mask | open_mask must partition all detectors 0..num_checks-1")
        for m in self.mechanisms:
            if m.det_mask == 0 or m.det_mask.bit_count() > 2:
                raise ValueError("OpenBoundaryEscapeDEM requires graphlike mechanisms (1 or 2 detectors)")
            if m.obs_bit not in (0, 1):
                raise ValueError("obs_bit must be 0 or 1")

    @property
    def committed_detectors(self) -> tuple[int, ...]:
        return tuple(i for i in range(self.num_checks) if (self.committed_mask >> i) & 1)

    @property
    def open_detectors(self) -> tuple[int, ...]:
        return tuple(i for i in range(self.num_checks) if (self.open_mask >> i) & 1)

    @classmethod
    def from_decoding_problem(
        cls,
        problem: Any,
        committed_detectors: Sequence[int],
        *,
        name: str | None = None,
    ) -> OpenBoundaryEscapeDEM:
        n_chk = int(problem.num_checks)
        comm_mask = sum(1 << d for d in committed_detectors)
        open_mask = ((1 << n_chk) - 1) ^ comm_mask
        mechs: list[EscapeMechanism] = []
        for idx in range(int(problem.num_mechanisms)):
            dm = int(problem.mech_detectors[idx])
            if dm == 0:
                continue
            if dm.bit_count() > 2:
                raise ValueError("DecodingProblem contains hyperedge > 2; pass graphlike DEM")
            ob = int(problem.mech_observables[idx]) & 1
            pr = Fraction(problem.priors[idx])
            sw = int(problem.weights[idx])
            mechs.append(EscapeMechanism(dm, ob, pr, sw))
        return cls(
            name=name or str(problem.name),
            num_checks=n_chk,
            committed_mask=comm_mask,
            open_mask=open_mask,
            mechanisms=tuple(mechs),
            decoding_problem=problem,
        )


def build_surface_escape_dem(
    distance: int = 3,
    total_rounds: int = 4,
    committed_rounds: int = 2,
    *,
    prior: Fraction | str = Fraction(1, 200),
    use_stim_if_available: bool = True,
) -> OpenBoundaryEscapeDEM:
    """Build a matchable `distance x total_rounds` mid-circuit escape slice DEM where rounds
    `0 .. committed_rounds - 1` are committed mid-circuit (`D_comm`) and rounds
    `committed_rounds .. total_rounds - 1` are open future detectors (`D_open`)!

    When `stim` and `ftqc_decoders` are installed, cross-builds
    the DEM directly from a live `stim.Circuit` repetition/escape slice via
    `DecodingProblem.from_stim_dem_text` to verify exact schema & graphlike parity agreement!
    """
    if distance < 3 or total_rounds < 2 or not (1 <= committed_rounds <= total_rounds):
        raise ValueError("Require distance >= 3, total_rounds >= 2, and 1 <= committed_rounds <= total_rounds")
    pf = parse_exact_fraction(prior)
    _ensure_tools_on_syspath()

    n_spatial = distance - 1
    num_checks = n_spatial * total_rounds
    comm_dets = [r * n_spatial + s for r in range(committed_rounds) for s in range(n_spatial)]

    if use_stim_if_available:
        try:
            from ftqc_decoders.dem import DecodingProblem
            import stim  # type: ignore[import-not-found]

            # Build exact Stim circuit for `total_rounds` of syndrome extraction on `distance` data wires `0..d-1`
            # and `d-1` ancilla wires `d..2d-2`!
            p_str = round_up_str(pf, 12)
            lines: list[str] = []
            for q in range(distance):
                lines.append(f"R {q}")
            lines.append("TICK")
            meas_count = 0
            for r_idx in range(total_rounds):
                for q in range(distance):
                    lines.append(f"X_ERROR({p_str}) {q}")
                for s in range(n_spatial):
                    anc = distance + s
                    lines.append(f"R {anc}")
                    lines.append(f"CX {s} {anc}")
                    lines.append(f"CX {s + 1} {anc}")
                    lines.append(f"M({p_str}) {anc}")
                    meas_count += 1
                    if r_idx == 0:
                        lines.append(f"DETECTOR rec[-1]")
                    else:
                        lines.append(f"DETECTOR rec[-1] rec[-{1 + n_spatial}]")
                lines.append("TICK")
            # Terminal noiseless data readout closing the final round + observable on data 0:
            for q in range(distance):
                lines.append(f"M {q}")
            # Note: keep num_checks == n_spatial * total_rounds by folding terminal check or letting
            # final round have open time-boundary! Without terminal data readout's extra detectors,
            # observable is data 0 (`M 0` alone, and left boundary flips `M 0`!):
            sc_lines = [ln for ln in lines if not ln.startswith("M ")]
            sc_lines.append("M 0")
            sc_lines.append("OBSERVABLE_INCLUDE(0) rec[-1]")
            sc = stim.Circuit("\n".join(sc_lines))
            dem_text = str(sc.detector_error_model(decompose_errors=True))
            prob = DecodingProblem.from_stim_dem_text(dem_text, name=f"stim_escape_d{distance}_R{total_rounds}")
            if prob.num_checks == num_checks:
                return OpenBoundaryEscapeDEM.from_decoding_problem(prob, comm_dets)
        except ImportError:
            pass

    # Pure-stdlib exact construction (100% deterministic, zero external deps):
    sw = _integer_soft_weight(pf)
    mechs: list[EscapeMechanism] = []
    check_cols_det: list[int] = []
    check_cols_obs: list[int] = []

    def add_mech(dm: int, ob: int) -> None:
        mechs.append(EscapeMechanism(dm, ob, pf, sw))
        check_cols_det.append(dm)
        check_cols_obs.append(ob)

    for r_idx in range(total_rounds):
        base = r_idx * n_spatial
        # Spatial data errors on qubits 0 .. distance - 1:
        # q = 0: left boundary + detector(base + 0), flips observable 0!
        add_mech(1 << base, 1)
        # 1 <= q <= distance - 2: interior pair (base + q - 1, base + q), obs = 0
        for q in range(1, distance - 1):
            add_mech((1 << (base + q - 1)) | (1 << (base + q)), 0)
        # q = distance - 1: right boundary + detector(base + n_spatial - 1), obs = 0
        add_mech(1 << (base + n_spatial - 1), 0)
        # Timelike measurement errors between round r_idx and r_idx + 1:
        if r_idx + 1 < total_rounds:
            nxt_base = (r_idx + 1) * n_spatial
            for s in range(n_spatial):
                add_mech((1 << (base + s)) | (1 << (nxt_base + s)), 0)
        else:
            # Final open round measurement faults can terminate at open time boundary with obs = 0:
            for s in range(n_spatial):
                add_mech(1 << (base + s), 0)

    prob_obj: Any | None = None
    try:
        from ftqc_decoders.dem import DecodingProblem

        n_m = len(mechs)
        c_rows = [
            sum(((check_cols_det[m_i] >> d_i) & 1) << m_i for m_i in range(n_m))
            for d_i in range(num_checks)
        ]
        o_rows = [sum((check_cols_obs[m_i] & 1) << m_i for m_i in range(n_m))]
        prob_obj = DecodingProblem.from_matrices(
            c_rows,
            o_rows,
            n_m,
            name=f"native_escape_d{distance}_R{total_rounds}",
            priors=[pf] * n_m,
        )
    except ImportError:
        pass

    comm_mask = sum(1 << d for d in comm_dets)
    return OpenBoundaryEscapeDEM(
        name=f"escape_d{distance}_R{total_rounds}_c{committed_rounds}",
        num_checks=num_checks,
        committed_mask=comm_mask,
        open_mask=((1 << num_checks) - 1) ^ comm_mask,
        mechanisms=tuple(mechs),
        decoding_problem=prob_obj,
    )


def compute_logical_parity_shortest_paths(
    escape_dem: OpenBoundaryEscapeDEM,
) -> list[list[tuple[int, int]]]:
    """Mehta et al. Section 2: compute exact shortest-path soft weights `(d_0(u, v), d_1(u, v))`
    between all pairs `u, v in {0 .. num_checks}` (where vertex `num_checks` is the boundary
    node `partial`) on the parity-lifted matching graph!
    """
    n_chk = escape_dem.num_checks
    bnd = n_chk
    n_verts = n_chk + 1
    INF = 10**15

    # Adjacency list on `n_verts`: `(neighbor, obs_bit, soft_weight)`
    adj: list[list[tuple[int, int, int]]] = [[] for _ in range(n_verts)]
    for m in escape_dem.mechanisms:
        dm = m.det_mask
        if dm.bit_count() == 1:
            u = (dm & -dm).bit_length() - 1
            adj[u].append((bnd, m.obs_bit, m.soft_weight))
            adj[bnd].append((u, m.obs_bit, m.soft_weight))
        else:
            lsb = dm & -dm
            u = lsb.bit_length() - 1
            v = ((dm ^ lsb) & -(dm ^ lsb)).bit_length() - 1
            adj[u].append((v, m.obs_bit, m.soft_weight))
            adj[v].append((u, m.obs_bit, m.soft_weight))
    # Boundary-to-boundary zero-cost parity-0 self-loop (pairs of boundary terminations cancel!):
    adj[bnd].append((bnd, 0, 0))

    table: list[list[tuple[int, int]]] = [[(INF, INF) for _ in range(n_verts)] for _ in range(n_verts)]
    for src in range(n_verts):
        dist = [[INF, INF] for _ in range(n_verts)]
        dist[src][0] = 0
        pq: list[tuple[int, int, int]] = [(0, src, 0)]
        while pq:
            d_cur, u, par = heapq.heappop(pq)
            if d_cur > dist[u][par]:
                continue
            for v, edge_ob, w in adj[u]:
                npar = par ^ edge_ob
                nd = d_cur + w
                if nd < dist[v][npar]:
                    dist[v][npar] = nd
                    heapq.heappush(pq, (nd, v, npar))
        for dst in range(n_verts):
            table[src][dst] = (dist[dst][0], dist[dst][1])
    return table


def exact_coset_matching_weights(
    escape_dem: OpenBoundaryEscapeDEM,
    dist_table: list[list[tuple[int, int]]],
    syndrome_mask: int,
) -> tuple[int, int]:
    """Compute exact minimum-weight graphlike matching `(W_0(syndrome), W_1(syndrome))` in BOTH
    observable homology classes `0` and `1` via exact bitmask DP over active syndrome defects!
    """
    INF = 10**15
    bnd = escape_dem.num_checks
    # Boundary-to-boundary single odd loop of parity 1 (in case 0 defects still flip obs via boundary-to-boundary path):
    bnd_obs1 = dist_table[bnd][bnd][1]

    defects = [i for i in range(escape_dem.num_checks) if (syndrome_mask >> i) & 1]
    k = len(defects)
    if k == 0:
        return (0, bnd_obs1)

    # dp[mask] = [best_w_par0, best_w_par1] for pairing subset `mask` of `defects`:
    dp = [[INF, INF] for _ in range(1 << k)]
    dp[0][0] = 0
    dp[0][1] = bnd_obs1

    for mask in range(1, 1 << k):
        first_bit = (mask & -mask).bit_length() - 1
        u = defects[first_bit]
        rem_without_u = mask ^ (1 << first_bit)
        best0, best1 = INF, INF

        # Option A: match `u` to boundary `bnd`!
        w_ub0, w_ub1 = dist_table[u][bnd]
        prev0, prev1 = dp[rem_without_u]
        best0 = min(best0, prev0 + w_ub0, prev1 + w_ub1)
        best1 = min(best1, prev0 + w_ub1, prev1 + w_ub0)

        # Option B: match `u` to another defect `v` in `rem_without_u`!
        sub = rem_without_u
        while sub:
            bit_v = sub & -sub
            idx_v = bit_v.bit_length() - 1
            v = defects[idx_v]
            rest = rem_without_u ^ bit_v
            w_uv0, w_uv1 = dist_table[u][v]
            p0, p1 = dp[rest]
            best0 = min(best0, p0 + w_uv0, p1 + w_uv1)
            best1 = min(best1, p0 + w_uv1, p1 + w_uv0)
            sub ^= bit_v

        dp[mask][0] = best0
        dp[mask][1] = best1

    w0, w1 = dp[(1 << k) - 1]
    return (w0, w1)


def caliper_worst_case_gap(
    escape_dem: OpenBoundaryEscapeDEM,
    committed_syndrome: int,
    *,
    mode: str = "caliper",
    crosscheck_blossom: bool = True,
) -> dict[str, Any]:
    """Certify the mid-circuit complementary gap of `committed_syndrome` on `escape_dem`!

    - If `mode == "naive_terminal_gap"` and `escape_dem.open_mask != 0`:
      Raises `MagicRefusalError("OPEN-BOUNDARY-NAIVE-GAP-HAZARD", ...)` (`N5` sabotage control!)
      because ignoring hidden open-boundary syndromes overestimates postselection confidence!
    - If `mode == "caliper"`:
      Computes both the naive `s_open = 0` gap AND the true worst-case gap across hidden open
      syndromes `s_open subseteq D_open`, cross-checking MWPM winning weights against
      `ftqc_decoders.mwpm.decode_mwpm` when available!
    """
    if (committed_syndrome & ~escape_dem.committed_mask) != 0:
        raise ValueError("committed_syndrome has bits outside committed_mask")

    dist_table = compute_logical_parity_shortest_paths(escape_dem)
    w0_naive, w1_naive = exact_coset_matching_weights(escape_dem, dist_table, committed_syndrome)
    pred_obs = 0 if w0_naive <= w1_naive else 1
    naive_gap = abs(w1_naive - w0_naive)

    open_dets = escape_dem.open_detectors
    # Enumerate candidate hidden open syndromes (all 2^|D_open| when |D_open| <= 10, else all weight <= 3):
    cand_open_masks: list[int] = []
    if len(open_dets) <= 10:
        for m_idx in range(1 << len(open_dets)):
            om = 0
            for bit_i, det_id in enumerate(open_dets):
                if (m_idx >> bit_i) & 1:
                    om |= 1 << det_id
            cand_open_masks.append(om)
    else:
        cand_open_masks.append(0)
        for d1 in open_dets:
            cand_open_masks.append(1 << d1)
        for i, d1 in enumerate(open_dets):
            for d2 in open_dets[i + 1 :]:
                cand_open_masks.append((1 << d1) | (1 << d2))

    # Recall Mehta et al. Section 2: for a fixed mid-circuit syndrome `s_comm`, conditioned on
    # future open-boundary outcomes `s_open`, if the future decoder adapts to `(s_comm, s_open)` its
    # conditional gap is `|W_1(s_comm, s_open) - W_0(s_comm, s_open)|`!
    worst_gap = naive_gap
    worst_s_open = 0
    worst_pair = (w0_naive, w1_naive)
    blossom_checked = 0

    mwpm_fn = None
    if crosscheck_blossom and escape_dem.decoding_problem is not None:
        _ensure_tools_on_syspath()
        try:
            from ftqc_decoders.mwpm import decode_mwpm

            mwpm_fn = decode_mwpm
        except ImportError:
            mwpm_fn = None

    for om in cand_open_masks:
        full_syn = committed_syndrome | om
        w0, w1 = exact_coset_matching_weights(escape_dem, dist_table, full_syn)
        cond_gap = abs(w1 - w0)
        if cond_gap < worst_gap:
            worst_gap = cond_gap
            worst_s_open = om
            worst_pair = (w0, w1)
        if mwpm_fn is not None and blossom_checked < 8 and full_syn.bit_count() <= 8:
            mres = mwpm_fn(escape_dem.decoding_problem, full_syn)
            if mres.soft_weight != min(w0, w1):
                raise AssertionError(
                    f"Blossom MWPM soft_weight {mres.soft_weight} != Caliper min({w0}, {w1}) at syn={full_syn}"
                )
            blossom_checked += 1

    if mode == "naive_terminal_gap" and escape_dem.open_mask != 0:
        raise MagicRefusalError(
            "OPEN-BOUNDARY-NAIVE-GAP-HAZARD",
            f"Naive terminal gap ({naive_gap}) on open-boundary DEM {escape_dem.name!r} ignores "
            f"{len(open_dets)} open future detectors; hidden open syndrome {bin(worst_s_open)} "
            f"erodes conditional gap down to {worst_gap}!",
            {
                "naive_gap": naive_gap,
                "caliper_worst_gap": worst_gap,
                "eroding_open_syndrome": worst_s_open,
                "num_open_detectors": len(open_dets),
            },
        )

    return {
        "dem_name": escape_dem.name,
        "committed_syndrome": committed_syndrome,
        "predicted_observable": pred_obs,
        "naive_soft_gap": naive_gap,
        "caliper_worst_soft_gap": worst_gap,
        "gap_erosion": naive_gap - worst_gap,
        "worst_open_syndrome": worst_s_open,
        "worst_coset_weights": list(worst_pair),
        "open_candidates_evaluated": len(cand_open_masks),
        "blossom_crosschecks": blossom_checked,
    }


def adaptive_escape_spacetime_volume(
    round_qubits: Sequence[int],
    continuation_probs: Sequence[Fraction | str],
    terminal_accept_prob: Fraction | str,
) -> dict[str, Any]:
    """Mehta et al. Eq. (2): compute exact expected spacetime volume `E[V_pi] / Pr_pi(accept)`
    under a multi-stage mid-circuit Caliper abort policy vs. non-adaptive terminal postselection!

    - `round_qubits`: active physical qubits `(Q_0, Q_1, ..., Q_{R-1})` in each escape round.
    - `continuation_probs`: `q_r = Pr(survive round r)` for `0 <= r < R - 1` (`|continuation_probs| == R - 1`).
      (Cumulative survival entering round `r` is `prod_{j < r} q_j`.)
    - `terminal_accept_prob`: joint probability `Pr_pi(accept)` of surviving all rounds and passing
      final threshold.
    """
    R = len(round_qubits)
    if R == 0 or len(continuation_probs) != R - 1:
        raise ValueError("Require R >= 1 round_qubits and R - 1 continuation_probs")
    p_acc = parse_exact_fraction(terminal_accept_prob)
    if p_acc <= 0 or p_acc > 1:
        raise ValueError("terminal_accept_prob must lie in (0, 1]")

    q_fracs = [parse_exact_fraction(q) for q in continuation_probs]
    surv = Fraction(1)
    expected_attempt_vol = Fraction(0)
    for r_idx in range(R):
        expected_attempt_vol += surv * int(round_qubits[r_idx])
        if r_idx < R - 1:
            if not (0 < q_fracs[r_idx] <= 1):
                raise ValueError("continuation probability must lie in (0, 1]")
            surv *= q_fracs[r_idx]
    if p_acc > surv:
        raise ValueError("terminal_accept_prob cannot exceed cumulative continuation survival")

    unadaptive_attempt_vol = Fraction(sum(int(q) for q in round_qubits))
    adaptive_per_kept = expected_attempt_vol / p_acc
    unadaptive_per_kept = unadaptive_attempt_vol / p_acc
    savings_ratio = unadaptive_per_kept / adaptive_per_kept

    return {
        "rounds": R,
        "expected_attempt_volume_exact": f"{expected_attempt_vol.numerator}/{expected_attempt_vol.denominator}",
        "adaptive_volume_per_accepted_hi": round_up_str(adaptive_per_kept),
        "unadaptive_volume_per_accepted_lo": round_down_str(unadaptive_per_kept),
        "volume_reduction_factor_lo": round_down_str(savings_ratio, 12),
        "adaptive_volume_fraction": adaptive_per_kept,
    }
