"""Multivariate renormalization contraction and universal threshold certificates.

Certifies exact rational contraction of extended-rectangle failure rates across
concatenation levels k >= 1 under either adversarial local-stochastic physical noise
or independent circuit-level depolarizing hardware noise, coupled to either
recursive noisy 15-to-1 Reed-Muller magic-state distillation or certified staged
non-Clifford gadgets.

All arithmetic uses exact integers and fractions.Fraction; multi-level iterates use
monotone upward rational rounding (tails.up) so every reported rate, threshold lower
bound, and total-variation upper bound is conservative.
"""

from fractions import Fraction
from math import comb, isqrt

from .filters import decimal
from .stochastic import probability
from .supports import fields
from .tails import binomial_tail, sci, up

SIG = 60


def _pos_rational(value, label):
    if type(value) not in (str, int):
        raise ValueError(f"{label} must be an exact integer or fraction/decimal string")
    try:
        q = Fraction(value)
    except (ValueError, ZeroDivisionError) as err:
        raise ValueError(f"invalid rational for {label}") from err
    if q <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return q


def _nonneg_rational(value, label):
    if type(value) not in (str, int):
        raise ValueError(f"{label} must be an exact integer or fraction/decimal string")
    try:
        q = Fraction(value)
    except (ValueError, ZeroDivisionError) as err:
        raise ValueError(f"invalid rational for {label}") from err
    if q < 0:
        raise ValueError(f"{label} must be nonnegative")
    return q


def _pos_int(value, label):
    if type(value) is not int or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _nonneg_int(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _grouped_e3(role_locs, rates):
    """Exact degree-3 elementary symmetric polynomial over location multiplicities."""
    s1 = Fraction(0)
    s2 = Fraction(0)
    s3 = Fraction(0)
    for r, n_r in role_locs.items():
        if n_r:
            q = rates[r]
            s1 += n_r * q
            s2 += n_r * q * q
            s3 += n_r * q * q * q
    return (s1 * s1 * s1 - 3 * s1 * s2 + 2 * s3) / 6


def _parse_role_pairs(pair_map, allowed_roles, label, rational=False):
    if not isinstance(pair_map, dict):
        raise ValueError(f"{label} must be a dict keyed by 'roleA,roleB'")
    out = {}
    for key, val in pair_map.items():
        a, sep, b = key.partition(",")
        if not sep or a not in allowed_roles or b not in allowed_roles:
            raise ValueError(f"{label}: unknown role pair {key!r}")
        rk = (a, b) if a <= b else (b, a)
        v = _nonneg_rational(val, f"{label}[{key}]") if rational else Fraction(_nonneg_int(val, f"{label}[{key}]"))
        out[rk] = out.get(rk, Fraction(0)) + v
    return out


def _eval_ls_clifford(roles_cfg, rates, clamp=True):
    """Exact all-q local-stochastic ExRec upper-bound map on Clifford roles."""
    nxt = {}
    for r, cfg in roles_cfg.items():
        quad = Fraction(0)
        for (a, b), coeff in cfg["supp_pairs"].items():
            quad += coeff * rates[a] * rates[b]
        cub = _grouped_e3(cfg["role_locs"], rates)
        val = quad + cub
        nxt[r] = min(Fraction(1), val) if clamp else val
    return nxt


def _eval_dep_level1_monotone(roles_cfg, p, clamp=True):
    """Globally monotone all-p upper bound for level 0 -> 1 under depolarizing noise."""
    out = {}
    for r, cfg in roles_cfg.items():
        val = cfg["dep_total"] * p * p + binomial_tail(cfg["locations"], 2, p)
        out[r] = min(Fraction(1), val) if clamp else val
    return out


def _eval_dep_level1_pointwise(roles_cfg, p):
    """Exact pointwise level 0 -> 1 upper bound at a fixed p in [0, 1]."""
    out = {}
    one_minus_p = 1 - p
    for r, cfg in roles_cfg.items():
        L = cfg["locations"]
        factor = one_minus_p ** max(0, L - 2)
        val = cfg["dep_total"] * p * p * factor + binomial_tail(L, 2, p)
        out[r] = min(Fraction(1), val)
    return out


def _eval_magic_step(nc_cfg, prev_rates, cur_cliff_rates, level_idx, noise, p0=None, clamp=True):
    if nc_cfg is None:
        return None
    scheme = nc_cfg["scheme"]
    if scheme == "rm15-distill":
        q = prev_rates["magic"]
        slots = nc_cfg["total_slots"]
        hp = nc_cfg["harmful_pairs"]
        dist_fail = hp * q * q + comb(slots, 3) * (q ** 3)
        cliff_fail = sum((cnt * cur_cliff_rates[r] for r, cnt in nc_cfg["clifford_counts"].items()), Fraction(0))
        val = dist_fail + cliff_fail
        return min(Fraction(1), val) if clamp else val
    # scheme == "staged-gadget"
    lin = nc_cfg["linear_floor"] * (p0 if level_idx == 1 else prev_rates["magic"])
    if level_idx == 1 and noise == "depolarizing":
        val = lin + nc_cfg["dep_pairs"] * p0 * p0 + binomial_tail(nc_cfg["locations"], 2, p0)
    else:
        quad = Fraction(0)
        for (a, b), coeff in nc_cfg["supp_pairs"].items():
            quad += coeff * prev_rates[a] * prev_rates[b]
        cub = _grouped_e3(nc_cfg["role_locs"], prev_rates)
        val = lin + quad + cub
    return min(Fraction(1), val) if clamp else val


def renormalize_result(spec):
    req = ["kind", "noise", "p", "roles", "cone_weights"]
    for opt in ("non_clifford", "max_levels", "workload"):
        if opt in spec:
            req.append(opt)
    fields(spec, *req)

    noise = spec["noise"]
    if noise not in ("local-stochastic", "depolarizing"):
        raise ValueError("noise must be 'local-stochastic' or 'depolarizing'")
    p = probability(spec["p"])
    max_levels = _pos_int(spec.get("max_levels", 6), "max_levels")
    if max_levels > 12:
        raise ValueError("max_levels exceeds 12")

    raw_roles = spec["roles"]
    if not isinstance(raw_roles, dict) or not raw_roles:
        raise ValueError("roles must be a nonempty mapping")
    role_names = sorted(raw_roles)

    verified_receipts = 0
    roles_cfg = {}
    for r in role_names:
        rc = raw_roles[r]
        opt_keys = [k for k in ("qubits", "one_rec_locations", "one_rec_role_counts", "receipts") if k in rc]
        base_keys = ["locations", "role_locations", "support_pairs_by_role"]
        if noise == "depolarizing" or "depolarizing_pairs_by_role" in rc:
            base_keys.append("depolarizing_pairs_by_role")
        fields(rc, *base_keys, *opt_keys)
        L = _pos_int(rc["locations"], f"roles[{r}].locations")
        rlocs = {k: _nonneg_int(v, f"roles[{r}].role_locations[{k}]") for k, v in rc["role_locations"].items()}
        if not set(rlocs) <= set(role_names) or sum(rlocs.values()) != L:
            raise ValueError(f"roles[{r}].role_locations must use declared roles and sum to locations={L}")
        spairs = _parse_role_pairs(rc["support_pairs_by_role"], set(role_names), f"roles[{r}].support_pairs_by_role")
        supp_tot = sum(spairs.values(), Fraction(0))
        if supp_tot > comb(L, 2):
            raise ValueError(f"roles[{r}] support pairs exceed C({L}, 2)")
        dpairs = {}
        dep_tot = Fraction(0)
        if "depolarizing_pairs_by_role" in rc:
            dpairs = _parse_role_pairs(
                rc["depolarizing_pairs_by_role"], set(role_names),
                f"roles[{r}].depolarizing_pairs_by_role", rational=True,
            )
            dep_tot = sum(dpairs.values(), Fraction(0))
            if dep_tot > supp_tot:
                raise ValueError(f"roles[{r}] depolarizing pair sum exceeds support pair count")
        elif noise == "depolarizing":
            raise ValueError(f"roles[{r}] missing depolarizing_pairs_by_role")
        if "receipts" in rc:
            from ftqc_gadgets import certify as gcert
            for g_rec in rc["receipts"]:
                gcert.verify(g_rec, full=False)
                verified_receipts += 1
        one_rec_counts = None
        if "one_rec_role_counts" in rc:
            one_rec_counts = {
                k: _nonneg_int(v, f"roles[{r}].one_rec_role_counts[{k}]")
                for k, v in rc["one_rec_role_counts"].items()
            }
            if not set(one_rec_counts) <= set(role_names):
                raise ValueError(f"roles[{r}].one_rec_role_counts references undeclared role")
        roles_cfg[r] = dict(
            locations=L,
            role_locs=rlocs,
            supp_pairs=spairs,
            supp_total=int(supp_tot),
            dep_pairs=dpairs,
            dep_total=dep_tot,
            qubits=_pos_int(rc.get("qubits", 7), f"roles[{r}].qubits"),
            one_rec_locations=_pos_int(rc.get("one_rec_locations", L), f"roles[{r}].one_rec_locations"),
            one_rec_role_counts=one_rec_counts,
        )

    nc_cfg = None
    all_vars = list(role_names)
    if "non_clifford" in spec and spec["non_clifford"] is not None:
        nc = spec["non_clifford"]
        scheme = nc.get("scheme")
        all_vars.append("magic")
        if scheme == "rm15-distill":
            fields(nc, "scheme", "harmful_pairs", "total_slots", "injection_multiplier", "clifford_counts")
            hp = _pos_int(nc["harmful_pairs"], "harmful_pairs")
            slots = _pos_int(nc["total_slots"], "total_slots")
            if slots < 3 or hp > comb(slots, 2):
                raise ValueError("invalid distillation slot/pair parameters")
            cc = {k: _nonneg_int(v, f"clifford_counts[{k}]") for k, v in nc["clifford_counts"].items()}
            if not set(cc) <= set(role_names):
                raise ValueError("non_clifford.clifford_counts references undeclared Clifford role")
            nc_cfg = dict(
                scheme=scheme,
                harmful_pairs=hp,
                total_slots=slots,
                inj_mult=_pos_rational(nc["injection_multiplier"], "injection_multiplier"),
                clifford_counts=cc,
                linear_floor=Fraction(0),
            )
        elif scheme == "staged-gadget":
            s_keys = ["scheme", "locations", "role_locations", "support_pairs_by_role", "linear_floor_coeff"]
            for opt in ("depolarizing_pair_coeff", "staged_receipt", "injection_multiplier"):
                if opt in nc:
                    s_keys.append(opt)
            fields(nc, *s_keys)
            if "staged_receipt" in nc:
                from ftqc_gadgets.staged import verify as verify_staged
                verify_staged(nc["staged_receipt"])
                verified_receipts += 1
            Ls = _pos_int(nc["locations"], "non_clifford.locations")
            rlocs_s = {k: _nonneg_int(v, f"non_clifford.role_locations[{k}]") for k, v in nc["role_locations"].items()}
            if not set(rlocs_s) <= set(all_vars) or sum(rlocs_s.values()) != Ls:
                raise ValueError("non_clifford.role_locations must sum to non_clifford.locations")
            sp_s = _parse_role_pairs(nc["support_pairs_by_role"], set(all_vars), "non_clifford.support_pairs_by_role")
            dep_s = _nonneg_rational(nc.get("depolarizing_pair_coeff", str(sum(sp_s.values(), Fraction(0)))),
                                     "non_clifford.depolarizing_pair_coeff")
            nc_cfg = dict(
                scheme=scheme,
                locations=Ls,
                role_locs=rlocs_s,
                supp_pairs=sp_s,
                dep_pairs=dep_s,
                linear_floor=_nonneg_rational(nc["linear_floor_coeff"], "linear_floor_coeff"),
                inj_mult=_pos_rational(nc.get("injection_multiplier", "1"), "injection_multiplier"),
            )
        else:
            raise ValueError(f"unknown non_clifford scheme {scheme!r}")

    cw_raw = spec["cone_weights"]
    if cw_raw == "auto":
        # Match cone weights to the normalized level-1 role profile M_r^(L1) / M_max^(L1)
        # so every Clifford role enters the contraction cone simultaneously at p_th without
        # bottlenecking small roles, then set c_magic by AM-GM + level-1 entry balance.
        l1_coeffs = {
            r: (cfg["dep_total"] if noise == "depolarizing" else Fraction(cfg["supp_total"]))
            for r, cfg in roles_cfg.items()
        }
        m_max = max(max(l1_coeffs.values()), Fraction(1))
        c_vec = {r: up(max(Fraction(1, 10000), v / m_max), 6) for r, v in l1_coeffs.items()}
        if nc_cfg is not None:
            if nc_cfg["scheme"] == "rm15-distill":
                B_ff = sum(
                    (cnt * sum((coeff * c_vec[a] * c_vec[b] for (a, b), coeff in roles_cfg[r]["supp_pairs"].items()),
                               Fraction(0))
                     for r, cnt in nc_cfg["clifford_counts"].items()),
                    Fraction(0),
                )
                ratio_scaled = (B_ff / nc_cfg["harmful_pairs"]) * (10 ** 12)
                sq_num = isqrt(max(1, ratio_scaled.numerator // ratio_scaled.denominator))
                c_amgm = Fraction(sq_num, 10 ** 6)
                c_l1_magic = (nc_cfg["harmful_pairs"] * (nc_cfg["inj_mult"] ** 2)) / m_max
                c_vec["magic"] = up(max(Fraction(1, 100), c_amgm, c_l1_magic), 6)
            else:
                l1_m = (nc_cfg["dep_pairs"] if noise == "depolarizing"
                        else sum(nc_cfg["supp_pairs"].values(), Fraction(0)))
                val_m = sum((coeff * c_vec.get(a, Fraction(1)) * c_vec.get(b, Fraction(1))
                             for (a, b), coeff in nc_cfg["supp_pairs"].items()), Fraction(1))
                mx_c = max((sum((coeff * c_vec[a] * c_vec[b] for (a, b), coeff in cfg["supp_pairs"].items()), Fraction(0))
                            / c_vec[r]) for r, cfg in roles_cfg.items())
                c_vec["magic"] = up(max(Fraction(1, 100), l1_m / m_max, (val_m / max(mx_c, Fraction(1))) * Fraction(6, 5)), 6)
    elif isinstance(cw_raw, dict) and set(cw_raw) == set(all_vars):
        c_vec = {k: _pos_rational(cw_raw[k], f"cone_weights[{k}]") for k in all_vars}
    else:
        raise ValueError("cone_weights must be 'auto' or provide one positive weight for every active role")

    # Check if a positive linear floor obstructs asymptotic contraction to zero.
    has_linear_floor = nc_cfg is not None and nc_cfg["linear_floor"] > 0

    def full_ls_unclamped(vec):
        nxt = _eval_ls_clifford(roles_cfg, vec, clamp=False)
        if nc_cfg is not None:
            nxt["magic"] = _eval_magic_step(nc_cfg, vec, nxt, 2, "local-stochastic", clamp=False)
        return nxt

    # Find exact fixed-cone contraction point u* = alpha* * c_vec with Phi_LS(u*) <= (999/1000) u*
    # when linear_floor == 0.
    threshold_lower = Fraction(0)
    contraction_witness = None
    if not has_linear_floor:
        lo, hi = Fraction(0), Fraction(1)
        while max(full_ls_unclamped({k: hi * v for k, v in c_vec.items()})[k] / (hi * c_vec[k]) for k in all_vars) > 1:
            hi /= 2
        # Refine hi upward once by binary search on [hi, 2*hi]
        a_lo, a_hi = hi, hi * 2
        for _ in range(32):
            mid = (a_lo + a_hi) / 2
            u_mid = {k: mid * v for k, v in c_vec.items()}
            img = full_ls_unclamped(u_mid)
            if all(img[k] <= u_mid[k] for k in all_vars) and all(u_mid[k] <= 1 for k in all_vars):
                a_lo = mid
            else:
                a_hi = mid
        # Step slightly inside the basin (999/1000 * a_lo) so strict contraction eta < 1 holds.
        alpha_star = up(a_lo * Fraction(999, 1000), 12)
        while True:
            u_star = {k: alpha_star * c_vec[k] for k, v in c_vec.items()}
            img_u = full_ls_unclamped(u_star)
            if all(img_u[k] < u_star[k] for k in all_vars) and all(u_star[k] <= 1 for k in all_vars):
                break
            alpha_star /= 2
        eta_cone = max(img_u[k] / u_star[k] for k in all_vars)

        # Binary search for the largest physical rate p_th where level-1 monotone bound <= u_star.
        def level1_monotone(p_cand):
            if noise == "local-stochastic":
                p0_vec = {r: p_cand for r in role_names}
                if nc_cfg is not None:
                    p0_vec["magic"] = nc_cfg["inj_mult"] * p_cand
                c1 = _eval_ls_clifford(roles_cfg, p0_vec, clamp=False)
            else:
                c1 = _eval_dep_level1_monotone(roles_cfg, p_cand, clamp=False)
            if nc_cfg is not None:
                p0_vec = {r: p_cand for r in role_names}
                p0_vec["magic"] = nc_cfg["inj_mult"] * p_cand
                c1["magic"] = _eval_magic_step(nc_cfg, p0_vec, c1, 1, noise, p0=p_cand, clamp=False)
            return c1

        p_hi = Fraction(1, 10)
        while any(level1_monotone(p_hi)[k] >= u_star[k] for k in all_vars):
            p_hi /= 2
        p_lo, p_top = p_hi, p_hi * 2
        for _ in range(36):
            pmid = (p_lo + p_top) / 2
            if pmid < 1 and all(level1_monotone(pmid)[k] < u_star[k] for k in all_vars):
                p_lo = pmid
            else:
                p_top = pmid
        # Round threshold DOWNWARD to 8 significant decimal digits so certified lower bound is clean and exact.
        e_dec = ((p_lo.numerator.bit_length() - p_lo.denominator.bit_length()) * 30103) // 100000
        scale = 8 - e_dec
        scaled = p_lo * (Fraction(10) ** scale)
        p_th = Fraction(scaled.numerator // scaled.denominator) / (Fraction(10) ** scale)
        l1_at_pth = level1_monotone(p_th)
        if not all(l1_at_pth[k] < u_star[k] for k in all_vars):
            raise ValueError("internal threshold certification check failed")
        eta_at_pth = max(l1_at_pth[k] / u_star[k] for k in all_vars)
        threshold_lower = p_th
        contraction_witness = dict(
            cone_weights={k: str(c_vec[k]) for k in all_vars},
            alpha_star=str(alpha_star),
            u_star={k: str(u_star[k]) for k in all_vars},
            cone_contraction_ratio=str(up(eta_cone, 12)),
            level1_to_ustar_ratio_at_pth=str(up(eta_at_pth, 12)),
        )

    # Compute exact outward-rounded trajectory p^(0), p^(1), ..., p^(max_levels) at spec["p"].
    p_cur = {r: p for r in role_names}
    if nc_cfg is not None:
        p_cur["magic"] = min(Fraction(1), nc_cfg["inj_mult"] * p)
    levels_trace = [{"level": 0, "rates": {k: decimal(up(p_cur[k], 16)) for k in all_vars},
                     "rates_sci": {k: sci(p_cur[k], 5) for k in all_vars}}]
    exact_levels = [dict(p_cur)]

    for lv in range(1, max_levels + 1):
        if lv == 1 and noise == "depolarizing":
            nxt_c = _eval_dep_level1_pointwise(roles_cfg, p)
        else:
            nxt_c = _eval_ls_clifford(roles_cfg, p_cur, clamp=True)
        nxt = {r: up(nxt_c[r], SIG) for r in role_names}
        if nc_cfg is not None:
            m_val = _eval_magic_step(nc_cfg, p_cur, nxt, lv, noise, p0=p, clamp=True)
            nxt["magic"] = min(Fraction(1), up(m_val, SIG))
        p_cur = nxt
        exact_levels.append(dict(p_cur))
        levels_trace.append({
            "level": lv,
            "rates": {k: decimal(up(p_cur[k], 16)) for k in all_vars},
            "rates_sci": {k: sci(p_cur[k], 5) for k in all_vars},
        })

    workload_report = None
    status = "proved"
    if "workload" in spec and spec["workload"] is not None:
        wl = spec["workload"]
        w_keys = ["logical_counts", "logical_qubits", "target_eps"]
        if "uncompute" in wl:
            w_keys.append("uncompute")
        fields(wl, *w_keys)
        eps = _pos_rational(wl["target_eps"], "workload.target_eps")
        if eps >= 1:
            raise ValueError("workload.target_eps must be in (0, 1)")
        q_log = _pos_int(wl["logical_qubits"], "workload.logical_qubits")
        lcounts = {k: _nonneg_int(v, f"logical_counts[{k}]") for k, v in wl["logical_counts"].items()}
        if not set(lcounts) <= set(all_vars):
            raise ValueError("workload.logical_counts references undeclared role")

        uncomp_info = None
        if "uncompute" in wl:
            import json
            from ftqc_uncompute.certify import RECEIPT_SCHEMA, build_receipt_for_kmx
            from ftqc_uncompute.ir_kickmix import KickMixCircuit

            ub = wl["uncompute"]
            fields(ub, "kmx", "receipt")
            urec = ub["receipt"]
            if urec.get("schema") != RECEIPT_SCHEMA or urec.get("verdict") != "CERTIFIED":
                raise ValueError("workload.uncompute requires a CERTIFIED ftqc_uncompute receipt")
            circuit = KickMixCircuit.from_kmx(ub["kmx"])
            if circuit.sha256() != urec.get("target_sha256"):
                raise ValueError("workload.uncompute kmx hash does not match receipt")
            opts = urec.get("verification_options", {})
            rebuilt = build_receipt_for_kmx(
                circuit,
                source_name=urec.get("source_name", "<memory>"),
                shots=int(opts.get("shots", 64)),
                salt=str(opts.get("salt", "ftqc_uncompute-fs-v1")),
                run_kraus=bool(opts.get("run_kraus", False)),
            )
            if json.dumps(rebuilt, sort_keys=True) != json.dumps(urec, sort_keys=True):
                raise ValueError("workload.uncompute receipt differs from recomputed evidence")
            verified_receipts += 1
            uncomp_info = {
                "target_sha256": urec["target_sha256"],
                "num_qubits": circuit.num_qubits,
                "branch_injective": bool(urec.get("hmr_verification", {}).get("branch_injective")),
            }

        chosen_k = None
        chosen_tv = None
        per_k_tv = []
        for lv, rates_lv in enumerate(exact_levels):
            tv = sum((cnt * rates_lv[k] for k, cnt in lcounts.items()), Fraction(0))
            per_k_tv.append({"level": lv, "tv_upper": decimal(up(tv, 16)), "tv_sci": sci(tv, 5)})
            if chosen_k is None and tv <= eps:
                chosen_k = lv
                chosen_tv = tv

        if chosen_k is None:
            status = "counterexample" if has_linear_floor else "unreached_within_max_levels"
            workload_report = {
                "achieved": False,
                "target_eps": str(eps),
                "logical_qubits": q_log,
                "logical_counts": lcounts,
                "per_level_tv": per_k_tv,
                "uncompute": uncomp_info,
                "reason": ("positive linear error floor blocks asymptotic suppression"
                           if has_linear_floor else f"target_eps not met within max_levels={max_levels}"),
            }
        else:
            # Compute physical qubits and recursive 1-Rec location volume at level chosen_k.
            # Data block qubits scale as 7^k; peak ancilla overhead of an ExRec adds PeakAnc * 7^(k-1).
            data_q = q_log * (7 ** chosen_k)
            anc_q = (max(c["qubits"] - 7 for c in roles_cfg.values()) * (7 ** (chosen_k - 1))
                     if chosen_k >= 1 else 0)
            # Recursive role-count expansion from level chosen_k down to level 0 physical locations.
            cur_counts = dict(lcounts)
            magic_factories_q = 0
            if "magic" in cur_counts and cur_counts["magic"] > 0 and nc_cfg is not None:
                if nc_cfg["scheme"] == "rm15-distill" and chosen_k >= 1:
                    magic_factories_q = (nc_cfg["total_slots"] ** chosen_k) * 7
            # Expand magic at level chosen_k down through levels chosen_k .. 1:
            role_vec = {r: cur_counts.get(r, 0) for r in role_names}
            magic_need = cur_counts.get("magic", 0)
            raw_inj_locs = 0
            for lv in range(chosen_k, 0, -1):
                if magic_need and nc_cfg is not None:
                    if nc_cfg["scheme"] == "rm15-distill":
                        for r, c_cnt in nc_cfg["clifford_counts"].items():
                            role_vec[r] += magic_need * c_cnt
                        magic_need *= nc_cfg["total_slots"]
                    else:
                        for r, c_cnt in nc_cfg["role_locs"].items():
                            if r in role_vec:
                                role_vec[r] += magic_need * c_cnt
                # Expand one concatenation step lv -> lv - 1 using one_rec_role_counts or role_locs
                nxt_vec = {r: 0 for r in role_names}
                for r, cnt in role_vec.items():
                    if not cnt:
                        continue
                    sub_map = roles_cfg[r]["one_rec_role_counts"] or roles_cfg[r]["role_locs"]
                    for sub_r, sub_cnt in sub_map.items():
                        nxt_vec[sub_r] += cnt * sub_cnt
                role_vec = nxt_vec
            if magic_need:
                # Base physical T-state preparations at level 0
                raw_inj_locs = magic_need * 4
            phys_volume = sum(role_vec.values()) + raw_inj_locs
            workload_report = {
                "achieved": True,
                "chosen_level": chosen_k,
                "target_eps": str(eps),
                "tv_upper": decimal(up(chosen_tv, 16)),
                "tv_sci": sci(chosen_tv, 5),
                "logical_qubits": q_log,
                "logical_counts": lcounts,
                "physical_data_qubits": data_q,
                "physical_ancilla_qubits": anc_q,
                "distillation_factory_qubits": magic_factories_q,
                "total_physical_qubits": data_q + anc_q + magic_factories_q,
                "physical_fault_volume": phys_volume,
                "level0_role_counts": role_vec,
                "level0_magic_injections": magic_need,
                "per_level_tv": per_k_tv,
                "uncompute": uncomp_info,
            }

    return {
        "status": status,
        "rule": "multivariate extended-rectangle renormalization with non-overlapping truncated regions",
        "noise": noise,
        "p": str(p),
        "below_certified_threshold": bool(not has_linear_floor and p <= threshold_lower),
        "threshold_lower": str(threshold_lower),
        "threshold_lower_sci": sci(threshold_lower, 5) if threshold_lower > 0 else "0",
        "has_linear_floor": has_linear_floor,
        "contraction_certificate": contraction_witness,
        "verified_subreceipts": verified_receipts,
        "levels": levels_trace,
        "workload": workload_report,
        "premises": [],
    }
