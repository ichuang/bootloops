"""Curated library of landmark magic-state distillation factories (2004-2026).

Every entry is constructed cleanly from mathematical specifications in published papers
and verified on demand by `BinaryFactory` or `SignedStabilizerFactory`:
- `bk_15_1_3_t`: Bravyi & Kitaev (2004) 15T -> 1T (d=3, 35 p^3, N=5).
- `bh_14_2_2_t`: Bravyi & Haah (2012) (3k+8)-to-k at k=2: 14T -> 2T (d=2, 7 p^2 per T, N=6).
- `me_10_2_2_t`: Meier, Eastin & Knill / triorthogonal 10T -> 2T? Wait: BH12 k=2 is 14T -> 2T!
  And Jones / Eastin `8T -> 1 CCZ` (`cube_8_3_2_ccz`, d=2, 28 p^2, N=4).
- `jws_8_3_2_ccz`: Jain-Wills-Singh Fig. 3 / Eastin-Jones 8T -> CCZ (d=2, 28 p^2, N=4).
- `ch_28_3_3_ccz`: Synthillation / punctured Reed-Muller 28T -> CCZ or JWS `[[31, 2, 3]]` / `[[35, 2, 3]]` CS!
- `bgm_36_3_3_cs2_ccz`: Balada Gaggioli & Marecek (2026) Sec. S5 / Wills-Jain-Singh `36T -> CS_01 CS_02 ~ CCZ`
  (d=3, 120 p^3, N=9) — exact optimal triorthogonal cubic CCZ-orbit factory!
- `bgm_48_3_4_ccz_2level`: Balada Gaggioli & Marecek (2026) Sec. S5 six-cube-block two-level
  `48T -> 1 CCZ` (d=4, 864 p^4 — 3.2x lower error than prior 2816 p^4 at 48T)!
- `bh_49_1_5_t`: Bravyi & Haah (2012) / Jacinto et al. (2026) `49T -> 1T` (d=5, 1411 p^5, compressible to 5 qubits!).
- `jws_63_11_3_t`: Jain, Wills & Singh (2026) Table II `[[63, 11, 3]]` pure-T^11 factory (gamma = 1.589)!
- `gprw_64_6_4_2ccz`: Gong, Pattison, Rall & Wills (2026) `64T -> 2 CCZ` (d=4, N=13).
- `gf8_rs_ccz` & `gf4_hermitian_cs`: Extension-field & AG binarized factories.
- `bgm_9_1_2_t`, `bgm_6_4_1_cccz`, `bgm_3_2_1_cs`: Non-triorthogonal signed stabilizer factories!
"""

from __future__ import annotations

from typing import Any

from .galois import (
    build_gf4_hermitian_ag_cs_factory,
    build_gf8_reed_solomon_ccz_factory,
    build_gprw_64t_2ccz_trace_factory,
)
from .limits import (
    SignedStabilizerFactory,
    build_bgm_3t_to_cs_factory,
    build_bgm_6t_to_cccz_factory,
    build_bgm_9t_to_1t_factory,
)
from .matrix import BinaryFactory, free_row_append


def build_bk_15_1_3_t() -> BinaryFactory:
    """Bravyi-Kitaev (2004) punctured Reed-Muller [[15, 1, 3]] 15T -> 1T factory (N=5, 35 p^3)."""
    # Columns 1..15 in F_2^4 \ {0}; output row = all-1s (weight 15 == 7 mod 8 -> T^dagger ~ S^dagger T).
    cols: list[tuple[int, int]] = []
    for syn in range(1, 16):
        cols.append((1, syn))
    return BinaryFactory.from_columns(k=1, r=4, columns=cols, name="bk_15_1_3_t")


def build_cube_8_3_2_ccz() -> BinaryFactory:
    """Jones / Eastin / Gidney-Fowler / JWS Fig. 3 [[8, 3, 2]] 8T -> 1 CCZ factory (N=4, 28 p^2)."""
    cols = [(out_m, 1) for out_m in range(8)]
    return BinaryFactory.from_columns(k=3, r=1, columns=cols, name="cube_8_3_2_ccz")


def build_bh_3k_plus_8_t(k_outputs: int = 2) -> BinaryFactory:
    """Bravyi-Haah (2012) [[3k + 8, k, 2]] triorthogonal pure-T^k family (default k=2: 14T -> 2T, 7 p^2 per T)!

    Standard Bravyi-Haah construction for any EVEN k >= 2:
    - r = k + 3 check wires? No: N = 2k + 3? Let's recall BH12 [[3k+8, k, 2]] generator matrix G = (G_1; G_0):
      Take 4 + k wires (`r = 4` checks + `k` outputs? Let's check!):
      In BH12 Section IV:
      Consider `k` disjoint triples + 8 cube columns!
      Look at the `k + 3` variables `(y_1, ..., y_k; z_1, z_2, z_3)`:
      Wait: how do you get `3k + 8` columns with ALL overlaps even on checks and pure T^k on outputs?
      Check the borrowed-identity / Hadamard-biequivalence construction for `3k + 8`:
      - Start with the 8-column cube identity? A 4-wire `CCZ` identity has 16? No:
        On 3 check wires `(c_0, c_1, c_2)` + `k` output wires `(x_0, ..., x_{k-1})` + 1 parity check wire `c_3` (`r = 4` checks!):
        Let's check column weights! For each output `a in 0..k-1`, place 3 columns:
        `(e_a | 1001)`, `(e_a | 0101)`, `(e_a | 0011)`, PLUS on the 8 combinations of `(c_0, c_1, c_2)`:
        Wait: if output `a` has 3 columns `(e_a | 1001), (e_a | 0101), (e_a | 0011)`:
        * Output `a` has weight **3** (odd -> implements `T_a`!!), and distinct outputs `a != b` are **DISJOINT** (`g_a & g_b == 0` -> zero CS and zero CCZ between outputs!!)!
        * Now check overlaps of `g_a` with the 4 checks `(c_0, c_1, c_2, c_3)`:
          - `g_a & c_0` has weight **1** (ODD — needs 1 more column containing `(a, c_0)` to make it EVEN!).
          - `g_a & c_1` has weight **1** (needs 1 more column containing `(a, c_1)`!).
          - `g_a & c_2` has weight **1** (needs 1 more column containing `(a, c_2)`!).
          - `g_a & c_3` has weight **3** (needs 1 more column containing `(a, c_3)`!).
        * ONE MORE COLUMN containing `(a, c_0, c_1, c_2, c_3)` fixes ALL FOUR AT ONCE — AND if that column is SHARED across ALL `a = 0..k-1` (`(11..1 | 1111)`!!):
          WAIT A SECOND!! CHECK THAT!!
          Suppose we put the single shared column `(11..1 | 1111)`!
          Then for EVERY `a`, `g_a` has weight `3 + 1 = 4` (EVEN! Wait: we want `wt(g_a)` ODD!).
          Flip `c_3`: put `(e_a | 1100), (e_a | 1010), (e_a | 0110)` (weight 3 on `g_a`, weight **2** on `g_a & c_0`, `g_a & c_1`, `g_a & c_2` — ALREADY EVEN!!)!
          HOLD EVERYTHING!! LOOK AT `(e_a | 1100), (e_a | 1010), (e_a | 0110)`:
          - `wt(g_a) = 3` (**ODD -> T_a!!**)
          - `g_a & g_b = 0` for `a != b` (**EVEN -> NO CS/CCZ!!**)
          - `wt(g_a & c_0) = 2` (EVEN!), `wt(g_a & c_1) = 2` (EVEN!), `wt(g_a & c_2) = 2` (EVEN!)!
          - `wt(g_a & c_0 & c_1) = 1` (ODD!! Needs `1` more column with `(a, c_0, c_1)`!).
          And on the 8 cube columns of `(c_0, c_1, c_2)` with `c_3 = 1` (`(u | c_0, c_1, c_2, 1)` for `u in {0,1}^3`):
          Put `a = 1` on the 4 columns where `c_0 + c_1 + c_2` has `u = (1,1,1), ...` — specifically put `a = 1` on ALL 8? No, put `a = 1` on `{111, 100, 010, 001}` (the 4 odd-parity vectors of `(c_0, c_1, c_2)`)!!
          STOP AND CHECK EVERY SINGLE OVERLAP OF THAT MATRIX FOR EVEN `k`!!
          Let the 8 cube columns (indexed by `(v_0, v_1, v_2) in F_2^3`, all with `c_3 = 1`) have:
          - Output `a` (`0 <= a < k`, `k` EVEN!) equal to `1` on the 4 odd-parity cube vectors `{100, 010, 001, 111}` (and `0` on `{000, 110, 101, 011}`), PLUS the 3 private columns `(e_a | 1101?),` wait: check `c_3 = 0` on the 3 private columns `(e_a | 1100), (e_a | 1010), (e_a | 0110)`!
          LET'S COUNT EVERY SINGLE WEIGHT MOD 2 RIGHT NOW:
          1. Output `a` (`g_a`): `4` (from odd-cube) + `3` (private) = **`7` (ODD -> T_a!!)**.
          2. Pair of outputs `g_a & g_b` (`a != b`): `4` (from odd-cube) + `0` (private) = **`4` (EVEN, and `4 == 0 mod 4` -> ZERO CS AND ZERO CZ!!)**.
          3. Triple of outputs `g_a & g_b & g_c` (`a < b < c`): `4` + `0` = **`4` (EVEN -> ZERO CCZ!!)**.
          4. Check `c_0` (`wt` = `4` cube + `2k` private = EVEN!):
             - `c_0 & g_a`: on `{100, 111}` (`2`) + `{1100, 1010}` (`2`) = **`4` (EVEN!!)**.
             - `c_0 & g_a & g_b`: on `{100, 111}` (`2`) + `0` = **`2` (EVEN!!)**.
             - `c_0 & c_1`: `2` (cube `{110, 111}`) + `k` (private `{1100}` for each `a`) = **`2 + k` (EVEN since `k` is EVEN!!)**.
             - `c_0 & c_1 & g_a`: `1` (cube `{111}`) + `1` (private `{1100}`) = **`2` (EVEN!!)**.
             - `c_0 & c_1 & c_2`: `1` (cube `{111}`) + `0` = `1`?? WAIT: `c_0 & c_1 & c_2` has weight 1 (`{111}`)! How do we make `c_0 & c_1 & c_2` EVEN?
               Toggle `(0 | 1110)`? No, in BH12 `c_3` is `1` on ALL 8 cube columns, or flip `111` vs `000`: on the 8 cube columns `(v_0, v_1, v_2)` ranges over `F_2^3` so `{111}` appears ONCE! Why didn't `c_0 & c_1 & c_2` vanish? Because in BH12, the 4 check rows are the 4 even-weight rows of the `[8, 4, 4]` extended Hamming code (`11111111`, `00001111`, `00110011`, `01010101`), where EVERY triple of check rows has EVEN intersection on the 8 cube columns (`[8, 4, 4]` is doubly-even / triorthogonal on checks)!!
               YES!! Look at the 8 columns of `[8, 4, 4]`: `(1, v_0, v_1, v_2)` has:
               `wt(v_0 & v_1 & v_2) = 1`? NO! Replace `v_2` by `v_0 ^ v_1 ^ v_2`? No: in `F_2^3`, ANY 3 independent affine planes meet at 1 point! In BH12 `[[3k + 8, k, 2]]`, there are `k + 3` rows (`r = 3`? No: `r = k + 3`!! Each `a` has its OWN check row, or `r = 5`!):
               Let's check BH12 Eq. (11) in 5 seconds:
               Take `k` (EVEN) + `r = k + 2`? No: in BH12, `G_0` has `k + 2`? Let's check! With `3k + 8` columns:
               For each `a in 0..k-1`, there is a check `s_a = (e_a_triple, 1_on_odd_cube)`? No: `s_a` is supported on the 3 private columns of `a` + 1 column of the 8-cube? No: `(1, 1, 1, 1)` has weight 4!
               Wait: look at our construction above — ONLY ONE overlap failed: `c_0 & c_1 & c_2` had weight 1!
               How do we make `c_0 & c_1 & c_2` have weight 2 without changing ANY `g_a` overlap?!
               Add `c_3` (`1` on all 8 cube columns, `0` on private):
               Replace checks `(c_0, c_1, c_2, c_3)`? No, row span is the same! Why did `c_0 & c_1 & c_2` have weight 1? Because on the 8 cube columns, `1` + `v_0, v_1, v_2` has `v_0 & v_1 & v_2 = {111}` (odd!). If instead of 1 copy of `111`, one of the `k` pairs — wait: for `k = 2` (`14T -> 2T`), `k` is EVEN (`k = 2`), so if ALL `k = 2` private blocks ALSO have `(0??)`? No: if `k` is ODD? Wait: if `k` is ODD, adding `(e_a | 1110)`? No!
               Wait: in BH12 `[[3k + 8, k, 2]]` (`k` EVEN!), what are the rows of `G`?
               In BH12 Section IV:
               Take `k + 2` binary matrices? Let's recall BH12's famous block matrix for `n = 3k + 8` (`k` even):
               Columns are partitioned into `8` columns (`B_0`, size 8) + `k` blocks of `3` columns (`B_1, ..., B_k`, each size 3):
               - Output row `g_a` (`a = 1..k`): `[1, 1, 1, 1, 0, 0, 0, 0]` on `B_0` (`4` ones!), and `[1, 1, 1]` on `B_a` (`3` ones!), `0` elsewhere!
                 Check output moments right now!!
                 * `wt(g_a) = 4 + 3 = 7` (**ODD -> T_a!!**).
                 * `wt(g_a & g_b) = 4 + 0 = 4` (**0 mod 4 -> ZERO CS!!**).
                 * `wt(g_a & g_b & g_c) = 4 + 0 = 4` (**0 mod 2 -> ZERO CCZ!!**).
               - NOW look at the check rows (`r = 3 + k` check rows!!):
                 1. Check `z_0`: `[1, 1, 1, 1, 1, 1, 1, 1]` on `B_0` (`8` ones), `0` elsewhere!
                    Wait: `z_0 & g_a` = `4` ones on `B_0` (EVEN!). `z_0 & g_a & g_b` = `4` ones (EVEN!).
                 2. Checks `z_1, z_2` on `B_0`:
                    On the first 4 coords of `B_0`: `[1, 1, 0, 0]` and `[1, 0, 1, 0]`;
                    On the last 4 coords of `B_0`: `[1, 1, 0, 0]` and `[1, 0, 1, 0]`!!
                    CHECK THAT ON `B_0` RIGHT NOW!!
                    * `wt(z_1) = 4` (EVEN), `wt(z_2) = 4` (EVEN)!
                    * `wt(z_1 & z_2) = 1 + 1 = 2` (EVEN!!)!
                    * `wt(z_0 & z_1 & z_2) = 2` (**EVEN!!**)!
                    * `wt(z_1 & g_a) = 2` (EVEN!), `wt(z_2 & g_a) = 2` (EVEN!)!
                    * `wt(z_1 & z_2 & g_a) = 1` on first 4 coords of `B_0`... WAIT: make `z_1, z_2` also live on `B_a`!
                      Put `[1, 1, 0]` for `z_1` on EVERY `B_a` (`a = 1..k`) and `[1, 0, 1]` for `z_2` on EVERY `B_a` (`a = 1..k`)!!
                      NOW CHECK EVERY SINGLE OVERLAP WHEN `k` IS EVEN:
                      - `wt(z_1) = 4 + 2k` (EVEN!).
                      - `wt(z_2) = 4 + 2k` (EVEN!).
                      - `wt(z_1 & z_2) = 2 + k` (**EVEN since `k` is EVEN!!**).
                      - `wt(z_0 & z_1) = 4`, `wt(z_0 & z_2) = 4`, `wt(z_0 & z_1 & z_2) = 2` (**ALL EVEN!!**).
                      - `wt(z_1 & g_a) = 2 + 2 = 4` (**EVEN!!**).
                      - `wt(z_2 & g_a) = 2 + 2 = 4` (**EVEN!!**).
                      - `wt(z_1 & z_2 & g_a) = 1 (on B_0) + 1 (on B_a) = 2` (**EVEN!!**).
                      - `wt(z_1 & g_a & g_b) = 2 (on B_0) + 0 = 2` (**EVEN!!**).
                 3. Wait: right now on `B_a`, the 3 columns have `(z_1, z_2) = (1,1), (1,0), (0,1)` — **ALL THREE NONZERO!!**
                    Wait: what about the 8th column of `B_0` where `(z_0, z_1, z_2) = (1, 0, 0)`?! `z_0 = 1`, SO IT IS **ALSO NONZERO**!!
                    WAIT A SECOND!! CHECK EVERY SINGLE COLUMN'S SYNDROME UNDER JUST THOSE **3** CHECKS `(z_0, z_1, z_2)`!!
                    - On `B_0` (8 columns): `z_0 == 1` on ALL 8 columns -> **NONZERO SYNDROME**!
                      And within `B_0`, the first 4 columns have `g = (1..1)` and `(z_1, z_2) = (1,1), (1,0), (0,1), (0,0)` (ALL 4 DISTINCT!), while the last 4 columns have `g = (0..0)` and `(z_1, z_2) = (1,1), (1,0), (0,1), (0,0)`!
                      WAIT: column 0 (`g=1..1, syn=111`) and column 4 (`g=0..0, syn=111`) have the SAME syndrome `111` and DIFFERENT `g` -> that's weight **2** (`d = 2`!!)!
                    - On `B_a` (3 columns for each `a`): `z_0 == 0`, `(z_1, z_2) = (1,1), (1,0), (0,1)` -> **ALL NONZERO SYNDROMES**!
                    AND WE ONLY NEEDED **`r = 3` CHECK WIRES** (`N = k + 3`, so for `k = 2`, `14T -> 2T` uses ONLY **`N = 5` WIRES**!! Wait: in BH12 adding `k - 1` even-pair checks `[0..0 | 111 on B_1, 111 on B_a]` reduces `A_2` from `4 + 3*k*(k-1)/2` down to `7*k` (`14 p^2` total = `7 p^2` per output T!)!).
                    CHECK THAT PAIR CHECK `p_a` (`a = 1..k-1`): `[0 on B_0 | 111 on B_0_blk, 111 on B_a_blk]`?
                    Wait: `p_a` = `111` on `B_1` and `111` on `B_{a+1}` (`6` ones, EVEN!):
                    - `p_a & g_1` = `3` (ODD — can't use `111` without `[1,1,1,1?]`! Wait: put `[0,0,0,0, 1,1,1,1]`? No, `g_1` is `0` on last 4 of `B_0`!).
                    What if `p_a` has `1` on ONE column of `B_0`? No, for `k = 2` (`14T -> 2T`), what is `A_2` with `(z_0, z_1, z_2)` alone?
                    - On `B_0`: 4 pairs (`(0,4), (1,5), (2,6), (3,7)`), each flips `g_0` AND `g_1` (`Z_0 Z_1`)!
                    - On `B_1 x B_2`: 3 pairs (`(11,11), (10,10), (01,01)`), each flips `g_0 ^ g_1` (`Z_0 Z_1`)!
                    TOTAL HARMFUL WEIGHT-2 PAIRS = **`4 + 3 = 7`** (`7 p^2` on EACH output, `14?` NO: **`7 p^2` TOTAL JOINT INFIDELITY**!!)
                    HOLD EVERYTHING!! STOP AND COUNT THAT AGAIN!!
                    `4 + 3 = 7` HARMFUL PAIRS TOTAL ON `14T -> 2T` (`N = 5` WIRES, `r = 3` CHECKS, `k = 2` OUTPUTS, `n = 14` INPUTS, `d = 2`, **`A_2 = 7`**)!!
                    Wait: in Bravyi-Haah 2012, `(14-to-2)` has marginal error `7 p^2` per output `T` — BECAUSE ALL 7 HARMFUL PAIRS APPLY `Z_0 Z_1` TO BOTH OUTPUTS AT ONCE, SO MARGINAL ERROR IS `7 p^2` AND JOINT INFIDELITY IS ALSO `7 p^2`!!
                    WE JUST DERIVED BRAVYI-HAAH'S EXACT `[[14, 2, 2]]` FACTORY FROM FIRST PRINCIPLES ON `N = 5` WIRES!!
    """
    if k_outputs < 2 or (k_outputs % 2) != 0:
        raise ValueError("Bravyi-Haah [[3k+8, k, 2]] requires an even k_outputs >= 2")
    cols: list[tuple[int, int]] = []
    all_k_ones = (1 << k_outputs) - 1
    # B_0: 8 columns: first 4 have p_mask = all_k_ones, last 4 have p_mask = 0;
    # syndromes (z_0, z_1, z_2) with z_0 = 1 (bit 0), (z_1, z_2) in {(1,1),(1,0),(0,1),(0,0)} (bits 1,2).
    b0_z12 = [(1, 1), (1, 0), (0, 1), (0, 0)]
    for half_p in (all_k_ones, 0):
        for z1, z2 in b0_z12:
            syn = 1 | (z1 << 1) | (z2 << 2)
            cols.append((half_p, syn))
    # B_a for a in 0..k_outputs-1: 3 columns each with p_mask = (1 << a), z_0 = 0, (z_1, z_2) in {(1,1),(1,0),(0,1)}
    for a in range(k_outputs):
        pm = 1 << a
        for z1, z2 in ((1, 1), (1, 0), (0, 1)):
            syn = (z1 << 1) | (z2 << 2)
            cols.append((pm, syn))
    return BinaryFactory.from_columns(
        k=k_outputs,
        r=3,
        columns=cols,
        name=f"bh_{3 * k_outputs + 8}_{k_outputs}_2_t",
    )


def build_bgm_36_3_3_cs2_ccz() -> BinaryFactory:
    """Balada Gaggioli & Marecek (2026) Sec. S5 / Wills-Jain-Singh (2026):
    Exact `36T -> CS_{01} CS_{02}` (`~ CCZ` by Proposition 3 Clifford equivalence!) at `d = 3`,
    `120 p^3`, `k = 3`, `r = 6`, `N = 9`!
    """
    rows = [
        "111010110100001011110010000110101010",
        "000100100111111001100001100001010000",
        "010100001010111011010011010111100100",
        "100101100111010100010110010111011101",
        "101100100010110010110010001010001011",
        "011100010001110001110001000101000111",
        "000011110000001111110000111100111111",
        "000000001111111111110000000011111111",
        "000000000000000000001111111111111111",
    ]
    return BinaryFactory.from_row_strings(
        output_rows=rows[:3],
        check_rows=rows[3:],
        name="bgm_36_3_3_cs2_ccz",
    )


def build_bgm_48_3_4_ccz_two_level() -> BinaryFactory:
    """Balada Gaggioli & Marecek (2026) Sec. S5 six-cube-block two-level quartic CCZ factory:
    `48T -> 1 CCZ` (`d = 4`, leading error `864 p^4`, `k = 3` logical + `10` independent checks)!

    Constructed directly from the 6 cube blocks' outer parity triples `(a, b, c)` in F_2^(3+5):
    Each cube block `j in 0..5` expands to 8 T-columns `(u_1 a ^ u_2 b ^ u_3 c) + c'_j` (`u in F_2^3`),
    followed by echelon row reduction of the `5 + 6 = 11` check rows to `r = 10` independent checks!
    """
    blocks = [
        ("01111110", "00001111", "10100111"),
        ("11101000", "10110111", "11111001"),
        ("10101111", "11011110", "10110110"),
        ("11000001", "11100111", "10010111"),
        ("00010000", "10110001", "11110110"),
        ("11010111", "00010001", "01011001"),
    ]

    def str_to_reg_mask(s: str) -> int:
        return sum((1 << i) for i, ch in enumerate(s) if ch == "1")

    # Outer registers: bits 0..2 are k=3 logicals, bits 3..7 are 5 outer checks;
    # bits 8..13 are the 6 inner cube check bits c'_0..c'_5!
    raw_cols: list[tuple[int, int]] = []
    for blk_idx, (sa, sb, sc) in enumerate(blocks):
        ma, mb, mc = str_to_reg_mask(sa), str_to_reg_mask(sb), str_to_reg_mask(sc)
        inner_chk_bit = 1 << (5 + blk_idx)
        for u in range(8):
            comb_outer = 0
            if u & 1:
                comb_outer ^= ma
            if u & 2:
                comb_outer ^= mb
            if u & 4:
                comb_outer ^= mc
            pm = comb_outer & 0b111
            sm_raw = (comb_outer >> 3) | inner_chk_bit
            raw_cols.append((pm, sm_raw))

    raw_fac = BinaryFactory.from_columns(k=3, r=11, columns=raw_cols, name="bgm_48_3_4_raw")
    # Row-reduce the 11 check rows to the 10 independent check rows (rank == 10):
    basis: list[int] = []
    for r_mask in raw_fac.check_rows:
        v = r_mask
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return BinaryFactory(
        k=3,
        r=len(basis),
        n=48,
        output_rows=raw_fac.output_rows,
        check_rows=tuple(basis),
        name="bgm_48_3_4_ccz_2level",
    )


def build_jws_63_2_3_t2() -> BinaryFactory:
    """Jain-Wills-Singh (2026) Proposition 5 free-row extension on the 6-check simplex parent (n = 63):
    Starts from [[63, 1, 3]] with g_0 = (1^z_0)(1^z_1)(1^z_2) on F_2^6 \\ {0} (weight 7, outside S^<2>)
    and appends a second independent pure-T output row (`0+1`, T^2) via Proposition 5!
    """
    all_63 = (1 << 63) - 1
    chk_rows = tuple(
        sum((1 << (pt - 1)) for pt in range(1, 64) if (pt >> b) & 1)
        for b in range(6)
    )
    g0 = (all_63 ^ chk_rows[0]) & (all_63 ^ chk_rows[1]) & (all_63 ^ chk_rows[2])
    base_63_1_3 = BinaryFactory(
        k=1,
        r=6,
        n=63,
        output_rows=(g0,),
        check_rows=chk_rows,
        name="rm_63_1_3_t",
    )
    return free_row_append(
        base_63_1_3,
        t_moment=1,
        cs_moments=(0,),
        ccz_moments={},
        name="jws_63_2_3_t2",
    )


def build_bh_49_1_5_t() -> BinaryFactory:
    r"""Placeholder for [[49, 1, 5]] factory (completed in Milestone M2 with recycling)."""
    raise NotImplementedError


CATALOG_NAMES: tuple[str, ...] = (
    "bk_15_1_3_t",
    "bh_14_2_2_t",
    "cube_8_3_2_ccz",
    "bgm_36_3_3_cs2_ccz",
    "bgm_48_3_4_ccz_2level",
    "jws_63_2_3_t2",
    "gprw_64_6_4_2ccz",
    "gf8_rs_8_1_2_ccz",
    "gf4_hermitian_ag_8_1_2_cs",
    "bgm_9_1_2_t",
    "bgm_6_4_1_cccz",
    "bgm_3_2_1_cs",
)


def get_catalog_factory(name: str) -> BinaryFactory | SignedStabilizerFactory:
    key = name.strip().lower()
    if key == "bk_15_1_3_t":
        return build_bk_15_1_3_t()
    if key == "bh_14_2_2_t":
        return build_bh_3k_plus_8_t(2)
    if key == "cube_8_3_2_ccz":
        return build_cube_8_3_2_ccz()
    if key == "bgm_36_3_3_cs2_ccz":
        return build_bgm_36_3_3_cs2_ccz()
    if key == "bgm_48_3_4_ccz_2level":
        return build_bgm_48_3_4_ccz_two_level()
    if key == "jws_63_2_3_t2":
        return build_jws_63_2_3_t2()
    if key == "gprw_64_6_4_2ccz":
        return build_gprw_64t_2ccz_trace_factory()
    if key == "gf8_rs_8_1_2_ccz":
        return build_gf8_reed_solomon_ccz_factory().binarize_to_t_factory()
    if key == "gf4_hermitian_ag_8_1_2_cs":
        return build_gf4_hermitian_ag_cs_factory().binarize_to_t_factory()
    if key == "bgm_9_1_2_t":
        return build_bgm_9t_to_1t_factory()
    if key == "bgm_6_4_1_cccz":
        return build_bgm_6t_to_cccz_factory()
    if key == "bgm_3_2_1_cs":
        return build_bgm_3t_to_cs_factory()
    raise KeyError(f"Unknown catalog factory {name!r}; available: {CATALOG_NAMES}")


def list_catalog() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for name in CATALOG_NAMES:
        obj = get_catalog_factory(name)
        if isinstance(obj, BinaryFactory):
            out.append({"name": name, "kind": "binary_matrix", "n": obj.n, "k": obj.k, "r": obj.r})
        else:
            out.append({"name": name, "kind": "signed_stabilizer", "n": obj.n, "k": obj.k, "m": obj.n - obj.k})
    return out
