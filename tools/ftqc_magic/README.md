# ftqc_magic

`ftqc_magic` verifies, optimizes, and composes magic-state distillation factories, magic-state cultivation schedules, open-boundary escape postselection, and single-shot qLDPC magic factories into replayable, zero-float JSON certificates (`ftqc_magic-receipt-1`).

All calculations run in exact cyclotomic ($\mathbb{Z}[\omega_8]$, $\mathbb{Z}[\sqrt{2}]$), field ($\mathbb{Q}(\sqrt{2})$, $\mathrm{GF}(2^s)$), and rational (`fractions.Fraction`) arithmetic with outward-rounded `Decimal` bounds (`ROUND_CEILING` / `ROUND_FLOOR`) at serialization boundaries. Optional external engines (**CaDiCaL** for LRAT optimality/no-go proofs and **Stim** for Clifford-proxy cross-checks) are probed dynamically when installed.

## Quickstart

From the BootLoops repository root:

```sh
export PYTHONPATH=tools
work=$(mktemp -d)

# 1. Run the cold-clone self-certification, refusal, and mutation battery
python3 tools/ftqc_magic/selftest.py

# 2. Compress Bravyi-Haah 49T -> 1T (d=5, 1411 p^5) from 14 logical wires down to 5 active qubits
python3 -m ftqc_magic check tools/ftqc_magic/examples/bh49_recycled_5qubits.json -o "$work/bh49.receipt.json"
python3 -m ftqc_magic verify tools/ftqc_magic/examples/bh49_recycled_5qubits.json "$work/bh49.receipt.json"

# 3. Certify the optimal non-triorthogonal [[9,1,2]] 9T -> 1T signed stabilizer factory over Z[omega_8]
python3 -m ftqc_magic check tools/ftqc_magic/examples/bgm_9_to_1_nontriorthogonal.json -o "$work/bgm9.receipt.json"

# 4. Run exact O(n^3) Clifford-hook fault enumeration on Z-flagged MSC-3 cultivation (d_T = 3)
python3 -m ftqc_magic cultivate --variant msc3_z_flagged --p 1/1000

# 5. Chain Z-flagged cultivation + Caliper open-boundary escape into a recycled Litinski 15-to-1 factory
python3 -m ftqc_magic pipeline tools/ftqc_magic/examples/flagged_cultivation_caliper_15to1.json \
  -o "$work/cult_15to1.receipt.json"

# 6. Query Pareto-optimal hybrid pipelines reaching 1e-10 conditional infidelity at p = 1/1000
python3 -m ftqc_magic suggest --target T --p 1/1000 --infidelity 1/10000000000
```

## Included Engines & Verifiers

| Capability | Module | Primary References |
|---|---|---|
| Exact $\mathbb{Z}[\omega_8]$, $\mathbb{Q}(\sqrt{2})$, `PolyQ`, & outward `Decimal` bounds | `rings.py`, `ir.py` | Clean-room cyclotomic & interval arithmetic |
| Binary / PPR / borrowed-identity factories, parity lift, puncturing, free-row extension | `matrix.py`, `catalog.py` | [Bravyi & Haah 2012](https://arxiv.org/abs/1209.2426) · [Campbell & Howard 2017](https://arxiv.org/abs/1606.01904) · [Jain, Wills & Singh 2026](https://arxiv.org/abs/2610.06535) |
| $\mathrm{GF}(2^s)$ self-dual trace binarization, Hermitian AG & 3-twirl concatenation | `galois.py` | [Gong, Pattison, Rall & Wills 2026](https://arxiv.org/abs/2610.06384) |
| Non-triorthogonal $\mathbb{Z}[\omega_8]$ stabilizer-collapse verifier & $n_{\min}(\tau, d)$ bounds | `limits.py` | [Balada Gaggioli & Mareček 2026](https://arxiv.org/abs/2610.31115) · [Beverland et al. 2020](https://arxiv.org/abs/1910.11827) |
| Mid-circuit check recycling ($15\to 1$ in 4 qubits, $49\to 1$ in 5 qubits) + LRAT lower bounds | `recycle.py`, `search.py` | [Jacinto et al. 2026](https://arxiv.org/abs/2606.07734) |
| Exact $O(n^3)$ Clifford-hook cultivation fault distance ($d_T$ vs $d_S$) & $Z$-flag restoration | `cultivate.py` | [Gidney, Shutty & Jones 2024](https://arxiv.org/abs/2409.17595) · [Chan, Strikis, Sun & Cai 2026](https://arxiv.org/abs/2609.17706) |
| Open-boundary mid-circuit escape postselection & Caliper worst-case soft-gap verifier | `caliper.py` | [Mehta et al. 2026](https://arxiv.org/abs/2610.06828) |
| Single-shot 3D tricycle $\to$ 4D quadcycle qLDPC `CCZ` factories & AOD shuttle budgets | `qldpc_factory.py` | [Menon et al. 2026](https://arxiv.org/abs/2610.10731) |
| Litinski asymmetric sizing, Gidney–Fowler `C2T` catalysis, & fail-closed pipeline composer | `pipeline.py`, `bridges/` | [Gidney & Fowler 2019](https://arxiv.org/abs/1812.01238) · [Litinski 2019](https://arxiv.org/abs/1905.06903) |

See [GUIDE.md](GUIDE.md) for CLI subcommands, receipt schema, fail-closed correlation/catalysis hazards (`N1`–`N6`), and full bibliography, and [../FTQC.md](../FTQC.md) for how `ftqc_magic` connects to `ftqc_codes`, `ftqc_faults`, `ftqc_decoders`, `ftqc_extractors`, `ftqc_compose`, and `ftqc_uncompute`.
