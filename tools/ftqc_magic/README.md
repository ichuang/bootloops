# ftqc_magic

Exact statements about magic-state factories, cultivation schedules,
open-boundary escape stages, qLDPC CCZ factories and the pipelines that chain
them, as `ftqc_magic-receipt-1` receipts. Distillation is analysed in the
standard model (independent input errors, ideal Clifford checks); cultivation
on its final check subcircuit with an ideal encoded input. All arithmetic is
exact (Z[ω₈], Z[√2], Q(√2), GF(2^s), rationals), with outward-rounded
decimals. CaDiCaL (LRAT proofs) and Stim (cross-checks) are optional.

## Quickstart

From the BootLoops repository root:

```sh
export PYTHONPATH=tools
work=$(mktemp -d)

# 1. Run the battery (about 8 s)
python3 tools/ftqc_magic/selftest.py

# 2. Compress Bravyi-Haah 49T -> 1T (d=5, 1411 p^5) from 14 logical wires down to 5 active qubits
python3 -m ftqc_magic check tools/ftqc_magic/examples/bh49_recycled_5qubits.json -o "$work/bh49.receipt.json"
python3 -m ftqc_magic verify tools/ftqc_magic/examples/bh49_recycled_5qubits.json "$work/bh49.receipt.json"

# 3. Certify the non-triorthogonal [[9,1,2]] 9T -> 1T signed stabilizer factory over Z[omega_8]
python3 -m ftqc_magic check tools/ftqc_magic/examples/bgm_9_to_1_nontriorthogonal.json -o "$work/bgm9.receipt.json"

# 4. Run exact O(n^3) Clifford-hook fault enumeration on Z-flagged MSC-3 cultivation (d_T = 3)
python3 -m ftqc_magic cultivate --variant msc3_z_flagged --p 1/1000

# 5. Chain Z-flagged cultivation, Caliper escape and a recycled 15-to-1 factory
#    (the 15-to-1 stage adds a heuristic surface-code error per rotation; see the guide)
python3 -m ftqc_magic pipeline tools/ftqc_magic/examples/flagged_cultivation_caliper_15to1.json \
  -o "$work/cult_15to1.receipt.json"

# 6. Search the catalog for pipelines reaching 1e-10 conditional infidelity at p = 1/1000 (within the models above)
python3 -m ftqc_magic suggest --target T --p 1/1000 --infidelity 1/10000000000
```

## Included Engines & Verifiers

| Capability | Module | Primary References |
|---|---|---|
| Exact $\mathbb{Z}[\omega_8]$, $\mathbb{Q}(\sqrt{2})$, `PolyQ`, & outward `Decimal` bounds | `rings.py`, `ir.py` | exact cyclotomic and interval arithmetic |
| Binary / PPR / borrowed-identity factories, parity lift, puncturing, free-row extension | `matrix.py`, `catalog.py` | [Bravyi & Haah 2012](https://arxiv.org/abs/1209.2426) · [Campbell & Howard 2017](https://arxiv.org/abs/1606.01904) · [Jain, Wills & Singh 2026](https://arxiv.org/abs/2610.06535) |
| $\mathrm{GF}(2^s)$ self-dual trace binarization, Hermitian AG & 3-twirl concatenation | `galois.py` | [Gong, Pattison, Rall & Wills 2026](https://arxiv.org/abs/2610.06384) |
| Non-triorthogonal $\mathbb{Z}[\omega_8]$ stabilizer-collapse verifier & $n_{\min}(\tau, d)$ bounds | `limits.py` | [Balada Gaggioli & Mareček 2026](https://arxiv.org/abs/2610.31115) · [Beverland et al. 2020](https://arxiv.org/abs/1910.11827) |
| Mid-circuit check recycling ($15\to 1$ in 4 qubits, $49\to 1$ in 5 qubits) + LRAT lower bounds | `recycle.py`, `search.py` | [Jacinto et al. 2026](https://arxiv.org/abs/2606.07734) |
| Exact $O(n^3)$ Clifford-hook cultivation fault distance ($d_T$ vs $d_S$) & $Z$-flag restoration | `cultivate.py` | [Gidney, Shutty & Jones 2024](https://arxiv.org/abs/2409.17595) · [Chan, Strikis, Sun & Cai 2026](https://arxiv.org/abs/2609.17706) |
| Open-boundary mid-circuit escape postselection & Caliper worst-case soft-gap verifier | `caliper.py` | [Mehta et al. 2026](https://arxiv.org/abs/2610.06828) |
| Single-shot 3D tricycle $\to$ 4D quadcycle qLDPC `CCZ` factories & AOD shuttle budgets | `qldpc_factory.py` | [Menon et al. 2026](https://arxiv.org/abs/2610.10731) |
| Litinski asymmetric sizing, Gidney–Fowler `C2T` catalysis, & pipeline composer | `pipeline.py`, `bridges/` | [Gidney & Fowler 2019](https://arxiv.org/abs/1812.01238) · [Litinski 2019](https://arxiv.org/abs/1905.06903) |

See [GUIDE.md](GUIDE.md) for subcommands, the receipt schema, the correlation and catalysis hazards, the known gaps and the bibliography, and [../FTQC.md](../FTQC.md) for how `ftqc_magic` connects to `ftqc_codes`, `ftqc_faults`, `ftqc_decoders`, `ftqc_extractors`, `ftqc_compose`, and `ftqc_uncompute`.
