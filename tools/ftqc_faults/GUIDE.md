# ftqc_faults — GUIDE

Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: package with two engines and a bridge. `audit/` is the location-fault
engine ([guide](../ftqc_audit/GUIDE.md)); `zx/` is the ZX fault-equivalence
engine ([guide](../ftqc_zx/GUIDE.md)); `steane.py` checks the two engines
against each other on the five-ancilla Steane extractor; `action.py` gives the
single-fault spread of an installed Clifford route; `staged.py` gives the
one-location responses of a verified-cat Pauli check for the staged gates of
`ftqc_gadgets`. The packages `tools/ftqc_audit` and `tools/ftqc_zx` are
compatibility entry points that forward here. Python standard library; Stim
and CaDiCaL are optional. Battery class: selftest (three entries:
`ftqc_audit`, `ftqc_zx`, and `ftqc_faults` for the bridge and the
compatibility paths).

## PURPOSE

What faults do inside a scheduled Clifford circuit, and whether two circuits
behave the same under faults.

- **Audit engine.** Input: one circuit of preparations, Clifford gates, idles,
  measurements and one-bit Z feed-forward, with declared detectors, logical
  observables and outputs. Output: for every location and every Pauli, which
  detectors fire, which observables flip and whether the output is wrong;
  circuit and fault distance as intervals with replayable witnesses; strict
  fault tolerance through a low order or its first counterexample; a census of
  locations. Fault weight is the number of faulty locations.
- **ZX engine.** Input: two Clifford circuits or ZX diagrams with explicit
  boundary, branch and noise conventions. Output: w-fault equivalence of the
  two, after Rodatz, Poór and Kissinger, checked by two independent
  computations; constructive ZX rewrite proofs checked in time linear in their
  size; synthesis of cat-state, Shor and Steane gadgets; flag-conditioned
  decoder tables; and exact-rational stochastic contracts for one scheduled
  gadget.
- **Steane bridge.** Runs both engines on the 15-CNOT, five-ancilla Steane
  extractor and requires them to agree on every single-location fault.

## SCOPE

Clifford circuits on qubits. Pauli faults at the listed locations: any
nonidentity Pauli after a preparation, gate or idle (all 15 two-qubit Paulis
after a CNOT are one location), a flipped outcome at a measurement. Statements
hold for every fault set up to a stated weight, usually small; the ZX engine
reaches large gadgets through rewrite proofs instead of enumeration.

## USE-WHEN

- A proof counts the locations of a concrete gadget that can leave a heavy
  error, or claims that a verification step catches every single fault.
- A CNOT schedule changed and hook errors may now align with a logical
  operator.
- A compact or rescheduled gadget should be shown fault-equivalent to its
  reference, or a flag gadget needs its decoder table.
- A logical construction in `ftqc_extractors`, `ftqc_gadgets` or
  `ftqc_uncompute` exports a Clifford circuit for a physical check.

## NOT-FOR

- Code distance (`ftqc_codes`), decoders (`ftqc_decoders`), or gadgets on
  unknown inputs with decoders (`ftqc_gadgets`).
- Non-Clifford gates, general adaptive programs, leakage, coherent errors and
  faults correlated across locations.
- Probabilities beyond the single-fault level, except the ZX contract, which
  charges every multi-fault event as a failure.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"`:

    python3 -m ftqc_faults audit analyze surface:3,2,z,bad -o hook.json
    python3 -m ftqc_faults audit verify hook.json --full
    python3 -m ftqc_faults audit records probe.json --stim
    python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
    python3 -m ftqc_faults zx contract steane_5anc --p 1/10000 --idle-ratio 1/10 --noise-model circuit_css --code steane_7_1_3
    python3 -m ftqc_faults steane -o steane.json
    python3 -m ftqc_faults steane --verify steane.json

`audit` and `zx` take the subcommands documented in their guides. The old
forms `python3 -m ftqc_audit ...` and `python3 -m ftqc_zx ...` are equivalent.
Library: `ftqc_faults.audit`, `ftqc_faults.zx`, `ftqc_faults.action`.

Exit codes: 0 ok; 1 a finding (counterexample, disagreement, receipt
mismatch); 2 usage or malformed input.

## OUTPUTS

`ftqc_audit-receipt-1` (`faultcert-receipt-1` receipts are also read),
`ftqc_zx-contract-v1` and the Steane bridge receipt, which binds both probe
circuits and the ordered response ledger by digest.

## GATES

Measured on a four-core Linux virtual machine with Python 3.10, without Stim
or CaDiCaL: `ftqc_audit` 5,637 checks in 0.7 s (two named skips),
`ftqc_zx` 10 legs in 5 s, `ftqc_faults` 5 tests in 0.5 s. The last runs the
Steane bridge on both sectors, checks that a changed output Pauli or a changed
measurement record is detected, and checks that the `ftqc_audit` and
`ftqc_zx` names resolve to the same modules and entry points. The engine guides list their own
gates.

## FOOTGUNS

- Both engines count one faulty location as one fault, but they measure the
  damage differently. The audit engine's strict fault tolerance uses the joint
  Pauli weight of the output error; the ZX `circuit_css` model counts the X and
  Z sectors separately. One CNOT fault that leaves X on one qubit and Z on the
  other is an output error of weight two in the first and of weight one per
  sector in the second. Compare results only under the same convention.
- Detectors decide what is undetected. With discovered detectors every
  deterministic parity counts as a check.
- Idles exist only where written. The Steane bridge writes every live idle;
  library circuits write the idles their builders insert.
- Audit numbers measurements in tick order and ZX in circuit order; the
  bridge matches labeled records before comparing them.

## KNOWN GAPS

- The Steane bridge covers single faults. Agreement at higher weight is
  tested only on random circuits inside each engine's battery.
- The single-fault spread of an installed route (`action.py`) is a count of
  data weights, not a circuit distance; that would need a detection or
  recovery protocol.

## Five-ancilla Steane extractor

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
python3 -m ftqc_faults steane -o "$out/steane.json"
python3 -m ftqc_faults steane --verify "$out/steane.json"
```

For each of the X and Z syntheses, the bridge checks the ZX rewrite proof and
equivalence in both directions on the Steane codespace, with the independent
ZX frame oracle enabled. The `w=3` comparison uses the CSS-sector convention
and tests fault classes of weight below three.

The physical comparison takes the synthesized 15-CNOT circuit and its
ten-tick schedule, with every live idle explicit. The audit probes begin with
an ideal encoded zero and end with ideal measurements of all seven data
qubits, once in Z and once in X. They declare the flag and the three syndrome
parities as detectors, and keep all five raw record displacements before
recovery, relative to the same fault-free branch.

Each sector has 95 noisy locations: five preparations, 15 CNOTs, five
measurements and 70 idles. The model tests all 15 nonidentity Paulis after
each CNOT, all three Paulis after each preparation or idle, and one outcome
flip per measurement. For each of the 455 modes, the ZX scheduled propagator
must agree with the audit adjoint signatures and with forward replay on the
four detectors, all raw record flips and the X and Z masks of the seven
outputs. A changed result fails replay.

## Import names and receipts

The canonical imports are `ftqc_faults.audit` and `ftqc_faults.zx`; the
imports and command lines `ftqc_audit` and `ftqc_zx` are equivalent. Circuit
receipts and stochastic contracts have their own schemas and verifiers.
Published injection cases and decoder sampling belong to `ftqc_extractors`
(`python3 -m ftqc_extractors injection bb-uniform`); the audit engine's
injection and sampling commands forward there. [Injection guide](../ftqc_extractors/INJECTION.md)

## Code automorphisms

The installed-action verifier uses `ftqc_faults.action` to propagate each
signed data X and Z generator through physical Clifford ticks. Its
single-fault model inserts each nonidentity Pauli after a gate and X, Y or Z
on every idle site per tick, and counts the resulting data weights, including
faults that spread from spectators. See the
[installation workflow](../FTQC.md#code-automorphisms). The cat-check
responses in `staged.py` serve the [staged logical gates](../ftqc_gadgets/STAGED.md).

CREDIT: the engines credit their sources in their own guides: Chen et al.
(CAV 2025) and Aaronson–Gottesman for the audit engine, Rodatz, Poór and
Kissinger for the ZX engine, Stim (Gidney) for the optional cross-checks.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
