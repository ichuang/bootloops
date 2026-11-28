# Direct logical gates with staged gadgets

The first study raises an unknown Steane qubit into the punctured Reed–Muller
code, applies transversal physical T-dagger and lowers it again. The second
implements logical CCZ directly on three Steane blocks. Both preserve
coherence with an external reference. Each check uses signed operators or
word amplitudes for every basis column, with a common normalization and phase
across columns. The examples run using only the Python standard library.

From the BootLoops root:

```sh
export PYTHONPATH=tools
out=$(mktemp -d)
python3 -m ftqc_gadgets staged tools/ftqc_gadgets/examples/staged-pieceable-ccz.json -o "$out/ccz.json"
python3 -m ftqc_gadgets verify "$out/ccz.json"
study=tools/ftqc_compose/examples/switched-t.json
python3 -m ftqc_compose check "$study" -o "$out/t.json"
python3 -m ftqc_compose verify "$study" "$out/t.json"
python3 -m ftqc_magic switched-source "$study"
python3 -m ftqc_uncompute price-staged "$study" --feedback-ticks 1 --terminal-wire-ticks 0 --seed-wire-ticks 0
```

## Codes and instruments

`ftqc_codes certify rm15` uses sites 0 through 6 in the existing Steane
order. Bit 0 of each mask denotes site 0. The extra eight sites carry a fixed
CSS state with X generators given by the three upper copies of the Steane
rows and the all-ones upper mask. A split computational-basis column is a
uniform coset of the seven split X generators. An RM column is a uniform
coset of the four RM X generators. Logical zero is the even column in both.
This choice accounts for the dagger in the physical transversal gate.

The raising instrument measures three RM Z checks dual to the three low
Steane X generators. Lowering measures those X generators. The supplied
eight-entry correction array is indexed by the raw gauge record, bit zero
first. Replay checks the projected words or signed X-projector expansion for
both logical columns in each of the eight equiprobable branches. The
lowering proof verifies factorization of the discarded eight sites.
[Anderson, Duclos-Cianci and Poulin](https://arxiv.org/abs/1403.2734)

In the CCZ study each triple addresses one site per block, in A,B,C order.
The four pieces have 9,6,6,6 gates. Each cut reads the three constant Z checks
on each of the three blocks until two consecutive full readings agree. The
X recovery uses only that reading. The terminal syndrome lookup retains the
three agreed readings to resolve the deferred phase errors. This is a finite
decoder synthesized from the stated fault set. A transcript unseen during
synthesis uses identity recovery. An adaptive cut that reaches four rounds
without agreement flags failure. The receipt checks that all admitted paths
settle. Removing a CCZ changes the relative phase of a basis column.
[Yoder, Takagi and Chuang, §IV.2](https://arxiv.org/html/1603.03948)

## Noisy locations and ideal services

A failure is an arbitrary nonidentity Pauli on a listed native gate after
its ideal action, or one listed preparation-output, idle or measurement
error. CCZ has all 63 correlated Pauli alternatives on its three legs.
Serial CCZ layers include an idle on every other data site. The transversal
phase is a parallel layer. At the input, each one-site Pauli is tested in a
separate fault-free run. The gauge and pieceable studies use cat-controlled
Pauli checks. Within a shot every cat leg couples to just one data site.
Each cat-output leg, two-qubit interaction, X readout, and uncovered data
idle is noisy. The switch repeats each gauge bit three times with majority
readout. The pieceable cut uses consecutive full-round agreement.

The cat seed, the entangled eight-site seed before its local output noise,
classical control, conditional correction, and terminal full syndrome and
recovery are ideal services where listed. A check is one modeled time step.
Physical circuits preparing those entangled seeds and extracting the terminal
syndrome would require additional fault certification and time. The reported
single-fault claim is conditional on these boundaries.

A propagated X on a CCZ leg induces a CZ on the opposite pair. Faults
retains the full quadratic phase until the terminal projection. The signed
Pauli expansion is combined modulo Z stabilizers before testing for
incompatible logical residuals at a common decoder record. Decoders exports
the table indexed by history, X-error syndrome and Z-error syndrome. To
test an installed decoder, supply its rows as `terminal_recovery` in the
gadget input. Each row is `[history,sx,sz,x_correction,z_correction]`.
Missing entries apply identity. Replay checks each admitted residual against
the supplied correction modulo the two endpoint stabilizer groups. Gadgets
returns a counterexample with the competing fault origins when these data
cannot distinguish two residuals. Receipt replay recalculates the complete
table, the ideal action and the location census.

## Composition and consumers

Compose requires consecutive ports to match in block order, code hash,
logical representatives, signs and frame. Instance IDs distinguish physical
locations. The claimed program action is I, T, CCZ or CCX. A CCX study wraps the CCZ
stage in two `target-h` stages on the third Steane block and charges both
transversal layers including idle controls. With clean encoded
input and all ideal services supplied, each stage restores its output through
one fault. If each pair of distinct reserved sites is jointly faulty with
probability at most p², a union bound gives the sum over stages of
binomial(max_locations, 2) p², capped at one. Adaptive slots reserve four
rounds per cut. This gives an upper bound under the stated noise model.

Extractors' `staged_install` study takes a Gadgets receipt, an injective
list of hardware sites (data first, reused cat meters second), undirected
cat-data `couplers` and unordered native `ccz_triples`. Replay checks the
couplers needed by each measured support and every scheduled triple.
Single-site native phases are assumed available.

Magic applies the closed switched T to a fresh ideal encoded plus. The output
has unit yield and a conditional infidelity bound under the Compose model.
Magic also accepts a `switched-source` spec with a nested `program` via its
ordinary `check` and `verify` commands. Uncompute prices the reserved online wire-ticks and charges the supplied
feedback delay and per-stage seed and terminal service budgets. Zero service
budgets show the online contribution. An end-to-end physical volume requires
budgets for the actual supplied services.

Each study lives in the existing package that owns its claim. The generic
Clifford audit, ZX verification and accepted diagonal-Pauli contracts retain
their respective fault models.
