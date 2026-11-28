# ftqc_compose — GUIDE

Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: package (Python 3.10+, standard library; imports the sibling FTQC
packages, whose receipts it re-verifies). Support families and contracts:
`supports.py`, `contracts.py`, `hgp.py`, `level.py`, `stochastic.py`.
Schedules: `tails.py` (cone bounds), `filters.py` and `filter_models.py`
(purification trees), `volumes.py`. Distillation and logical T: `distill.py`,
`logical_t.py`, `magic_*.py`, `addressable_magic.py`. Measurement residuals
and channels: `gross_residual.py`, `surgery.py`, `injection_batch.py`.
In-situ factory: `insitu.py`, `bicycle_*.py`, `pivot.py`. Staged gates:
`staged.py`. Optional comparisons: `compare_autdec.py`,
`compare_bicycle_compiler.py`. Tests `test_*.py`, `examples/`. Battery class:
selftest.

## PURPOSE

Join operations whose receipts come from the other tools and compute what
follows for the whole: whether each output meets the next input's
requirements, and what failure probability follows under a stated noise law.
Every study is one JSON input to `check`, which writes a receipt that `verify`
replays. The framework is the composable fault tolerance of He, Nguyen and
Pattison (bad-support families, typed contracts, level reduction,
local-stochastic noise), extended by studies for fault-cone tail bounds,
purification-tree contracts, schedule volumes, distillation into a logical T,
residuals of logical measurements, magic-state consumers and staged gates.

## SCOPE

Results are exact rational computations on finite objects. Their meaning
depends on premises that the receipt lists: noise laws (independence or the
local-stochastic subset condition, declared by the user and checked only for
an explicit finite law on at most 16 sites), quoted theorems (Theorem 4.11 and
Proposition 4.9 of He, Nguyen and Pattison, the cone theorem and the filter
rule of `ftqc_gadgets`), supplied contracts, clean encoded inputs and ideal
classical processing. A result that rests on a premise is labelled
`conditional`.

## USE-WHEN

- Two or more certified operations are used in sequence or in parallel and
  their interfaces (wires, codes, logical bases, decoders, frames) must match.
- A fault-count receipt needs a probability: union bounds over bad supports,
  tail bounds over the fault cones of a schedule, a factory bound carried into
  a logical T.
- A schedule needs a cost in expected physical locations.

## NOT-FOR

- Discharging the quoted theorems or the noise hypotheses: they remain
  premises.
- Tight logical error rates. Most bounds are union bounds; see FOOTGUNS.
- Simulation or sampling, except where a study says so.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"`:

    python3 -m ftqc_compose check STUDY.json -o RECEIPT.json
    python3 -m ftqc_compose verify STUDY.json RECEIPT.json

`check` exits 0 for a proof, 1 for a checked counterexample (the receipt is
still written) and 2 for invalid input, a failed join or a refused resource.
`verify` exits 0 when the receipt replays faithfully and 2 otherwise.

## INPUTS

The `kind` field selects the study. Each has a section below.

| Kind | What it computes |
|---|---|
| `family` | minimal bad supports of a union, product or substitution of families, checked against a Boolean interpreter on up to 16 sites |
| `contract`, `hgp-repeat` | typed sequential and parallel joins under the accepted-diagonal-Pauli scope; the one-edge HGP measurement pair |
| `level` | level reduction of a DAG of gadgets (Theorem 4.11), optionally with a classical output check or `distill15` |
| `local-stochastic` | the probability bound of a checked support cover at an exact rate |
| `cone-bound`, `filter-contract`, `volume` | tail bounds over fault cones, purification-tree contracts, expected physical locations of a schedule |
| `distilled-logical-t`, `magic-logical-t`, `addressable-magic-t` | a factory bound carried through an encoder or transfer into a logical T |
| `magic-batch`, `magic-interleaved` | joint error laws of correlated factory outputs at their consumers |
| `gross-residual`, `block-reading-residual` | residual data and frame errors after a gross-code measurement or a two-block reading |
| `decoded-surgery`, `surgery-channel` | decoded or exact output law of an installed surgery instrument |
| `injection-batch` | joint accepted error law of simultaneous Y sources at ideal S consumers (needs Stim) |
| `insitu-rm15`, `insitu-retries`, `bb144-rm15`, `bb144-rm15-price` | the in-situ BB144 factory: source assignment, retries, installation and pricing |
| `installed-action`, `rm15-action` | a signed code action on couplers; RM15 axes in a verified frame |
| `staged-program` | joined stages of the switched T or pieceable CCZ under a pair-noise premise |

Inputs are trusted JSON; do not run files from sources you do not trust.

## OUTPUTS

A receipt with the input digest, `status` (`proved`, `conditional`,
`computed`, `counterexample`, `no bound`), the computed values as exact
rationals (with outward-rounded decimals where useful), and the list of
premises and ideal boundaries the result depends on.

## GATES

`python3 -m unittest discover -s tools/ftqc_compose -t tools` (the
registered battery): 90 tests, 3 skipped, in 22 s on a four-core Linux virtual
machine without Stim, just above the toolkit's 20 s standard. The cone and
filter tests compare the bounds with exact enumeration, exact accepted-output
laws and sampled physical failures; 17 planted defects in the cone bounds and
42 in the filter calculus each fail them. Receipt tests refuse tampered
evidence, stale digests and failed joins.

## FOOTGUNS

- **Union bounds are loose.** The decoded surface-code ZZ example has no
  single-fault failure, yet at p = 10⁻⁵ its bound is 1.2 × 10⁻⁴ because every
  pair of its 1,542 locations is charged. Read a bound as a certified ceiling,
  not an estimate.
- **The noise premise is per subset.** Local-stochastic noise bounds the joint
  probability of every set of locations, which per-site rates alone do not
  give.
- **Conditional rates need acceptance.** A bound on "accepted and wrong" is not
  a bound on "wrong given accepted" until an acceptance lower bound is
  supplied.
- **Contracts are hypotheses in `cone-bound`.** Derive them with
  `filter-contract` or say where they come from.
- **Correlated outputs.** Outputs of one factory batch or one catalysed
  stage are not independent; the magic studies refuse to treat them so unless
  distinct instances and an independence premise are declared.
- **Locations once.** Shared padding or a reused ancilla couples neighbouring
  leaves; each physical location must have one identity in the enclosing
  circuit.

## KNOWN GAPS

- No end-to-end theorem: the quoted theorems and the cone and filter
  arguments are not proved here. The `closure` study removes the declared
  contracts and the independence premises from a two-level bound; the
  deterministic statements it lists are still quoted.
- The bundle filter for arbitrary damaged quantum inputs is open for the HGP
  leaves, so friendliness cannot be required of them.
- A general noncommuting adaptive surgery contract and automatic use of
  surgery channels in the filter calculus are open.
- Physical circuits for the ideal services of the staged gates (cat seeds,
  extended seeds, terminal recovery) are not certified.

## Bad-support families

Build and replay finite bad-support calculations from He, Nguyen and Pattison, §2. A support is bad when it contains a listed generator. The receipt counts inclusion-minimal generators by size. This polynomial gives a union bound when the event family covers physical failures and the supplied noise law obeys the local-stochastic subset condition.

With `tools` on PYTHONPATH, run `python3 -m ftqc_compose check SPEC -o RECEIPT` and `python3 -m ftqc_compose verify SPEC RECEIPT`. The `family` input takes an `expression` whose operations are `family` (sites and minimal supports), `union`, `product` (left and right), or `substitute` (outer expression and an expression for each outer site under `at`). Products and substituted children require disjoint physical domains. Every generated family receipt checks every subset against an independent Boolean interpreter, up to 16 physical sites. Larger domains are refused. [Paper](https://arxiv.org/abs/2508.08246)

## Typed contracts

A `contract` study supplies an ordered tree under `program` using `sequence` and `parallel`. The `accepted-diagonal-pauli` scope tracks correctness of accepted branches under stochastic Pauli faults. It permits correlated fault supports. Every join matches the ordered wire coordinates, classical or quantum kind, code, logical basis and decoder. Avoidance of the output's bad residual supports must imply avoidance of the next input's bad supports. Sequential and parallel composition take the union of the physical bad families on disjoint location domains.

The built-in classical leaves are three fallible preparations of zero (`prepare0`), an ideal classical `majority` decoder and an ideal `bit-wire`. A `hypothesis` leaf can name other simulation obligations and one factor of a friendly filter per input bundle, or null if friendliness is unproved. These declarations yield a conditional receipt with its unresolved premises. Set `require_closed` to require checked leaves, and `require_friendly` to require a filter derivation. Parallel composition needs both filters. Sequential composition inherits the first gadget's filter under a valid join (Prop. 4.9). A friendly conditional receipt depends on its declared hypotheses. Fully coherent faults and unconditional postselection probabilities require separate proofs.

## One-edge HGP contract

Use `{"kind":"hgp-repeat","variant":"edge-extractor"}` to check two nondestructive measurements of logical Z1 on the same eight-qubit memory. Use `edge-hook` to replay the negative control. The underlying `hgp-z` leaf takes an instance name and a variant. It can also appear in an ordered contract tree. Its output bundles are the updated quantum memory and one classical measurement bit. A second measurement in parallel with a classical wire preserves the first reported bit.

The promise assumes an encoded input with no physical error, the schedule's noisy gates and classical feedback operation, ideal syndrome tests at each output boundary, and ideal storage of the reported bits. The check conditions on acceptance. It covers diagonal Pauli alternatives at the listed locations. The bundle filter for arbitrary damaged quantum inputs is open, so requesting friendliness is refused. Basis masks (X1,Z1)=(3,5) and (X2,Z2)=(80,48) use qubit zero as the low bit. The signed Heisenberg flows check the ideal instrument on arbitrary encoded states. After projection onto a measurement branch, a stabilizer or the measured Z contributes only a phase. The output test reduces the data Pauli modulo these branch stabilizers and checks the reported bit and spectator. This certifies a clean conditional output at each successful join.

Each good leaf rejects 121 single-location modes and accepts 12 with the correct readout and residual. Taking every two-location subset within either leaf gives 342 sufficient minimal generators for the composed cover, with enumerator 342 x². Some pairs in that cover are harmless. A separate flat simulator rebases the second measurement's record controls and compares both boundary states with independent audit transport for every zero-, one- and two-location assignment. Its witness establishes distance two. In the hook circuit a Z fault on the target of the first extra CNOT (location 19, mode 8) escapes all checks and flips the spectator X test in both measurements.

A checked counterexample makes `check` exit 1 and still writes its receipt. A successful `verify` exits 0 for a faithfully replayed proof or counterexample. Invalid inputs, failed joins and receipt mismatches exit 2. A conditional receipt lists the hypotheses on which its conclusion depends.

## Define and replay a contract

This closed classical example recovers a zero bit after up to one faulty preparation. The three physical location IDs must be unique in the enclosing circuit.

```json
{
  "kind": "contract",
  "scope": "accepted-diagonal-pauli",
  "require_closed": true,
  "require_friendly": true,
  "program": {
    "op": "sequence",
    "left": {"op": "prepare0", "sites": ["p0", "p1", "p2"]},
    "right": {"op": "majority"}
  }
}
```

An external hypothesis has fields `op` (hypothesis), `id`, `ideal`, `inputs`, `outputs`, `faults` (a family leaf) and `filters`. Each ordered port supplies `kind` (Q or C), `wires` (ordered local coordinates), `code`, `basis`, `decoder`, and `bad` (minimal forbidden residual supports in these coordinates). For a nontrivial code those identifiers refer to the mathematical data of the premise. Checking that two identifiers match does not discharge the premise. A closed proof uses checked leaves only.

The canonical battery is registered as `selftest`. It requires Python 3.10+ and the sibling packages in this checkout. It covers empty and overlapping families, physical aliasing, all 512 two-level majority inputs, incompatible joins, conditional filters, tampered evidence, and the HGP proof and hook control. Replaying a receipt rebuilds the derivation against the input SHA-256. A scope or resource refusal produces no new receipt.

Credit: Zhiyang He, Quynh T. Nguyen and Christopher A. Pattison, *Composable Quantum Fault-Tolerance*, §§2 and 4. The support operations and typed composition rules follow their framework. The classical fixtures and finite diagonal-Pauli HGP checks are BootLoops implementations. The HGP ideal-flow checks come from ftqc_extractors and its forward transporter from ftqc_faults.audit.

## Level reduction

A `level` study constructs the physical circuit by expanding an ordered DAG of gadgets. Name its external input bundles, output wires, physical sites and vertices. Each vertex binds a unique ID to a gadget and ordered input/output wire names. Quantum and classical wires are consumed once. All vertex location domains, including ideal sites, must partition the stated physical universe. The `outer` family is over all contracted vertex IDs. Each of its generators is replaced with the disjoint product of the local bad families. The result records the minimal physical support cover.

Theorem 4.11 needs trivial external codes and output-adjacent gadgets with empty bad families. The implementation supports individual unencoded bit and qubit bundles at the boundaries. Native prepare0, majority, vote3 and bit-wire discharge unconditional local obligations. For another gadget, an `assume-4.8` node explicitly wraps a hypothesis for unconditional simulation with environment, non-diagonal Pauli terms and bundle-wise friendly filters. The result remains conditional. Accepted-run contracts cannot serve as this premise.

The theorem gives an avoiding logical fault, potentially a general operation with environment. To conclude correctness one must also establish that every such logical fault preserves the intended output. Set `output_check` to `classical` for exhaustive logical truth-table verification on the native classical network, or `none` for the simulation conclusion alone. Use `require_closed` to refuse declared premises. The classical checker limits external bits plus vertices to 16; the explicit family lift limits generated terms to 200,000.

## Local-stochastic bounds

A `local-stochastic` study embeds a checked `source` study and an exact rate `p`, for example `"0.01"` or `"1/100"`. The checker regenerates the support cover and evaluates its nonnegative weight enumerator using rational arithmetic. It reports the raw sum and the probability bound capped at one. The premise bounds the joint inclusion probability for **every** subset of locations. It allows correlations. Per-site marginals alone are insufficient.

Optionally provide a finite `law`, a list of objects with fields `support` and `mass`. The checker verifies normalization and every subset inequality on up to 16 physical sites, then calculates the actual covered-event probability. Without a law, the distribution condition is an explicit assumption. A family source gives an event bound. A level source with a checked logical output proof also bounds the unconditional output distance. An accepted-run contract bounds the joint probability of acceptance and harm under its input promise. A conditional error rate on accepted runs would also require an acceptance lower bound.

## Tail bounds over fault cones

A `cone-bound` study takes a block-level `program` in the format of `ftqc_gadgets cones`, an exact rate `p`, the `inner` code's `n` and `t`, and `contracts` for the ancilla kinds the program uses. It recomputes the cones, checks that linear propagation stays inside them, and sums a bound over the decoding events. A decoding event is wrong only if more than t of its n coordinates are dirty, and a coordinate is dirty when one of the cone's L locations there is faulty or an ancilla in the cone carries an error there.

When `outer` is given every block of the program is an outer block, also when the program delivers no outer resource, and its outputs are decoded at both levels. Set `noise` to `independent` or `local-stochastic`. For a fixed independent mode with location rate p and ancilla rates q₁,…,q_k, a coordinate is dirty with probability d = 1 − (1 − p)^L (1 − q₁)⋯(1 − q_k). For n independent coordinates at that rate, the binomial tail Pr[Binomial(n, d) > t] is the probability of exceeding the decoder radius. It upper-bounds decoder failure because decoding can still succeed above t. When the supplied rates are upper bounds, the same calculation bounds the radius-exceedance probability. A contract is then `{"pi": ..., "modes": [[w, q], ...]}` per sector: apart from a junk event of probability `pi`, charged where the ancilla is delivered and at each half of a correction, every increasing event of the ancilla's dirt is at most the sum over modes of w times its probability under independent dirt of rate q. Several ancillas in one cone multiply out, with identical ones grouped. Under the local-stochastic subset condition, more than t dirty coordinates needs one faulty location in each of t + 1 coordinates, which gives C(n, t+1)(Lp)^(t+1). A contract is then `{"equivalents": m}`, a number of location equivalents per coordinate added to L. Only single-level programs are bounded under this condition.

For a two-level program supply `outer` as well. The sum of the inner-event failure bounds in an outer cone bounds the probability of a dirty column. Under the declared independence of columns, the binomial tail for exceeding the outer radius then upper-bounds outer decoder failure. Junk probabilities of outer resources are top-level terms.

The result is `conditional`. Its premises are the cone theorem for the program, the declared contracts, the noise condition, and for two levels the independence of columns. The contracts are hypotheses: this study does not derive them. Rates and final values are exact rationals rounded upward to 50 significant digits. Mixtures are multiplied out in fixed point with 200 digits, clean probabilities rounded downward and weights upward, so every rounding raises the bound. When a cone holds many ancillas with several modes, the terms that matter most for the tail are kept and the rest are merged at their dirtiest rate; the result counts how often that happened.

The tests compare the binomial tail with enumeration, the threshold count with an explicit bad-support family of this package, mixtures with their exact expansion, every inner and outer term with the formula written out, and sampled failures of a physical circuit (through `ftqc_gadgets.cones.Direct`) with the bound. Seventeen planted defects each fail them; one needed a new check first.

Two levels can also be bounded under `local-stochastic` noise, without independence of columns. A set of faulty sources that makes a column dirty lies in that column, and different columns have disjoint sources, so the probability that t₂ + 1 given columns are all dirty is at most the product of their union weights W. The outer event then fails with probability at most C(n₂, t₂+1)·W^(t₂+1). Inner ancillas keep their contract in location equivalents. An outer ancilla needs `{"column_weight": w, "junk": j}` per sector: apart from an event of probability j, charged where it is delivered, the columns of any set A are all dirty with probability at most w^|A|, jointly with the physical locations. The premise list names the disjointness of column sources in place of the independence of columns. The tests compare every term with the formula written out and with the independent-noise result for the same model.

## Regions of several noise rates

A `regions` study takes a checked `source` (as for `local-stochastic`) and `kinds`, a map from every non-ideal physical site to a named kind of location. When a location of kind k is faulty with probability c_k·ε^(w_k), the union bound over the minimal bad supports has leading power ε^τ with τ(w) the minimum of w·n over the vectors n of fault counts by kind. The result lists these vectors with their multiplicities, the componentwise-minimal ones, and for two kinds the vertices of the lower hull with the exponent ratios at which the dominant vertex changes. Each entry of the optional `exponents` list gets its τ, the dominant vectors and the number of minimal supports behind them. With independent faults the bad-support probability has the same exponent and leading coefficient; the tests check this against exact enumeration.

```json
{"kind": "regions", "source": {"kind": "family", "expression": {"op": "family",
   "sites": ["d1", "d2", "m1"], "minimal": [["d1", "d2"], ["d1", "m1"]]}},
 "kinds": {"d1": "data", "d2": "data", "m1": "meas"},
 "exponents": [{"data": "1", "meas": "1/2"}]}
```

## Filter calculus for purification trees

A `filter-contract` study derives the contract that a `cone-bound` study takes as a hypothesis: what an accepted purification tree delivers, given that every check on the way accepted. A contract of one sector is a junk probability `pi` and modes `[[w, q], ...]`. Apart from the junk event, the block's error is equivalent to one supported on a random set of dirty coordinates, and every increasing event of that set is at most the sum over modes of w times its probability under independent dirt of rate q.

A leaf with J junk-capable locations and at most L local locations at one coordinate has pi = 1 − (1 − p)^J and one mode (1, 1 − (1 − p)^L). One check of sector s on a survivor A and a partner B follows the filter rule. The survivor is junk only if both blocks are, or if the dirty sets of the two sides together reach the minimum weight of the trivial patterns. Otherwise its old dirt survives only where the partner's side is dirty too, so rates multiply. When one side is junk the other side's dirt passes unmultiplied, with the junk probability as its weight. In the other sector junk probabilities and rates add. Every coefficient is divided by a lower bound on the probability of acceptance: (1 − ν_A)(1 − ν_B) times the probability that no step location that can cause rejection is faulty, with ν the probability that a block's pattern is not trivial. A reset step unites the fresh dirt of the check with the inherited part without letting the leading coefficient grow from check to check. Both statements are printed in the receipt (`check`, `reset`).

Set `level` to `physical` or `encoded`.

```json
{"kind": "filter-contract", "level": "physical", "p": "1/10000",
 "tree": {"n": 23, "wmin": {"X": 8, "Z": 7}, "sequence": "XZXZXZ"},
 "leaf": {"locations": 151, "X": {"junk": 70, "local": 5}, "Z": {"junk": 71, "local": 5}}}
```

A physical tree is given by counts, as above, or by `factory`, a receipt of `ftqc_gadgets --tree`. The receipt is verified again, and the classes of the leaf locations and the gate, measurement and idle locations of every check are read from the recognised tree. Rates of disjoint sets of locations combine as 1 − (1 − a)(1 − b). The [banked example](examples/filter-contract-golay.json) is the 64-leaf Golay tree at p = 10⁻⁴.

An encoded tree has code blocks for qubits, and its "locations" are decoding events of the level below. Supply `inner` (n and t), the `contracts` of the inner ancillas, and `raw`: the block-level program of the raw preparation with its corrections, and `last_layers`, the number of final layers that form its last round. The tree is given by `tree` or by a `factory` receipt with strict classes. The study bounds every decoding event of the raw preparation and of one check of each sector by the tail bound over its cone. Failures before the last round count as junk. Failures of the last round are local to a block, which is checked: their cones must share no source across blocks, and every block must end as a block just corrected.

Decoding events of one column can share sources, so two things change. Rates are only added, never combined as independent. And the dirt a block received in its latest step is carried through the next check of its sector without being multiplied by the partner's rate, because it can depend on the same inner ancillas as the partner's measured syndrome. The study does not take this on trust. It builds one column of the whole tree, computes the cones, and reports for every level whether the survivor's older dirt is disjoint from the partner's side (required, or the study refuses) and whether its latest dirt shares sources with it (then it is carried).

The result is `conditional`, or `no bound` when the acceptance bound of some check is not positive. It lists the contract in the form `cone-bound` reads, the acceptance bound of every level, the expected cost counting repetitions, and the premises: the filter rule for the tree, the noise condition, and for an encoded tree the cone theorem, the inner contracts and that dirt indicators are increasing functions of independent sources. The calculus itself is quoted, not proved here. `modes` (default 4) and `order` set how the final mixture is summarised: merging is upward and loses least for events that need `order` dirty coordinates. Bounds are rounded upward and acceptance downward at 50 significant digits.

The tests compare the calculus with two exact computations. One follows the filter rule on two coordinates and sums over everything, for three acceptance rules, including one that accepts exactly the harmful outcomes; with shared sources it shows that multiplying the latest dirt gives bounds that are too small by a factor of two, and that carrying it is sound. The other is the exact law of the accepted output of Steane trees under independent faults (`ftqc_gadgets.factories.accepted_law`), for uniform and one-sided fault labels. Forty-two planted defects each fail the tests; six needed a new check first.

### Local-stochastic noise

Add `"noise": "local-stochastic"` to a `filter-contract` study of either level. The premise on the noise is then only that the probability of all locations of a set S being faulty is at most p^|S|; faults need not be independent. Every number is the total weight of a family of fault sets that covers an event of the filter rule. Weights of covers on disjoint sources multiply, weights of alternatives add (rates unite as a + b), and the tail of "at least w dirty coordinates" is C(n, w) times a moment. Nothing is conditioned on acceptance: a rejected check is repeated with fresh inputs, the delivered block is the first accepted attempt, and its weights are those of one attempt divided by one minus the weight of rejection. The contract is `{"junk": j, "modes": [[w, q], ...]}` per sector: junk is covered with weight j, and "not junk, and every coordinate of a set A dirty" with weight Σ w·q^|A|. `cone-bound` reads this form under `local-stochastic` noise for inner and outer ancillas; `{"equivalents": k}` and `{"column_weight", "junk"}` are its one-mode special cases. For an encoded tree the inner `contracts` must be in this form. Module `filters_ls.py` states the calculus and the argument for each line.

The tests compare junk, every moment and the acceptance bound with exact laws of the small models, for independent sources, for shared sources, and for two leaf laws that satisfy the noise condition without independence. One of those laws violates the acceptance bound of the calculus for independent faults. On the Golay example the junk weights are within three percent of the independent-noise values.

## Closure of a two-level program

A `closure` study chains the local-stochastic studies so that no contract is declared: the contracts of the inner ancillas from their physical trees, the contracts of the outer ancillas from their trees of encoded blocks, and the tail bound over the cones of the two-level `program`. Give `p`, `inner` and `outer` (n and t), the `program`, and under `inner_ancillas` and `outer_ancillas` one tree per ancilla kind in the fields of a `filter-contract` study (`tree` and `leaf`, or `factory`; for an outer ancilla `tree` or `factory`, and `raw`).

The result is `conditional`. `upper` bounds the probability that a delivered ancilla is junk or that some decoding event of the program has more than its radius of dirty coordinates. `noise_premise` is the single probabilistic premise. `deterministic_premises` lists the quoted statements that the other studies check on the instance: the cone theorem, the filter rule, the acceptance rule with repetition, and that distinct ancillas, attempts and columns occupy distinct locations. `no_longer_premises` lists what the chain replaces: declared contracts, independence, conditioning. The [banked example](examples/closure-two-level-steane.json) has two Steane levels and gives 5.121e-13 at p = 10⁻⁶, scaling as p⁴.

## Volume of a schedule

A `volume` study prices what a schedule spends, in expected physical fault locations. It takes a `program` in the format of `ftqc_gadgets cones`, the `inner` code's `n`, for a two-level program the `outer` code's `n`, and `costs`: for each ancilla kind the expected number of physical locations spent to deliver one accepted ancilla, repetitions after rejection included. The counts come from `ftqc_gadgets.cones.ledger`.

The total has three parts. Online: n₁·n₂ locations for every live block in every layer. Inner corrections: each consumes one |+⟩ and one |0⟩ ancilla in each of the n₂ columns. Resources: every delivered resource of a costed kind, the two ancillas of every outer correction included. Kinds listed in `uncosted` are counted and left out of the total, and the result says so; a kind that is neither priced nor listed is refused. Blocks of kind `std` are existing blocks and cost nothing. The locations an ancilla occupies while a correction uses it are not counted: its cost covers its preparation only.

Costs are inputs. A physical `filter-contract` study reports the cost of its tree as `expected_cost`. An encoded `filter-contract` study given `inner_costs` (`{"plus": ..., "zero": ...}`) prices its raw preparation and one check of each sector with the same ledger, applies the expected-cost recursion (2·cost + check)/acceptance at every level, and reports the cost of the outer ancilla in physical locations, with the two priced pieces under `priced`. The [banked example](examples/volume-two-level.json) is a small two-level schedule with one outer correction.

`ftqc_gadgets.cones.correction_pattern` measures how a schedule follows its gates and idle ticks by inner corrections: how many such steps are followed directly by a correction, by a measurement, or by another step, and the lengths of the runs of steps without a correction. A schedule that corrects after every gate and every idle tick has runs of length one only.

## Correcting 15-to-1 distillation

Use `output_check: "distill15"` on a level-reduction circuit of 15 distinct `magic-source` preparations and one `rm15-ideal` processing vertex. The latter denotes a fault-free encoded-plus preparation, 15 coordinate-wise T injections, syndrome extraction with minimum-weight correction, and logical decoding. The output is unconditional. The sources are the fallible physical locations. An explicit support law must assign zero mass to faults on ideal sites.

The checker constructs the punctured Reed–Muller CSS code on columns 1 through 15, with the four coordinate functions as X checks. The odd coset encodes logical zero and the even coset logical one. In this basis transversal T gives phases ω⁷ and 1, which differ by the required logical T phase (ω = exp(iπ/4)). It enumerates every weight-one syndrome, checks syndrome correction for left/right Pauli pairs and checks the corrected logical Kraus amplitudes for injection on every coordinate with integer arithmetic in the eighth roots of unity. The equality of the two logical amplitudes proves that a single faulty resource changes only the scalar and syndrome on the corrected branches. Linearity covers a general one-resource error.

The logical cover must include every pair of sources and a fault at the distillation vertex, or a stronger cover. Substitution removes the ideal vertex's singleton. The supplied study gives 105 physical pairs and the bound min(1, 105p²) on output trace distance. It also checks that two distinct resource Z faults fool the decoder into a third Z correction and leave logical Z. In the banked correlated example at p = 1/100 each pair has mass 1/10000, with the remaining mass on no faults. These disjoint atoms attain the bound if the listed sources suffer Z errors, since Z|T⟩ is orthogonal to |T⟩. Independent noise generally gives a different result.

[He, Nguyen and Pattison, Theorem 4.11, Proposition 5.2 and Proposition 7.9](https://arxiv.org/abs/2508.08246)

## Gross-code residual data and frame

A `gross-residual` input embeds a native `gross_circuit` specification under `extractor` and its complete `distance_receipt`. Verification replays the distance computation and binds both complementary probes by their circuit hashes. It then transports every location and every allowed Pauli or measurement flip to the 144 output data qubits and the 18 X split records. The counterfactual ideal data readouts reveal the raw X and Z displacements. Their outcomes can be individually random. The fault-free checks and logical tests must still be deterministic.

For split displacement s, port q gets a Z update equal to the parity of s on its spanning-tree path. On each fault the checker evaluates all terminal checks and logicals on (raw X, raw Z XOR frame Z) and compares these signs and the split cycles to the audit model. It checks that the kernel of all terminal tests is generated by the X stabilizers in one sector and the Z stabilizers plus the measured Z in the other. Multiplying by the measured Z only contributes a phase on a projected branch. All results use clean encoded input and the native stochastic Pauli location model.

History and split-cycle acceptance leave an open data block. The receipt counts its accepted single-location modes with a nontrivial residual and gives forward-replayed examples. The reported representative weights are constructive upper bounds. A zero representative is exact by row reduction. For t locations, XOR of representatives gives the stated linear spread bound even when multiple detector flips cancel. Under the additional ideal terminal syndrome test the native distance lower bound rules out an incorrect bit or a nontrivial residual below that weight. A downstream circuit that only observes the history must carry the open residual.

## From a distilled qubit to logical T

A `distilled-logical-t` study embeds a `local-stochastic` factory with the `distill15` logical output check, a native extractor `t_gate` and its complete `gate_receipt`. Set `handoff` to `fresh-product-css`. The logical input may be entangled with a reference. The factory supplies a fresh resource independent of that joint input and subsequent faults. Factory processing uses the fault-free operations in the RM15 contract. The gate uses the native flags, repeated parities, and ideal terminal syndrome checks.

The encoder takes the consumer's actual HGP checks and positive logical X and Z representatives. It prepares the uniform positive superposition over the X stabilizers, copies the raw basis bit b into the logical X string, and erases the raw wire with the logical Z parity. Binary propagation through its CNOT list verifies both cosets and disentanglement for arbitrary amplitudes. The receipt binds the ordered schedule by hash. The other resource logicals start in a fixed state. Data-block spectators may be entangled with the addressed data qubit.

For parity p and X result m the ideal consumer measures Z⊗Z and resource X, then uses the frame Sᵖ Zᵐ. The checker constructs all four Kraus operators on the joint data and resource input, verifies their completeness with exact eighth-root phases, and matches their action to the physical consumer's verified branches. With |T⟩ the action is T. With Z|T⟩ it is ZT. Encoding by an isometry preserves trace distance. Appending a fresh resource and applying this trace-preserving instrument contracts it on data plus any reference, including when the transcript is retained. Consequently half the diamond distance on logical data with fault-free encoding and consumption is at most the factory bound δ.

The two exact rates `encoder_p` and `online_p` specify separate noise assumptions. For N encoder preparations and CNOTs at per-location rate e, the probability of any encoder fault is at most Ne. A Z on the raw control after its first fanout CNOT propagates to the encoded logical Z. This single fault passes the otherwise ideal consumer, so the deterministic encoder has accepted fault distance one. The receipt gives its physical location. The encoder model includes the listed operations, with ideal discard and classical frame processing.

For a native consumer case with L locations and verified lower distance d, a joint wrong-and-accepted test requires at least d online faults with the ideal resource. Under a local-stochastic online rate q, the union bound for this event is C(L,d)qᵈ. Adding δ and Ne bounds the same test with the supplied factory and encoder. All probabilities are capped at one. A test here is a verifier that would pass a correct output with certainty. For the complete adaptive signed-X readout the checker adds the two possible arms' fault tails, counting the resource contribution once. Noise bounds apply to potential locations before selection of an arm. A rate conditioned on acceptance would additionally need an acceptance lower bound.

## Two-block reading residuals

A `block-reading-residual` input embeds the `gross_block_reading` specification as `extractor`, its replayable `distance_receipt`, and a nonnegative integer `faults`. The result sums the two component spread bounds after the split frame and caps the resulting physical support by the 288 output data qubits. It applies to the open output under the declared clean-input and gate-location model. For one fault the X-reading fixture gives an upper bound of 74 and the Z-reading fixture gives 90. Both count residual support on the 288 data outputs after the frame update.

## Simultaneous Y source and ideal consumers

The [BB144 consumer example](examples/bb-injection-consumers.json) takes logical slots 2 and 3 of the four-output Y injection. In each ideal S gadget the data controls a CNOT into the resource, its Z measurement selects a Z correction on the data, and a wrong Y sign induces a Z after S. Replay checks both measurement branches on both input basis columns, which determines the linear map on an arbitrary data state. The interface assumes ideal terminal decoding and transfer. Its frame gives one detector-parity bit mask per resource in circuit detector order.

The `injection-batch` verifier rebuilds the actual location responses. Within one physical channel its Pauli alternatives are exclusive. Distinct listed channels are independent and unlisted channels are ideal. A finite convolution tracks the complete postselection syndrome and all output signs, then conditions on their shared acceptance. Set locations to the literal string `all` to request every noisy channel. Exceeding the state budget leaves that guarantee unavailable. The example enables three channels including the noisy CX on physical sites 38 and 254. Conditional on acceptance the two consumers both acquire a Z with probability 1/14988. The receipt supplies the full joint distribution. The empirical BP+OSD histograms in the sampling study are a separate evidence type.

## Signed surgery instruments

See [the surgery interface and replay examples](../ftqc_gadgets/SURGERY.md) for physical measurement receipts, installed feed-forward schedules and conditional residual laws.

### Magic-backed logical T

Run `python3 -m ftqc_compose check study.json -o proof.json` and replay with
`verify study.json proof.json`. A `magic-logical-t` study embeds the Magic
factory spec and strict receipt, an attempt ID, a native extractor gate and
receipt, and a Gadgets fresh-transfer receipt. Build the latter with
`magic_handoff.transfer_for(gate, port)`. The port names its block, ordered
logical X/Z pair, applied input frame and matching correction (I/I or Z/Z).
The source reader re-evaluates the target phases and exact joint accepted
Z-mask law. The consumer uses its conditional error once and reports the
joint accepted failure per original factory attempt. An optional
`encoded_origin` can forward the same verified encoded output with its ordered
data wires. The original noisy transfer is charged once. Separate factory,
transfer and online locations and a resource independent of the logical data,
reference and later faults are required. A final conditional logical error
requires a lower bound on the final acceptance probability.

### Correlated logical consumers

A `magic-batch` binds one verified factory attempt to ordered resource-wire
groups and distinct logical data ports. For a binary producer with an exact
mask law, each atom lists all affected data ports and its original accepted
mass and conditional mass. The T and CCZ injection checks enumerate every
measurement transcript and residual Z mask. A CCZ correction uses only
diagonal Clifford gates. A qLDPC pair uses six wires in two groups and
carries the native joint state bound. Its detailed Pauli law is unknown.

For an exact CCZ source, `twirl: joint_3twirl` samples an independent uniform
GL(3,2) change of basis per CCZ with diagonal Clifford compensation. Replay
transports the full mask distribution through the consumers. A generic
native cascade with only a normalized kernel has no verified downstream
mask transfer. For separate factory lanes, use `magic-interleaved` with
distinct instance IDs, a verified member from each, and an explicit
independent-noise premise. A terminal qLDPC/C2T experiment requires a
fresh catalyst with empty history and retirement at this batch. Lifetime
declarations apply to this workload. Physical non-Clifford preparation and
factory-to-memory distance require a separate circuit audit.

The [in-situ BB144 factory guide](../ftqc_extractors/BICYCLE_FACTORY.md)
covers injection batch replay, RM15 input assignment, signed native pivot
measurements, physical scheduling and full-attempt pricing.

## Code automorphisms

The `installed-action` study binds physical couplers, ticks and fault propagation to a verified signed code action. Use the ordinary `check` and `verify` commands. The `rm15-action` study transports the 15 rotation axes into a verified frame. The `bb144-rm15` installation checks the signed Gross orbit and replays physical return routes. See [code automorphisms](../FTQC.md#code-automorphisms) for examples and an optional Qiskit comparison.

## Staged programs

A `staged-program` study joins the stages of the switched T or the pieceable
CCZ, checks that consecutive ports match in block order, code digest, logical
representatives, signs and frame, and bounds the failure by the sum over
stages of the location pairs times p², under the premise that every pair of
distinct locations is jointly faulty with probability at most p². The banked
`examples/switched-t.json` gives 44431/10¹² at p = 10⁻⁶. See the
[staged logical gate guide](../ftqc_gadgets/STAGED.md).

## Discharged level reduction and multivariate renormalization

In a `level` study, a vertex with `"op": "verified-gadget"` takes a
`ftqc_gadgets-receipt-1` carrying a passing `coherent` certificate, replays its
hash and evidence, and emits zero undischarged premises so `require_closed: true`
returns `status: "proved"`. Reach for it when replacing an `assume-4.8`
hypothesis by a checked zero-postselection gadget; use `contract` instead when
working on postselected diagonal-Pauli branches.

A `renormalize` study certifies all-level concatenation contraction and workload
compilation from exact extended-rectangle malignant-pair tables under
`"noise": "local-stochastic"` or `"noise": "depolarizing"`, coupled to recursive
15-to-1 Reed–Muller distillation (`"scheme": "rm15-distill"`) or staged
non-Clifford gates (`"scheme": "staged-gadget"`), with optional `ftqc_uncompute`
HMR receipt binding under `workload.uncompute`. Reach for it when proving an
explicit positive rational threshold lower bound (`threshold_lower`) and
selecting a concatenation depth `chosen_level` with certified output variation
`tv_upper <= target_eps`; use `closure` instead for fixed two-level transversal
Clifford schedules with purification trees. When a staged non-Clifford component
carries a positive `linear_floor_coeff`, `threshold_lower` is `"0"` and
requesting `target_eps` below the floor returns `status: "counterexample"`.
Believe it when `verify SPEC.json RECEIPT.json` reproduces every rational field
and a positive linear-floor control is rejected.

CREDIT: Zhiyang He, Quynh T. Nguyen and Christopher A. Pattison, *Composable
Quantum Fault-Tolerance* (arXiv:2508.08246), §§2 and 4, Theorem 4.11,
Propositions 4.9, 5.2 and 7.9. The cone and filter studies use the cone and
factory certificates of `ftqc_gadgets`; the distillation check follows Bravyi
and Kitaev's 15-to-1 protocol.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
