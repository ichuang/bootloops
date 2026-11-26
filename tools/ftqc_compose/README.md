# ftqc_compose

Build and replay finite bad-support calculations from He, Nguyen and Pattison, §2. A support is bad when it contains a listed generator. The receipt counts inclusion-minimal generators by size. This polynomial gives a union bound when the event family covers physical failures and the supplied noise law obeys the local-stochastic subset condition.

With `tools` on PYTHONPATH, run `python3 -m ftqc_compose check SPEC -o RECEIPT` and `python3 -m ftqc_compose verify SPEC RECEIPT`. The `family` input takes an `expression` whose operations are `family` (sites and minimal supports), `union`, `product` (left and right), or `substitute` (outer expression and an expression for each outer site under `at`). Products and substituted children require disjoint physical domains. Every generated family receipt checks every subset against an independent Boolean interpreter, up to 16 physical sites. Larger domains are refused. [Paper](https://arxiv.org/abs/2508.08246)
