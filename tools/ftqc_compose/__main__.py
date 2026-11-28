"""Generate and replay content-bound deterministic composition receipts."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

from ftqc_codes.certify import Refusal as CodeRefusal
from .supports import family_result


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def receipt(spec):
    from .contracts import Violation

    try:
        result = evaluate(spec)
    except Violation as failure:
        result = {
            "status": "counterexample",
            "scope": "accepted-diagonal-pauli",
            "failed_obligation": "leaf simulation",
            "witness": failure.evidence,
        }
    return {"schema": "ftqc_compose-1", "input_sha256": digest(spec), "result": result}


def evaluate(spec):
    if spec["kind"] == "staged-program":
        from .staged import program
        return program(spec)
    if spec["kind"] == "rm15-action":
        from ftqc_magic.action_program import transport
        return transport(spec)
    if spec["kind"] == "installed-action":
        from ftqc_gadgets.installed_action import analyze
        return analyze(spec)
    if spec["kind"] == "decoded-surgery":
        from .surgery import decoded_channel
        return decoded_channel(spec)
    if spec["kind"] == "surgery-channel":
        from .surgery import channel
        return channel(spec)
    if spec["kind"] == "injection-batch":
        from .injection_batch import result

        return result(spec)
    if spec["kind"] == "local-stochastic":
        from .stochastic import stochastic_result

        return stochastic_result(spec)
    if spec["kind"] == "cone-bound":
        from .tails import cone_bound_result

        return cone_bound_result(spec)
    if spec["kind"] == "volume":
        from .volumes import volume_result

        return volume_result(spec)
    if spec["kind"] == "filter-contract":
        from .filters import filter_contract_result

        return filter_contract_result(spec)
    if spec["kind"] == "magic-interleaved":
        from .magic_correlations import interleaved

        return interleaved(spec)
    if spec["kind"] == "magic-batch":
        from .magic_batch import consume_batch

        return consume_batch(spec)
    if spec["kind"] == "addressable-magic-t":
        from .addressable_magic import consume_addressable_lp_t

        return consume_addressable_lp_t(spec)
    if spec["kind"] == "insitu-retries":
        from .insitu import retries
        return retries(spec)
    if spec["kind"] == "bb144-rm15-price":
        from .bicycle_price import price
        return price(spec)
    if spec["kind"] == "bb144-rm15":
        from .bicycle_install import install
        return install(spec)
    if spec["kind"] == "insitu-rm15":
        from .insitu import factory
        return factory(spec)
    if spec["kind"] == "magic-logical-t":
        from .magic_handoff import consume_t

        return consume_t(spec)
    if spec["kind"] == "distilled-logical-t":
        from .logical_t import logical_t_result

        return logical_t_result(spec)
    if spec["kind"] == "block-reading-residual":
        from ftqc_extractors.__main__ import verify

        source, proof = spec["extractor"], spec["distance_receipt"]
        if (
            source["kind"] != "gross_block_reading"
            or proof["result"]["status"] != "proved"
        ):
            raise ValueError("expected verified BB144 block reading")
        verify(source, proof)
        t = spec["faults"]
        if type(t) is not int or t < 0:
            raise ValueError("fault count must be a nonnegative integer")
        probes = proof["result"]["circuit"]["probes"]
        spread = sum(p["residual_component_weight_per_fault"] for p in probes.values())
        return dict(
            status="proved",
            axis=source["axis"],
            faults=t,
            residual_data_weight_upper=min(288, t * spread),
            assumption="ideal encoded inputs; apply split frame; open output",
            probe_hashes={b: p["circuit_sha256"] for b, p in probes.items()},
        )
    if spec["kind"] == "gross-residual":
        from .gross_residual import gross_result

        return gross_result(spec)
    if spec["kind"] == "level":
        from .level import level_result

        return level_result(spec)
    if spec["kind"] == "hgp-repeat":
        from .hgp import repeated_result

        return repeated_result(spec)
    if spec["kind"] == "family":
        result = family_result(spec)
    elif spec["kind"] == "contract":
        from .contracts import contract_result

        result = contract_result(spec)
    else:
        raise ValueError("unknown study kind")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check")
    check.add_argument("spec", type=Path)
    check.add_argument("-o", type=Path)
    verify = commands.add_parser("verify")
    verify.add_argument("spec", type=Path)
    verify.add_argument("receipt", type=Path)
    args = parser.parse_args()
    try:
        result = receipt(json.loads(args.spec.read_text()))
        if args.command == "verify":
            if canonical(result) != canonical(json.loads(args.receipt.read_text())):
                raise ValueError("receipt differs from recomputed evidence")
            print("verified receipt:", result["result"]["status"])
            return 0
        content = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.o:
            args.o.write_text(content)
        else:
            print(content, end="")
        return int(result["result"]["status"] == "counterexample")
    except (ValueError, CodeRefusal, TypeError, KeyError, AttributeError, OSError) as error:
        print(f"ftqc_compose: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
