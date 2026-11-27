"""Regenerate Liu et al. benchmark inputs from their pinned research checkout.

Run as a script with the optional pinned reference environment (see INJECTION.md).
"""

import argparse
import json
from pathlib import Path
import subprocess

REVISION = "09b005ef4f2038baaa5fea565468ec18044567e7"


def export(root, out, circuits):
    import qldpc_msi
    from qldpc_msi import tour_bb, xu_hgp
    from qldpc_msi.data_mgr import InitCfgMgrV2
    from qldpc_msi.noise_spec import get_noise_model

    if (
        Path(qldpc_msi.__file__).resolve().parents[2] != root.resolve()
        or subprocess.check_output(
            ["git", "-C", root, "rev-parse", "HEAD"], text=True
        ).strip()
        != REVISION
        or subprocess.run(
            [
                "git",
                "-C",
                root,
                "diff",
                "--quiet",
                "HEAD",
                "--",
                "src",
                "data/init_cfgs_v2",
            ]
        ).returncode
    ):
        raise ValueError("use an unmodified reference checkout at " + REVISION)

    def mask(row):
        return sum(1 << i for i, v in enumerate(row) if v)

    for n, facade, ids in [
        (144, tour_bb, [491]),
        (225, xu_hgp, [32, 45, 38, 35, 4, 5]),
    ]:
        code = facade.make_code(n=n)
        if list(code.all_dat_uids) != list(range(n)):
            raise ValueError("unexpected physical coordinate order")
        mem = dict(n=n, hx=[mask(r) for r in code.Hx], hz=[mask(r) for r in code.Hz])
        xs, zs = code.get_logical_ops_stim()
        for idx in ids:
            spec = InitCfgMgrV2().load(
                code_id=code.identifier,
                noise_model="li_uniform",
                init_method="max_supp_protection",
                remove_dummy_cnots=False,
                _id=idx,
            )
            cfg = spec["init_cfg"]
            pairs = [
                [mask(xs[i].to_numpy()[0]), mask(zs[i].to_numpy()[1]), 0, 0]
                for i in cfg.inj_lids
            ]
            prep = dict(
                carriers=cfg.inj_uids,
                plus=cfg.dat_uids_init_in_x,
                zero=cfg.dat_uids_init_in_z,
                bells=cfg.ghz_uids_lst,
            )
            # Stim's initial state is |0> on sites omitted by the reference builder.
            present = set(
                prep["carriers"] + prep["plus"] + prep["zero"] + sum(prep["bells"], [])
            )
            prep["zero"] = sorted(prep["zero"] + list(set(range(n)) - present))
            fixed = dict(X=sorted(cfg.fixed_x_stab_ids), Z=sorted(cfg.fixed_z_stab_ids))
            spec = dict(
                schema=1,
                kind="injection",
                method="selected",
                memory=mem,
                logicals=pairs,
                preparation=prep,
                postselect=fixed,
            )
            name = f"{'bb144' if n == 144 else 'hgp225'}-injection-{len(pairs)}"
            (out / (name + ".json")).write_text(json.dumps(spec, indent=2) + "\n")
            print(
                name,
                idx,
                cfg.inj_lids,
                {a: len(v) for a, v in fixed.items()},
                flush=True,
            )
            if circuits is None or idx not in (491, 5):
                continue
            cases = (
                [
                    ("uniform", 0.001, 2),
                    ("uniform-three", 0.001, 3),
                    ("asymmetric", 0.001, 3),
                ]
                if n == 144
                else [("uniform", 0.0001, 3)]
            )
            for model, p, rounds in cases:
                if n == 144:
                    noise = (
                        "li_uniform"
                        if model.startswith("uniform")
                        else "li_twoq_dominant"
                    )
                    actual = InitCfgMgrV2().load(
                        code_id=code.identifier,
                        noise_model=noise,
                        init_method="max_supp_protection",
                        remove_dummy_cnots=False,
                        _id=idx,
                    )["init_cfg"]
                    if actual != cfg:
                        raise ValueError("asymmetric source preparation differs")
                    ideal = facade.build_msi_circuit(
                        code, init_cfg=cfg, n_cycles=rounds
                    )
                    circ = get_noise_model(
                        noise, p
                    ).noisy_circuit_skipping_mpp_boundaries(ideal)
                else:
                    params = dict.fromkeys(
                        ("p_i", "p_state_p", "p_m", "p_CX", "p_idling_gate"), p
                    )
                    circ = facade.build_msi_circuit(
                        code, init_cfg=cfg, n_cycles=rounds, error_params=params
                    )
                (circuits / (name + "-" + model + ".stim")).write_text(
                    str(circ.flattened()) + "\n"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream", type=Path)
    parser.add_argument("examples", type=Path)
    parser.add_argument("--circuits", type=Path)
    args = parser.parse_args()
    export(args.upstream, args.examples, args.circuits)
