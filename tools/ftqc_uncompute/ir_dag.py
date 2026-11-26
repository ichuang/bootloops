"""Qualtran-isomorphic typed arithmetic & uncomputation DAG (`ArithDAG`).

Captures high-level compute/uncompute chirality, spooky-pebble ghost states, and
superposition-masked accumulator recurrences, with deterministic lowering to
scoped `.kmx` (`KickMixCircuit`) and structural lifting from `.kmx`.

100% Python standard library, MIT-licensed, zero float literals.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import json
from typing import Any, Sequence

from .gf2 import GF2Field
from .ir_kickmix import CircuitBuilder, KickMixCircuit, TOFFOLI_LIKE_OPS


VALID_REG_KINDS = frozenset({"quint", "qgf2", "qmodint", "mask", "classical"})
VALID_ROLES = frozenset({"input", "output", "io", "scratch", "mask", "classical"})
VALID_CHIRALITIES = frozenset({"compute", "uncompute", "inplace", "phase", "stream_probe"})
VALID_UNCOMPUTE_MODES = frozenset({"none", "reverse", "hmr", "ghost", "masked"})


@dataclass(frozen=True)
class RegDecl:
    name: str
    kind: str  # "quint", "qgf2", "qmodint", "mask", "classical"
    width: int
    role: str = "io"
    modulus: int | None = None
    gf_modulus: int | None = None
    mask_kind: str | None = None  # e.g. "mult_mod_n", "uniform_mod_n", "phase_grad", "carry_runway"

    def __post_init__(self) -> None:
        if self.kind not in VALID_REG_KINDS:
            raise ValueError(f"Invalid RegDecl kind {self.kind!r}")
        if self.role not in VALID_ROLES:
            raise ValueError(f"Invalid RegDecl role {self.role!r}")
        if self.width < 1:
            raise ValueError(f"Register {self.name!r} width must be >= 1")

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "name": self.name,
            "kind": self.kind,
            "width": self.width,
            "role": self.role,
        }
        if self.modulus is not None:
            out["modulus"] = self.modulus
        if self.gf_modulus is not None:
            out["gf_modulus"] = self.gf_modulus
        if self.mask_kind is not None:
            out["mask_kind"] = self.mask_kind
        return out

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "RegDecl":
        return cls(
            name=str(d["name"]),
            kind=str(d["kind"]),
            width=int(d["width"]),
            role=str(d.get("role", "io")),
            modulus=int(d["modulus"]) if d.get("modulus") is not None else None,
            gf_modulus=int(d["gf_modulus"]) if d.get("gf_modulus") is not None else None,
            mask_kind=str(d["mask_kind"]) if d.get("mask_kind") is not None else None,
        )


@dataclass(frozen=True)
class DAGNode:
    node_id: str
    op_kind: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    chirality: str = "compute"
    uncompute_mode: str = "none"
    pair_with: str | None = None
    params: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.chirality not in VALID_CHIRALITIES:
            raise ValueError(f"Invalid chirality {self.chirality!r}")
        if self.uncompute_mode not in VALID_UNCOMPUTE_MODES:
            raise ValueError(f"Invalid uncompute_mode {self.uncompute_mode!r}")

    def to_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "op_kind": self.op_kind,
            "inputs": list(self.inputs),
            "outputs": list(self.outputs),
            "chirality": self.chirality,
            "uncompute_mode": self.uncompute_mode,
            "pair_with": self.pair_with,
            "params": self.params,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "DAGNode":
        return cls(
            node_id=str(d["node_id"]),
            op_kind=str(d["op_kind"]),
            inputs=tuple(str(x) for x in d.get("inputs", ())),
            outputs=tuple(str(x) for x in d.get("outputs", ())),
            chirality=str(d.get("chirality", "compute")),
            uncompute_mode=str(d.get("uncompute_mode", "none")),
            pair_with=str(d["pair_with"]) if d.get("pair_with") is not None else None,
            params=dict(d.get("params", {})),
        )


@dataclass(frozen=True)
class ArithDAG:
    name: str
    registers: tuple[RegDecl, ...]
    nodes: tuple[DAGNode, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def reg_map(self) -> dict[str, RegDecl]:
        return {r.name: r for r in self.registers}

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "ftqc_uncompute-dag-1",
            "name": self.name,
            "registers": [r.to_dict() for r in self.registers],
            "nodes": [n.to_dict() for n in self.nodes],
            "metadata": self.metadata,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "ArithDAG":
        if d.get("schema") not in (None, "ftqc_uncompute-dag-1"):
            raise ValueError(f"Unsupported ArithDAG schema: {d.get('schema')!r}")
        return cls(
            name=str(d.get("name", "arith_dag")),
            registers=tuple(RegDecl.from_dict(r) for r in d.get("registers", ())),
            nodes=tuple(DAGNode.from_dict(n) for n in d.get("nodes", ())),
            metadata=dict(d.get("metadata", {})),
        )

    @classmethod
    def from_json(cls, text: str) -> "ArithDAG":
        def _reject_float(val: str) -> None:
            raise ValueError(f"Float literal {val!r} forbidden in ArithDAG JSON; use integer or string Fraction")

        return cls.from_dict(json.loads(text, parse_float=_reject_float))

    def estimate_costs(self) -> dict[str, int]:
        """Fast Bloq-style analytical cost query without flattening large registers."""
        rmap = self.reg_map()
        toffolis = 0
        t_gates = 0
        rotations = 0
        hmr_count = 0
        depth_steps = 0
        for node in self.nodes:
            depth_steps += 1
            kind = node.op_kind
            mode = node.uncompute_mode
            if kind == "and_gate":
                if node.chirality == "uncompute" and mode in ("hmr", "ghost"):
                    hmr_count += 1
                else:
                    toffolis += 1
            elif kind == "gf2_mul":
                out_w = rmap[node.outputs[0]].width
                if node.chirality == "uncompute" and mode in ("hmr", "ghost"):
                    hmr_count += out_w
                else:
                    toffolis += int(node.params.get("toffolis", out_w * out_w))
            elif kind == "qrom_lookup":
                table = node.params.get("table", [])
                n_rows = int(node.params.get("num_rows", len(table)))
                out_w = rmap[node.outputs[0]].width
                if node.chirality == "uncompute" and mode == "hmr":
                    hmr_count += out_w
                    toffolis += max(0, (n_rows + 1) // 2 - 1)
                else:
                    toffolis += max(0, n_rows - 2)
            elif kind in ("int_add", "int_sub"):
                out_w = rmap[node.outputs[0]].width
                use_hmr = bool(node.params.get("use_hmr", mode == "hmr"))
                if use_hmr:
                    toffolis += max(0, out_w - 1)
                    hmr_count += max(0, out_w - 1)
                else:
                    toffolis += max(0, 2 * out_w - 2)
            elif kind in ("mod_mul", "mod_sqr", "fib_step", "fib_inv_step"):
                if node.chirality == "uncompute" and mode == "ghost":
                    out_w = rmap[node.outputs[0]].width
                    hmr_count += out_w
                elif node.chirality == "uncompute" and mode == "masked":
                    pass  # eliminated by superposition masking!
                else:
                    out_w = rmap[node.outputs[0]].width
                    toffolis += int(node.params.get("toffolis", 2 * out_w * out_w))
            elif kind == "phase_by_bits":
                in_w = rmap[node.inputs[0]].width
                if mode == "masked":
                    # Phase-gradient addition: (in_w - 1) Toffolis, 0 Z_POW syntheses!
                    toffolis += max(0, in_w - 1)
                    hmr_count += max(0, in_w - 1)
                else:
                    # Bitwise Z_POW(2^(j - in_w)): top 2 bits are Z and S, 3rd is T, rest rotations
                    t_gates += 1 if in_w >= 3 else 0
                    rotations += max(0, in_w - 3)
        total_qubits = sum(r.width for r in self.registers if r.kind != "classical")
        return {
            "toffolis": toffolis,
            "t_gates": t_gates,
            "rotations": rotations,
            "hmr_measurements": hmr_count,
            "total_declared_qubits": total_qubits,
            "dag_steps": depth_steps,
        }

    def lower_to_kmx(self) -> KickMixCircuit:
        """Lower concrete finite-width nodes (`and_gate`, `gf2_mul`, `qrom_lookup`, `int_add/sub`, `phase_by_bits`) to `.kmx`."""
        b = CircuitBuilder()
        b.set_meta("dag_name", self.name)
        q_regs: dict[str, list[Any]] = {}
        c_regs: dict[str, list[Any]] = {}
        rmap = self.reg_map()

        for r in self.registers:
            if r.kind == "classical":
                c_regs[r.name] = b.create_classical_register(r.width, name=r.name)
            elif r.role == "scratch":
                q_regs[r.name] = b.alloc_qubits(r.width)
            else:
                q_regs[r.name] = b.create_quantum_register(r.width, name=r.name)
                if r.kind == "mask" or r.mask_kind:
                    mk = r.mask_kind or "uniform_mod_n"
                    kw: dict[str, int | str] = {"width": r.width}
                    if r.modulus is not None:
                        kw["modulus"] = r.modulus
                    b.declare_mask(r.name, mk, **kw)

        for node in self.nodes:
            tag = f"dag:{node.node_id}:{node.op_kind}:{node.chirality}:{node.uncompute_mode}"
            with b.mark(tag):
                kind = node.op_kind
                if kind == "and_gate":
                    c1 = q_regs[node.inputs[0]][int(node.params.get("bit1", 0))]
                    c2 = q_regs[node.inputs[1]][int(node.params.get("bit2", 0))]
                    tgt = q_regs[node.outputs[0]][int(node.params.get("target_bit", 0))]
                    if node.chirality == "uncompute" and node.uncompute_mode == "hmr":
                        b.del_and(c1, c2, tgt)
                    else:
                        b.ccx(c1, c2, tgt)
                elif kind == "gf2_mul":
                    lhs = q_regs[node.inputs[0]]
                    rhs = q_regs[node.inputs[1]]
                    tgt = q_regs[node.outputs[0]]
                    r_spec = rmap[node.outputs[0]]
                    fld = GF2Field(r_spec.width, modulus=r_spec.gf_modulus)
                    if node.chirality == "uncompute" and node.uncompute_mode == "hmr":
                        b.del_gf2_mul(lhs, rhs, tgt, fld)
                    else:
                        b.init_gf2_mul(lhs, rhs, tgt, fld)
                elif kind == "qrom_lookup":
                    addr = q_regs[node.inputs[0]]
                    tgt = q_regs[node.outputs[0]]
                    tbl = [int(v) for v in node.params["table"]]
                    if node.chirality == "uncompute" and node.uncompute_mode == "hmr":
                        b.del_lookup(tbl, addr, tgt)
                    else:
                        b.init_lookup(tbl, addr, tgt)
                elif kind == "int_add":
                    lhs = q_regs[node.inputs[0]]
                    tgt = q_regs[node.outputs[0]]
                    use_hmr = bool(node.params.get("use_hmr", node.uncompute_mode == "hmr"))
                    b.iadd(lhs, tgt, use_hmr=use_hmr)
                elif kind == "int_sub":
                    lhs = q_regs[node.inputs[0]]
                    tgt = q_regs[node.outputs[0]]
                    use_hmr = bool(node.params.get("use_hmr", node.uncompute_mode == "hmr"))
                    b.isub(lhs, tgt, use_hmr=use_hmr)
                elif kind == "phase_by_bits":
                    reg = q_regs[node.inputs[0]]
                    f_bits = len(reg)
                    if node.uncompute_mode == "masked":
                        grad_name = node.params["gradient_reg"]
                        b.isub(reg, q_regs[grad_name], use_hmr=True)
                    else:
                        for j in range(f_bits):
                            # 2^(j - f_bits) turns == 2^(j - f_bits + 1) half-turns
                            exp = Fraction(1, 1 << (f_bits - 1 - j)) if (f_bits - 1 - j) >= 0 else Fraction(2, 1)
                            b.z_pow(reg[j], exp)
                else:
                    raise ValueError(f"DAG node {node.node_id!r} ({kind}) is abstract/pebbling-level and needs a concrete subroutine lowerer")

        return b.finish_circuit()


def lift_scoped_kmx_to_dag(circuit: KickMixCircuit) -> list[dict[str, Any]]:
    """Extract top-level `@enter:`/`@exit:` blocks and symmetric CCX patterns from a `KickMixCircuit`."""
    blocks: list[dict[str, Any]] = []
    stack: list[tuple[str, int]] = []
    for idx, op in enumerate(circuit.ops):
        if op.kind == "INFO_TAG" and op.info_tag:
            if op.info_tag.startswith("@enter:"):
                stack.append((op.info_tag[7:], idx))
            elif op.info_tag.startswith("@exit:") and stack:
                name, start_idx = stack.pop()
                if not stack:
                    sub_ops = circuit.ops[start_idx + 1 : idx]
                    toff = sum(1 for x in sub_ops if x.kind in TOFFOLI_LIKE_OPS)
                    hmr = sum(1 for x in sub_ops if x.kind == "HMR")
                    blocks.append(
                        {
                            "scope": name,
                            "start_op": start_idx,
                            "end_op": idx,
                            "toffoli_like": toff,
                            "hmr": hmr,
                        }
                    )
    return blocks
