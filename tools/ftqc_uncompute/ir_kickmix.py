"""Clean-room KickMix `.kmx` intermediate representation, parser, emitter, and builder.

Implements the published KickMix instruction set and subroutine patterns from:
  - Khattar et al., arXiv:2610.23832v2 (Section II.A-II.B)
  - Gidney, Quantum 2, 74 (2018) [temporary logical-AND & ripple-carry adder]

100% Python standard library, MIT-licensed, zero float literals or float calls.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
import re
from typing import Iterable, Iterator, Sequence

from .gf2 import GF2Field, gf2_lup_decompose


@dataclass(frozen=True, order=True)
class Qubit:
    id: int

    def __post_init__(self) -> None:
        if type(self.id) is not int or self.id < 0:
            raise ValueError(f"Qubit id must be a non-negative int, got {self.id!r}")

    def __str__(self) -> str:
        return f"q{self.id}"


@dataclass(frozen=True, order=True)
class Bit:
    id: int

    def __post_init__(self) -> None:
        if type(self.id) is not int or self.id < 0:
            raise ValueError(f"Bit id must be a non-negative int, got {self.id!r}")

    def __str__(self) -> str:
        return f"b{self.id}"


@dataclass(frozen=True)
class RegisterSpec:
    reg_id: int
    name: str | None
    kind: str  # "q" or "b"
    wires: tuple[Qubit | Bit, ...]

    def __post_init__(self) -> None:
        if self.kind not in ("q", "b"):
            raise ValueError(f"Register kind must be 'q' or 'b', got {self.kind!r}")
        expected_type = Qubit if self.kind == "q" else Bit
        if not all(isinstance(w, expected_type) for w in self.wires):
            raise ValueError(f"Register {self.name!r} wires must all be {expected_type.__name__}")


REVERSIBLE_OPS = frozenset({"X", "CX", "CCX", "SWAP", "CSWAP"})
DIAGONAL_OPS = frozenset({"Z", "CZ", "CCZ", "Z_POW", "NEG"})
MEAS_RESET_OPS = frozenset({"HMR", "R"})
BIT_OPS = frozenset({"BIT_STORE0", "BIT_STORE1", "BIT_INVERT"})
COND_STACK_OPS = frozenset({"PUSH_CONDITION", "POP_CONDITION"})
ANNOTATION_OPS = frozenset({"INFO_TAG", "DEBUG_PRINT"})

ALL_OPCODES = (
    REVERSIBLE_OPS
    | DIAGONAL_OPS
    | MEAS_RESET_OPS
    | BIT_OPS
    | COND_STACK_OPS
    | ANNOTATION_OPS
)

TOFFOLI_LIKE_OPS = frozenset({"CCX", "CCZ", "CSWAP"})


def normalize_half_turns(angle: Fraction | int | str) -> Fraction:
    """Normalize a half-turn phase exponent modulo 2 (1 full turn = 2 half-turns)."""
    if isinstance(angle, float):
        raise TypeError("float phase exponents are forbidden; pass Fraction, int, or decimal/rational str")
    val = Fraction(angle)
    reduced = val % Fraction(2, 1)
    if reduced > Fraction(1, 1):
        reduced -= Fraction(2, 1)
    return reduced


def format_half_turns(angle: Fraction) -> str:
    """Format a Fraction half-turn exponent deterministically (decimal if power-of-2 denom, else p/q)."""
    if angle.denominator == 1:
        return str(angle.numerator)
    den = angle.denominator
    # Check if denominator divides 10^k (exact terminating decimal, e.g. dyadic fractions 1/2, 1/4, 1/8...)
    tmp = den
    pow2 = 0
    pow5 = 0
    while tmp % 2 == 0:
        tmp //= 2
        pow2 += 1
    while tmp % 5 == 0:
        tmp //= 5
        pow5 += 1
    if tmp == 1 and max(pow2, pow5) <= 24:
        k = max(pow2, pow5)
        scaled = angle.numerator * (2 ** (k - pow2)) * (5 ** (k - pow5))
        sign = "-" if scaled < 0 else ""
        digits = str(abs(scaled)).rjust(k + 1, "0")
        whole, frac = digits[:-k], digits[-k:].rstrip("0")
        return f"{sign}{whole}.{frac}"
    return f"{angle.numerator}/{angle.denominator}"


def z_pow_hierarchy_level(half_turns: Fraction) -> int:
    """Clifford hierarchy level of `Z ** half_turns` (0=I, 1=Z, 2=S, 3=T, >=4=finer rotation)."""
    red = half_turns % Fraction(2, 1)
    if red == 0:
        return 0
    if red == 1:
        return 1
    num, den = red.numerator, red.denominator
    if num % 2 == 0:
        return 4
    # If den == 2^(k-1) for k >= 2, level is k (den=2 -> S -> level 2; den=4 -> T -> level 3)
    if den > 0 and (den & (den - 1)) == 0:
        return den.bit_length()
    return 4


@dataclass(frozen=True)
class KickMixOp:
    kind: str
    qubits: tuple[int, ...] = ()
    bits: tuple[int, ...] = ()
    angle: Fraction | None = None
    condition: int | None = None
    info_tag: str | None = None

    def __post_init__(self) -> None:
        if self.kind not in ALL_OPCODES:
            raise ValueError(f"Unknown KickMix opcode: {self.kind!r}")
        if any(q < 0 for q in self.qubits):
            raise ValueError("Negative qubit index")
        if any(b < 0 for b in self.bits):
            raise ValueError("Negative classical bit index")
        if self.condition is not None and self.condition < 0:
            raise ValueError("Negative condition bit index")
        if len(set(self.qubits)) != len(self.qubits):
            raise ValueError(f"{self.kind} cannot act on duplicate qubits {self.qubits}")
        # Arity checks
        q_arity = {
            "NEG": 0, "POP_CONDITION": 0, "PUSH_CONDITION": 0, "INFO_TAG": 0,
            "BIT_STORE0": 0, "BIT_STORE1": 0, "BIT_INVERT": 0,
            "X": 1, "Z": 1, "R": 1, "HMR": 1, "Z_POW": 1,
            "CX": 2, "CZ": 2, "SWAP": 2,
            "CCX": 3, "CCZ": 3, "CSWAP": 3,
        }
        b_arity = {
            "HMR": 1,
            "BIT_STORE0": 1, "BIT_STORE1": 1, "BIT_INVERT": 1,
        }
        if self.kind in q_arity and len(self.qubits) != q_arity[self.kind]:
            raise ValueError(f"{self.kind} expects {q_arity[self.kind]} qubits, got {len(self.qubits)}")
        if self.kind != "DEBUG_PRINT":
            expected_b = b_arity.get(self.kind, 0)
            if len(self.bits) != expected_b:
                raise ValueError(f"{self.kind} expects {expected_b} target bits, got {len(self.bits)}")
        if self.kind == "Z_POW":
            if self.angle is None:
                raise ValueError("Z_POW requires an angle in half-turns")
            object.__setattr__(self, "angle", normalize_half_turns(self.angle))
        elif self.angle is not None:
            raise ValueError(f"{self.kind} does not take an angle")
        if self.kind == "PUSH_CONDITION" and self.condition is None:
            raise ValueError("PUSH_CONDITION requires `condition` bit")
        if self.kind in ("POP_CONDITION", "INFO_TAG") and self.condition is not None:
            raise ValueError(f"{self.kind} cannot carry a per-op condition")
        if self.kind == "INFO_TAG" and self.info_tag is None:
            raise ValueError("INFO_TAG requires info_tag text")

    def to_kmx_line(self) -> str:
        if self.kind == "INFO_TAG":
            escaped = (self.info_tag or "").replace("\\", "\\\\").replace('"', '\\"')
            return f'INFO_TAG "{escaped}"'
        if self.kind == "POP_CONDITION":
            return "POP_CONDITION"
        if self.kind == "PUSH_CONDITION":
            return f"PUSH_CONDITION if b{self.condition}"
        parts = [self.kind]
        parts.extend(f"q{q}" for q in self.qubits)
        parts.extend(f"b{b}" for b in self.bits)
        if self.angle is not None:
            parts.append(format_half_turns(self.angle))
        if self.condition is not None:
            parts.append(f"if b{self.condition}")
        return " ".join(parts)


@dataclass(frozen=True)
class OperationStatistics:
    counts: dict[str, int]
    num_toffoli_like: int
    num_t_like: int
    num_rotations: int
    num_hmr: int
    num_resets: int
    num_clifford_2q: int

    def to_dict(self) -> dict[str, object]:
        return {
            "counts": dict(sorted(self.counts.items())),
            "num_toffoli_like": self.num_toffoli_like,
            "num_t_like": self.num_t_like,
            "num_rotations": self.num_rotations,
            "num_hmr": self.num_hmr,
            "num_resets": self.num_resets,
            "num_clifford_2q": self.num_clifford_2q,
        }


@dataclass(frozen=True)
class MaskDirective:
    """Structured `# @mask` annotation for superposition-masked registers."""

    register: str
    kind: str  # "mult_mod_n", "uniform_mod_n", "add_mod_2w", "phase_grad", "carry_runway"
    params: dict[str, int | str] = field(default_factory=dict)


@dataclass(frozen=True)
class KickMixCircuit:
    ops: tuple[KickMixOp, ...]
    registers: tuple[RegisterSpec, ...] = ()
    num_qubits: int = 0
    num_bits: int = 0
    masks: tuple[MaskDirective, ...] = ()
    metadata: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        names = [r.name if r.name is not None else f"r{r.reg_id}" for r in self.registers]
        ids = [r.reg_id for r in self.registers]
        wires = [w for r in self.registers for w in r.wires]
        if len(set(names)) != len(names) or len(set(ids)) != len(ids) or len(set(wires)) != len(wires):
            raise ValueError("Register names, IDs and wire assignments must be unique")
        max_q = -1
        max_b = -1
        for r in self.registers:
            for w in r.wires:
                if isinstance(w, Qubit):
                    max_q = max(max_q, w.id)
                else:
                    max_b = max(max_b, w.id)
        for op in self.ops:
            for q in op.qubits:
                max_q = max(max_q, q)
            for b in op.bits:
                max_b = max(max_b, b)
            if op.condition is not None:
                max_b = max(max_b, op.condition)
        nq = max(self.num_qubits, max_q + 1)
        nb = max(self.num_bits, max_b + 1)
        object.__setattr__(self, "num_qubits", nq)
        object.__setattr__(self, "num_bits", nb)

    def __len__(self) -> int:
        return len(self.ops)

    @property
    def register_data(self) -> dict[str, tuple[Qubit | Bit, ...]]:
        out: dict[str, tuple[Qubit | Bit, ...]] = {}
        for r in self.registers:
            key = r.name if r.name is not None else f"r{r.reg_id}"
            out[key] = r.wires
        return out

    def public_qubits(self) -> frozenset[int]:
        if not self.registers:
            return frozenset(range(self.num_qubits))
        return frozenset(
            w.id for r in self.registers if r.kind == "q" for w in r.wires if isinstance(w, Qubit)
        )

    def scratch_qubits(self) -> tuple[int, ...]:
        pub = self.public_qubits()
        return tuple(q for q in range(self.num_qubits) if q not in pub)

    def operation_statistics(self) -> OperationStatistics:
        counts: dict[str, int] = {}
        toff = 0
        t_like = 0
        rots = 0
        hmr_c = 0
        rst_c = 0
        c2q = 0
        for op in self.ops:
            counts[op.kind] = counts.get(op.kind, 0) + 1
            if op.kind in TOFFOLI_LIKE_OPS:
                toff += 1
            elif op.kind == "Z_POW" and op.angle is not None:
                lvl = z_pow_hierarchy_level(op.angle)
                if lvl == 3:
                    t_like += 1
                elif lvl >= 4:
                    rots += 1
            elif op.kind == "HMR":
                hmr_c += 1
            elif op.kind == "R":
                rst_c += 1
            elif op.kind in ("CX", "CZ", "SWAP"):
                c2q += 1
        return OperationStatistics(
            counts=counts,
            num_toffoli_like=toff,
            num_t_like=t_like,
            num_rotations=rots,
            num_hmr=hmr_c,
            num_resets=rst_c,
            num_clifford_2q=c2q,
        )

    def reaction_depth(self) -> int:
        """Causal classical-feedback depth across HMR measurements and conditioned gates."""
        q_ready = [0] * self.num_qubits
        b_ready = [0] * self.num_bits
        cond_stack: list[int] = []
        max_depth = 0
        for op in self.ops:
            if op.kind == "INFO_TAG" or op.kind == "DEBUG_PRINT":
                continue
            if op.kind == "PUSH_CONDITION":
                assert op.condition is not None
                cond_stack.append(op.condition)
                continue
            if op.kind == "POP_CONDITION":
                if cond_stack:
                    cond_stack.pop()
                continue
            active_conds = list(cond_stack)
            if op.condition is not None:
                active_conds.append(op.condition)
            cond_lvl = max((b_ready[c] for c in active_conds), default=0)
            if op.kind == "HMR":
                q = op.qubits[0]
                b = op.bits[0]
                lvl = max(q_ready[q], cond_lvl) + 1
                b_ready[b] = lvl
                q_ready[q] = lvl
                max_depth = max(max_depth, lvl)
            elif op.kind in BIT_OPS:
                b = op.bits[0]
                prev = b_ready[b] if op.kind == "BIT_INVERT" else 0
                b_ready[b] = max(prev, cond_lvl)
            else:
                q_lvl = max((q_ready[q] for q in op.qubits), default=0)
                lvl = max(q_lvl, cond_lvl)
                for q in op.qubits:
                    q_ready[q] = lvl
        return max_depth

    def scope_costs(self) -> list[dict[str, object]]:
        """Hierarchical cost breakdown from `INFO_TAG "@enter:<name>"` / `"@exit:<name>"`."""
        stack: list[tuple[str, int, set[int]]] = []
        records: list[dict[str, object]] = []
        for op in self.ops:
            if op.kind == "INFO_TAG" and op.info_tag:
                if op.info_tag.startswith("@enter:"):
                    stack.append((op.info_tag[7:], 0, set()))
                elif op.info_tag.startswith("@exit:") and stack:
                    name, toffs, qs = stack.pop()
                    records.append({"scope": name, "toffoli_like": toffs, "qubits_touched": len(qs)})
                continue
            is_toff = 1 if op.kind in TOFFOLI_LIKE_OPS else 0
            for idx in range(len(stack)):
                nm, tf, qs = stack[idx]
                qs.update(op.qubits)
                stack[idx] = (nm, tf + is_toff, qs)
        return records

    def to_kmx(self, include_directives: bool = True) -> str:
        lines: list[str] = []
        if include_directives:
            for k, v in sorted(self.metadata.items()):
                lines.append(f"# @meta {k}={v}")
            for m in self.masks:
                kv = " ".join(f"{pk}={pv}" for pk, pv in sorted(m.params.items()))
                suffix = f" {kv}" if kv else ""
                lines.append(f"# @mask {m.register} kind={m.kind}{suffix}")
        for r in self.registers:
            if r.name is not None:
                escaped = r.name.replace("\\", "\\\\").replace('"', '\\"')
                lines.append(f'REGISTER r{r.reg_id} "{escaped}"')
            else:
                lines.append(f"REGISTER r{r.reg_id}")
            for w in r.wires:
                lines.append(f"APPEND_TO_REGISTER {w} r{r.reg_id}")
        for op in self.ops:
            lines.append(op.to_kmx_line())
        return "\n".join(lines) + ("\n" if lines else "")

    def canonical_bytes(self) -> bytes:
        """Comment-independent canonical serialization for SHA-256 Fiat-Shamir seeding."""
        lines: list[str] = []
        for r in self.registers:
            if r.name is not None:
                lines.append(f'REGISTER r{r.reg_id} "{r.name}"')
            else:
                lines.append(f"REGISTER r{r.reg_id}")
            for w in r.wires:
                lines.append(f"APPEND_TO_REGISTER {w} r{r.reg_id}")
        for mask in self.masks:
            lines.append("# @mask " + json.dumps([mask.register, mask.kind, mask.params], sort_keys=True))
        for op in self.ops:
            if op.kind not in ANNOTATION_OPS:
                lines.append(op.to_kmx_line())
        return ("\n".join(lines) + "\n").encode("utf-8")

    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_bytes()).hexdigest()

    @classmethod
    def from_kmx(cls, text: str) -> "KickMixCircuit":
        reg_names: dict[int, str | None] = {}
        reg_wires: dict[int, list[Qubit | Bit]] = {}
        reg_order: list[int] = []
        ops: list[KickMixOp] = []
        masks: list[MaskDirective] = []
        metadata: dict[str, str] = {}

        reg_re = re.compile(r'^REGISTER\s+r(\d+)(?:\s+"((?:\\.|[^"\\])*)")?$')
        app_re = re.compile(r"^APPEND_TO_REGISTER\s+([qb])(\d+)\s+r(\d+)$")
        info_re = re.compile(r'^INFO_TAG\s+"((?:\\.|[^"\\])*)"$')

        for lineno, raw_line in enumerate(text.splitlines(), start=1):
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                body = stripped[1:].strip()
                if body.startswith("@mask "):
                    toks = body[6:].split()
                    if len(toks) < 2:
                        raise ValueError(f"Line {lineno}: incomplete mask directive")
                    params = {}
                    for token in toks[1:]:
                        key, sep, value = token.partition("=")
                        if not sep or not value or not key or key in params:
                            raise ValueError(f"Line {lineno}: invalid mask parameter {token!r}")
                        params[key] = int(value) if value.lstrip("-").isdigit() else value
                    masks.append(MaskDirective(toks[0], params.pop("kind", "uniform_mod_n"), params))
                elif body.startswith("@meta ") and "=" in body[6:]:
                    k, v = body[6:].split("=", 1)
                    metadata[k.strip()] = v.strip()
                continue

            # Strip inline `#` comments outside quotes
            if not stripped.startswith(("REGISTER", "INFO_TAG")) and "#" in stripped:
                stripped = stripped.split("#", 1)[0].strip()

            m_reg = reg_re.match(stripped)
            if m_reg:
                rid = int(m_reg.group(1))
                rname = m_reg.group(2)
                if rname is not None:
                    rname = rname.replace('\\"', '"').replace("\\\\", "\\")
                if rid in reg_names:
                    raise ValueError(f"Line {lineno}: duplicate REGISTER r{rid}")
                reg_order.append(rid)
                reg_wires[rid] = []
                reg_names[rid] = rname
                continue

            m_app = app_re.match(stripped)
            if m_app:
                wkind, wid, rid = m_app.group(1), int(m_app.group(2)), int(m_app.group(3))
                if rid not in reg_names:
                    raise ValueError(f"Line {lineno}: APPEND_TO_REGISTER to undeclared r{rid}")
                wire: Qubit | Bit = Qubit(wid) if wkind == "q" else Bit(wid)
                reg_wires[rid].append(wire)
                continue

            m_info = info_re.match(stripped)
            if m_info:
                tag = m_info.group(1).replace('\\"', '"').replace("\\\\", "\\")
                ops.append(KickMixOp("INFO_TAG", info_tag=tag))
                continue

            tokens = stripped.split()
            cond: int | None = None
            if len(tokens) >= 3 and tokens[-2] == "if":
                ctok = tokens[-1]
                if not ctok.startswith("b") or not ctok[1:].isdigit():
                    raise ValueError(f"Line {lineno}: invalid condition token {ctok!r}")
                cond = int(ctok[1:])
                tokens = tokens[:-2]

            opcode = tokens[0]
            if opcode not in ALL_OPCODES:
                raise ValueError(f"Line {lineno}: unknown operation {opcode!r}")

            qs: list[int] = []
            bs: list[int] = []
            angle: Fraction | None = None
            for tok in tokens[1:]:
                if tok.startswith("q") and tok[1:].isdigit():
                    qs.append(int(tok[1:]))
                elif tok.startswith("b") and tok[1:].isdigit():
                    bs.append(int(tok[1:]))
                elif opcode == "Z_POW" and angle is None:
                    angle = normalize_half_turns(tok)
                else:
                    raise ValueError(f"Line {lineno}: unexpected operand {tok!r} for {opcode}")

            ops.append(
                KickMixOp(
                    kind=opcode,
                    qubits=tuple(qs),
                    bits=tuple(bs),
                    angle=angle,
                    condition=cond,
                )
            )

        registers: list[RegisterSpec] = []
        for rid in reg_order:
            ws = reg_wires[rid]
            if not ws:
                raise ValueError(f"Register r{rid} has no wires")
            rk = "q" if isinstance(ws[0], Qubit) else "b"
            registers.append(RegisterSpec(reg_id=rid, name=reg_names[rid], kind=rk, wires=tuple(ws)))

        return cls(ops=tuple(ops), registers=tuple(registers), masks=tuple(masks), metadata=metadata)


class CircuitBuilder:
    """Ergonomic clean-room builder for KickMix circuits and MBUC arithmetic subroutines."""

    def __init__(self, max_qubits: int | None = None) -> None:
        self._max_qubits = max_qubits
        self._next_q = 0
        self._next_b = 0
        self._free_qs: list[int] = []
        self._free_bs: list[int] = []
        self._live_qs: set[int] = set()
        self._live_bs: set[int] = set()
        self._registers: list[RegisterSpec] = []
        self._ops: list[KickMixOp] = []
        self._masks: list[MaskDirective] = []
        self._metadata: dict[str, str] = {}

    def set_meta(self, key: str, value: str) -> None:
        self._metadata[key] = value

    def declare_mask(self, register: str, kind: str, **params: int | str) -> None:
        self._masks.append(MaskDirective(register=register, kind=kind, params=dict(params)))

    def alloc_qubits(self, count: int) -> list[Qubit]:
        if count < 0:
            raise ValueError("count must be >= 0")
        out: list[Qubit] = []
        for _ in range(count):
            if self._free_qs:
                qid = self._free_qs.pop(0)
            else:
                if self._max_qubits is not None and self._next_q >= self._max_qubits:
                    raise RuntimeError(f"Exceeded max_qubits budget ({self._max_qubits})")
                qid = self._next_q
                self._next_q += 1
            self._live_qs.add(qid)
            out.append(Qubit(qid))
        return out

    def alloc_bits(self, count: int) -> list[Bit]:
        if count < 0:
            raise ValueError("count must be >= 0")
        out: list[Bit] = []
        for _ in range(count):
            if self._free_bs:
                bid = self._free_bs.pop(0)
            else:
                bid = self._next_b
                self._next_b += 1
            self._live_bs.add(bid)
            out.append(Bit(bid))
        return out

    def free(self, wires: Qubit | Bit | Iterable[Qubit | Bit]) -> None:
        items = [wires] if isinstance(wires, (Qubit, Bit)) else list(wires)
        for w in items:
            if isinstance(w, Qubit):
                if w.id not in self._live_qs:
                    raise ValueError(f"Double or invalid free of {w}")
                self._live_qs.remove(w.id)
                self._free_qs.append(w.id)
                self._free_qs.sort()
            else:
                if w.id not in self._live_bs:
                    raise ValueError(f"Double or invalid free of {w}")
                self._live_bs.remove(w.id)
                self._free_bs.append(w.id)
                self._free_bs.sort()

    def create_quantum_register(self, width: int, name: str | None = None) -> list[Qubit]:
        qs = self.alloc_qubits(width)
        rid = len(self._registers)
        self._registers.append(RegisterSpec(reg_id=rid, name=name, kind="q", wires=tuple(qs)))
        return qs

    def create_classical_register(self, width: int, name: str | None = None) -> list[Bit]:
        bs = self.alloc_bits(width)
        rid = len(self._registers)
        self._registers.append(RegisterSpec(reg_id=rid, name=name, kind="b", wires=tuple(bs)))
        return bs

    @contextmanager
    def mark(self, scope_name: str) -> Iterator[None]:
        self._ops.append(KickMixOp("INFO_TAG", info_tag=f"@enter:{scope_name}"))
        try:
            yield
        finally:
            self._ops.append(KickMixOp("INFO_TAG", info_tag=f"@exit:{scope_name}"))

    @contextmanager
    def condition(self, cond: Bit | bool) -> Iterator[None]:
        if isinstance(cond, bool):
            if cond:
                yield
            return
        self._ops.append(KickMixOp("PUSH_CONDITION", condition=cond.id))
        try:
            yield
        finally:
            self._ops.append(KickMixOp("POP_CONDITION"))

    @staticmethod
    def _cond_id(condition: Bit | bool | None) -> tuple[bool, int | None]:
        if condition is None or condition is True:
            return True, None
        if condition is False:
            return False, None
        if isinstance(condition, Bit):
            return True, condition.id
        raise TypeError(f"Unsupported condition {condition!r}")

    # Primitive instructions with automatic classical downgrading
    def x(self, target: Qubit | Sequence[Qubit], *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if not active:
            return
        ts = [target] if isinstance(target, Qubit) else list(target)
        for q in ts:
            self._ops.append(KickMixOp("X", qubits=(q.id,), condition=cid))

    def z(self, target: Qubit | Sequence[Qubit], *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if not active:
            return
        ts = [target] if isinstance(target, Qubit) else list(target)
        for q in ts:
            self._ops.append(KickMixOp("Z", qubits=(q.id,), condition=cid))

    def neg(self, *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if active:
            self._ops.append(KickMixOp("NEG", condition=cid))

    def z_pow(
        self,
        target: Qubit,
        half_turns: Fraction | int | str,
        *,
        condition: Bit | bool | None = None,
    ) -> None:
        active, cid = self._cond_id(condition)
        if not active:
            return
        ang = normalize_half_turns(half_turns)
        if ang == 0:
            return
        if ang == 1 or ang == -1:
            self.z(target, condition=condition)
            return
        self._ops.append(KickMixOp("Z_POW", qubits=(target.id,), angle=ang, condition=cid))

    def cx(
        self,
        control: Qubit | Bit | bool,
        target: Qubit,
        *,
        condition: Bit | bool | None = None,
    ) -> None:
        active, cid = self._cond_id(condition)
        if not active or control is False:
            return
        if control is True:
            self.x(target, condition=condition)
        elif isinstance(control, Bit):
            if cid is None:
                self.x(target, condition=control)
            else:
                with self.condition(control):
                    self.x(target, condition=condition)
        else:
            self._ops.append(KickMixOp("CX", qubits=(control.id, target.id), condition=cid))

    def cz(
        self,
        c1: Qubit | Bit | bool,
        c2: Qubit | Bit | bool,
        *,
        condition: Bit | bool | None = None,
    ) -> None:
        active, cid = self._cond_id(condition)
        if not active or c1 is False or c2 is False:
            return
        if c1 is True:
            if isinstance(c2, Qubit):
                self.z(c2, condition=condition)
            elif isinstance(c2, Bit):
                with self.condition(c2):
                    self.neg(condition=condition)
            else:
                self.neg(condition=condition)
            return
        if c2 is True:
            self.cz(True, c1, condition=condition)
            return
        if isinstance(c1, Bit):
            with self.condition(c1):
                self.cz(True, c2, condition=condition)
            return
        if isinstance(c2, Bit):
            with self.condition(c2):
                self.cz(True, c1, condition=condition)
            return
        self._ops.append(KickMixOp("CZ", qubits=(c1.id, c2.id), condition=cid))

    def ccx(
        self,
        c1: Qubit | Bit | bool,
        c2: Qubit | Bit | bool,
        target: Qubit,
        *,
        condition: Bit | bool | None = None,
    ) -> None:
        active, cid = self._cond_id(condition)
        if not active or c1 is False or c2 is False:
            return
        if c1 is True:
            self.cx(c2, target, condition=condition)
            return
        if c2 is True:
            self.cx(c1, target, condition=condition)
            return
        if isinstance(c1, Bit):
            with self.condition(c1):
                self.cx(c2, target, condition=condition)
            return
        if isinstance(c2, Bit):
            with self.condition(c2):
                self.cx(c1, target, condition=condition)
            return
        self._ops.append(KickMixOp("CCX", qubits=(c1.id, c2.id, target.id), condition=cid))

    def ccz(
        self,
        c1: Qubit,
        c2: Qubit,
        c3: Qubit,
        *,
        condition: Bit | bool | None = None,
    ) -> None:
        active, cid = self._cond_id(condition)
        if active:
            self._ops.append(KickMixOp("CCZ", qubits=(c1.id, c2.id, c3.id), condition=cid))

    def swap(self, q1: Qubit, q2: Qubit, *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if active and q1 != q2:
            self._ops.append(KickMixOp("SWAP", qubits=(q1.id, q2.id), condition=cid))

    def cswap(
        self,
        control: Qubit | Bit | bool,
        q1: Qubit,
        q2: Qubit,
        *,
        condition: Bit | bool | None = None,
    ) -> None:
        active, cid = self._cond_id(condition)
        if not active or control is False or q1 == q2:
            return
        if control is True:
            self.swap(q1, q2, condition=condition)
        elif isinstance(control, Bit):
            with self.condition(control):
                self.swap(q1, q2, condition=condition)
        else:
            self._ops.append(KickMixOp("CSWAP", qubits=(control.id, q1.id, q2.id), condition=cid))

    def reset(self, target: Qubit | Sequence[Qubit], *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if not active:
            return
        ts = [target] if isinstance(target, Qubit) else list(target)
        for q in ts:
            self._ops.append(KickMixOp("R", qubits=(q.id,), condition=cid))

    def hmr(self, target: Qubit, out_bit: Bit, *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if active:
            self._ops.append(KickMixOp("HMR", qubits=(target.id,), bits=(out_bit.id,), condition=cid))

    def bit_store(self, target: Bit, value: bool, *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if active:
            op = "BIT_STORE1" if value else "BIT_STORE0"
            self._ops.append(KickMixOp(op, bits=(target.id,), condition=cid))

    def bit_invert(self, target: Bit, *, condition: Bit | bool | None = None) -> None:
        active, cid = self._cond_id(condition)
        if active:
            self._ops.append(KickMixOp("BIT_INVERT", bits=(target.id,), condition=cid))

    # ------------------------------------------------------------------
    # Measurement-based uncomputation (MBUC) & arithmetic subroutines
    # ------------------------------------------------------------------

    def init_and(
        self,
        c1: Qubit | Bit | bool,
        c2: Qubit | Bit | bool,
        target: Qubit,
        *,
        explicit_reset: bool = False,
    ) -> None:
        """Compute `target := c1 & c2` on a clean `|0>` ancilla (1 Toffoli)."""
        with self.mark("init_and"):
            if explicit_reset:
                self.reset(target)
            self.ccx(c1, c2, target)

    def del_and(self, c1: Qubit | Bit | bool, c2: Qubit | Bit | bool, target: Qubit) -> None:
        """Uncompute `target == c1 & c2` back to `|0>` via HMR + conditioned CZ (0 Toffolis)."""
        with self.mark("del_and"):
            b = self.alloc_bits(1)[0]
            self.hmr(target, b)
            self.cz(c1, c2, condition=b)
            self.bit_store(b, False)
            self.free(b)

    def iadd(
        self,
        lhs: Sequence[Qubit] | int,
        target: Sequence[Qubit],
        *,
        control: Qubit | Bit | bool = True,
        use_hmr: bool = True,
    ) -> None:
        """In-place modular addition `target := (target + lhs) mod 2^w` (Gidney 2018).

        Uses `w - 1` Toffolis when uncontrolled with `use_hmr=True` (vs `2w - 2` symmetric).
        """
        w = len(target)
        if w == 0 or control is False:
            return
        with self.mark("iadd_gidney" if use_hmr else "iadd_symmetric"):
            if isinstance(lhs, int):
                val = lhs % (1 << w)
                if val == 0:
                    return
                tmp = self.alloc_qubits(w)
                for i in range(w):
                    if (val >> i) & 1:
                        self.x(tmp[i])
                self.iadd(tmp, target, control=control, use_hmr=use_hmr)
                for i in range(w):
                    if (val >> i) & 1:
                        self.x(tmp[i])
                self.free(tmp)
                return

            if len(lhs) != w:
                raise ValueError(f"iadd expects equal register widths, got {len(lhs)} and {w}")
            if w == 1:
                if control is True:
                    self.cx(lhs[0], target[0])
                else:
                    self.ccx(control, lhs[0], target[0])
                return

            carries = self.alloc_qubits(w - 1)
            # Forward carry ripple (Gidney 2018 Fig. 2)
            self.ccx(lhs[0], target[0], carries[0])
            for i in range(1, w - 1):
                self.cx(carries[i - 1], lhs[i])
                self.cx(carries[i - 1], target[i])
                self.ccx(lhs[i], target[i], carries[i])
                self.cx(carries[i - 1], carries[i])

            # Top bit sum
            if control is True:
                self.cx(carries[w - 2], target[w - 1])
                self.cx(lhs[w - 1], target[w - 1])
            else:
                self.ccx(control, carries[w - 2], target[w - 1])
                self.ccx(control, lhs[w - 1], target[w - 1])

            # Reverse carry sweep & sum extraction
            for i in range(w - 2, 0, -1):
                self.cx(carries[i - 1], carries[i])
                if use_hmr:
                    self.del_and(lhs[i], target[i], carries[i])
                else:
                    self.ccx(lhs[i], target[i], carries[i])
                self.cx(carries[i - 1], lhs[i])
                if control is True:
                    self.cx(lhs[i], target[i])
                else:
                    self.cx(carries[i - 1], target[i])
                    self.ccx(control, lhs[i], target[i])
                    self.ccx(control, carries[i - 1], target[i])

            if use_hmr:
                self.del_and(lhs[0], target[0], carries[0])
            else:
                self.ccx(lhs[0], target[0], carries[0])
            if control is True:
                self.cx(lhs[0], target[0])
            else:
                self.ccx(control, lhs[0], target[0])
            self.free(carries)

    def isub(
        self,
        lhs: Sequence[Qubit] | int,
        target: Sequence[Qubit],
        *,
        control: Qubit | Bit | bool = True,
        use_hmr: bool = True,
    ) -> None:
        """In-place modular subtraction `target := (target - lhs) mod 2^w` via two's complement."""
        w = len(target)
        if w == 0 or control is False:
            return
        with self.mark("isub"):
            if isinstance(lhs, int):
                neg_val = (-lhs) % (1 << w)
                self.iadd(neg_val, target, control=control, use_hmr=use_hmr)
                return
            self.x(target)
            self.iadd(lhs, target, control=control, use_hmr=use_hmr)
            self.x(target)

    def gf2_imul_const(self, factor: int, target: Sequence[Qubit], field: GF2Field) -> None:
        """In-place reversible constant multiplication in GF(2^m) using CNOT/SWAP LUP synthesis."""
        if factor == 0:
            raise ValueError("GF(2^m) in-place multiplication requires nonzero factor")
        if len(target) != field.degree:
            raise ValueError("target width must match field degree")
        rows = field.mul_const_matrix(factor)
        for kind, a, b in gf2_lup_decompose(rows, field.degree):
            if kind == "SWAP":
                self.swap(target[a], target[b])
            else:
                self.cx(target[a], target[b])

    def init_gf2_mul(
        self,
        lhs: Sequence[Qubit],
        rhs: Sequence[Qubit],
        target: Sequence[Qubit],
        field: GF2Field,
    ) -> None:
        """Clean out-of-place GF(2^m) product `target := (lhs * rhs) mod modulus`."""
        m = field.degree
        if not (len(lhs) == len(rhs) == len(target) == m):
            raise ValueError("GF(2^m) registers must all have width field.degree")
        with self.mark("init_gf2_mul"):
            tensor = field.mul_bilinear_tensor()
            for k in range(m):
                for i in range(m):
                    mask_j = tensor[k][i]
                    for j in range(m):
                        if (mask_j >> j) & 1:
                            self.ccx(lhs[i], rhs[j], target[k])

    def del_gf2_mul(
        self,
        lhs: Sequence[Qubit],
        rhs: Sequence[Qubit],
        target: Sequence[Qubit],
        field: GF2Field,
    ) -> None:
        """Zero-Toffoli MBUC cleanup of `target == (lhs * rhs) mod modulus` via HMR + CZ matrix."""
        m = field.degree
        if not (len(lhs) == len(rhs) == len(target) == m):
            raise ValueError("GF(2^m) registers must all have width field.degree")
        with self.mark("del_gf2_mul"):
            meas = self.alloc_bits(m)
            for k in range(m):
                self.hmr(target[k], meas[k])
            tensor = field.mul_bilinear_tensor()
            par_bit = self.alloc_bits(1)[0]
            for i in range(m):
                for j in range(m):
                    k_contributors = [k for k in range(m) if (tensor[k][i] >> j) & 1]
                    if not k_contributors:
                        continue
                    if len(k_contributors) == 1:
                        self.cz(lhs[i], rhs[j], condition=meas[k_contributors[0]])
                    else:
                        self.bit_store(par_bit, False)
                        for k in k_contributors:
                            self.bit_invert(par_bit, condition=meas[k])
                        self.cz(lhs[i], rhs[j], condition=par_bit)
                        self.bit_store(par_bit, False)
            for k in range(m):
                self.bit_store(meas[k], False)
            self.free(par_bit)
            self.free(meas)

    def init_lookup(
        self,
        table: Sequence[int],
        address: Sequence[Qubit],
        target: Sequence[Qubit],
    ) -> None:
        """Clean QROM table lookup `target ^= table[address]` via segment-tree unary iteration."""
        num_rows = len(table)
        w = len(target)
        n = len(address)
        if num_rows < 1 or num_rows > (1 << n):
            raise ValueError(f"table length {num_rows} invalid for {n} address qubits")
        with self.mark("init_lookup"):
            if n == 0:
                val = int(table[0])
                for b in range(w):
                    if (val >> b) & 1:
                        self.x(target[b])
                return

            padded = [int(table[r]) if r < num_rows else 0 for r in range(1 << n)]
            bits_msb = list(reversed(address))
            ancs = self.alloc_qubits(n - 1) if n > 1 else []

            def load_val(val: int, ctrl: Qubit | bool) -> None:
                for b in range(w):
                    if (val >> b) & 1:
                        self.cx(ctrl, target[b])

            def walk(ctrl: Qubit | bool, level: int, lo: int, hi: int) -> None:
                if all(padded[r] == 0 for r in range(lo, hi)):
                    return
                first = padded[lo]
                if all(padded[r] == first for r in range(lo, hi)):
                    load_val(first, ctrl)
                    return
                if hi - lo == 1:
                    load_val(first, ctrl)
                    return
                mid = (lo + hi) // 2
                left_zero = all(padded[r] == 0 for r in range(lo, mid))
                right_zero = all(padded[r] == 0 for r in range(mid, hi))
                abit = bits_msb[level]
                if ctrl is True:
                    if not left_zero:
                        self.x(abit)
                        walk(abit, level + 1, lo, mid)
                        self.x(abit)
                    if not right_zero:
                        walk(abit, level + 1, mid, hi)
                else:
                    anc = ancs[level - 1]
                    if right_zero:
                        self.x(abit)
                        self.init_and(ctrl, abit, anc)
                        walk(anc, level + 1, lo, mid)
                        self.del_and(ctrl, abit, anc)
                        self.x(abit)
                    elif left_zero:
                        self.init_and(ctrl, abit, anc)
                        walk(anc, level + 1, mid, hi)
                        self.del_and(ctrl, abit, anc)
                    else:
                        self.x(abit)
                        self.init_and(ctrl, abit, anc)
                        self.x(abit)
                        walk(anc, level + 1, lo, mid)
                        self.cx(ctrl, anc)
                        walk(anc, level + 1, mid, hi)
                        self.del_and(ctrl, abit, anc)

            walk(True, 0, 0, 1 << n)
            if ancs:
                self.free(ancs)

    def del_lookup(
        self,
        table: Sequence[int],
        address: Sequence[Qubit],
        target: Sequence[Qubit],
    ) -> None:
        """MBUC cleanup of `target == table[address]` via HMR + LSB-folded 1-bit sign lookup.

        Cuts Toffoli count by >= 2x relative to reversing `init_lookup` (Berry-Craig-Babbush 2019).
        """
        num_rows = len(table)
        w = len(target)
        n = len(address)
        if num_rows < 1 or num_rows > (1 << n):
            raise ValueError(f"table length {num_rows} invalid for {n} address qubits")
        with self.mark("del_lookup"):
            meas = self.alloc_bits(w)
            for b in range(w):
                self.hmr(target[b], meas[b])

            def emit_row_parity(row_idx: int, parity_bit: Bit) -> None:
                self.bit_store(parity_bit, False)
                val = int(table[row_idx])
                for b in range(w):
                    if (val >> b) & 1:
                        self.bit_invert(parity_bit, condition=meas[b])

            if n == 0:
                p0 = self.alloc_bits(1)[0]
                emit_row_parity(0, p0)
                self.neg(condition=p0)
                self.bit_store(p0, False)
                self.free(p0)
            elif n == 1:
                p0, p1 = self.alloc_bits(2)
                emit_row_parity(0, p0)
                if num_rows > 1:
                    emit_row_parity(1, p1)
                else:
                    self.bit_store(p1, False)
                # (-1)^(p0 * (1 - a0) + p1 * a0) = (-1)^p0 * (-1)^((p0 ^ p1) * a0)
                self.neg(condition=p0)
                self.bit_invert(p1, condition=p0)
                self.z(address[0], condition=p1)
                self.bit_store(p0, False)
                self.bit_store(p1, False)
                self.free([p0, p1])
            else:
                # Fold LSB `address[0]` into CZ/Z at leaves of an `(n-1)`-bit unary tree!
                lsb = address[0]
                high_msb = list(reversed(address[1:]))
                high_n = n - 1
                num_pairs = (num_rows + 1) // 2
                ancs = self.alloc_qubits(high_n - 1) if high_n > 1 else []
                p0, p1 = self.alloc_bits(2)

                def apply_pair_phase(pair_idx: int, ctrl: Qubit) -> None:
                    r0 = 2 * pair_idx
                    r1 = 2 * pair_idx + 1
                    emit_row_parity(r0, p0)
                    if r1 < num_rows:
                        emit_row_parity(r1, p1)
                    else:
                        self.bit_store(p1, False)
                    # Phase on ctrl==1 is (-1)^(p0 ^ ((p0 ^ p1) & lsb))
                    self.z(ctrl, condition=p0)
                    self.bit_invert(p1, condition=p0)
                    self.cz(ctrl, lsb, condition=p1)
                    self.bit_store(p0, False)
                    self.bit_store(p1, False)

                def walk_pairs(ctrl: Qubit | bool, level: int, lo: int, hi: int) -> None:
                    if lo >= num_pairs:
                        return
                    if hi - lo == 1:
                        assert isinstance(ctrl, Qubit)
                        apply_pair_phase(lo, ctrl)
                        return
                    mid = (lo + hi) // 2
                    abit = high_msb[level]
                    if ctrl is True:
                        self.x(abit)
                        walk_pairs(abit, level + 1, lo, mid)
                        self.x(abit)
                        if mid < num_pairs:
                            walk_pairs(abit, level + 1, mid, hi)
                    else:
                        if mid >= num_pairs:
                            anc = ancs[level - 1]
                            self.x(abit)
                            self.init_and(ctrl, abit, anc)
                            walk_pairs(anc, level + 1, lo, mid)
                            self.del_and(ctrl, abit, anc)
                            self.x(abit)
                            return
                        anc = ancs[level - 1]
                        self.x(abit)
                        self.init_and(ctrl, abit, anc)
                        self.x(abit)
                        walk_pairs(anc, level + 1, lo, mid)
                        self.cx(ctrl, anc)
                        walk_pairs(anc, level + 1, mid, hi)
                        self.del_and(ctrl, abit, anc)

                walk_pairs(True, 0, 0, 1 << high_n)
                self.free([p0, p1])
                if ancs:
                    self.free(ancs)

            for b in range(w):
                self.bit_store(meas[b], False)
            self.free(meas)

    def finish_circuit(self) -> KickMixCircuit:
        return KickMixCircuit(
            ops=tuple(self._ops),
            registers=tuple(self._registers),
            num_qubits=self._next_q,
            num_bits=self._next_b,
            masks=tuple(self._masks),
            metadata=dict(self._metadata),
        )
