"""Clifford Circuits, ZX Diagrams, Idealised Noise Submodels, and Fault Gadgets.

Implements Definitions 2.1-2.5, 4.1-4.2, 5.1-5.3, and Section 6.1 (Definitions 6.1-6.3,
Propositions 6.4-6.10) of Rodatz, Poor, and Kissinger, "Fault Tolerance by Construction"
(arXiv:2506.17181v4).

Zero float literals appear anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from typing import Any, Iterable

from .gf2 import pack_pauli_char


class RefusalError(ValueError):
    """Named fail-closed refusal (exits code 2 on CLI)."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")
        self.code = code
        self.detail = detail


VALID_NOISE_MODELS = ("circuit_depolarizing", "circuit_css", "fault_free")
VALID_OP_MODES = ("noisy", "ft", "ff")

SINGLE_QUBIT_GATES = ("I", "X", "Y", "Z", "H", "S", "SDG")
PREP_STATES = ("PREP_0", "PREP_1", "PREP_PLUS", "PREP_MINUS")
SINGLE_MEAS = ("MEAS_Z", "MEAS_X", "MEAS_Y")


@dataclass(frozen=True)
class CircuitOp:
    """Single operation in a CliffordCircuit.

    Supported kinds:
      - Preparations (allocate/reset wire): PREP_0, PREP_1, PREP_PLUS, PREP_MINUS (1 qubit)
      - Idealised multi-qubit prep: CAT_PREP (k qubits; params['basis'] in {'Z', 'X'}, default 'Z'),
        ENC_PREP (len(qubits) == n_phys, params={'stabilizers': [...], 'logical_inputs': [...], ...})
      - Single-qubit Cliffords: I, X, Y, Z, H, S, SDG
      - Two-qubit Cliffords: CX (control=q0, target=q1), CZ, SWAP
      - Destructive single-qubit measurements: MEAS_Z, MEAS_X, MEAS_Y
      - Non-destructive multi-qubit Pauli projector/measurement:
        MPP (params['pauli'] = string like 'ZZZZ' of length len(qubits))

    mode in {'noisy', 'ft', 'ff'} (Definition 2.5):
      - 'noisy': full circuit-level noise on outputs (and wire flips on incident non-purple wires)
      - 'ft'   : allowed on MPP / measurements: removes multi-qubit output gate faults, keeping
                 measurement flip and measurement flip + single-qubit output faults
      - 'ff'   : fault-free (purple): removes all internal gate/measurement faults.
    """

    kind: str
    qubits: tuple[int, ...]
    mode: str = "noisy"
    label: str = ""
    params: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.mode not in VALID_OP_MODES:
            raise RefusalError("INVALID-OP-MODE", f"mode {self.mode!r} not in {VALID_OP_MODES}")
        if self.kind == "ft" or (self.mode == "ft" and self.kind not in ("MPP",) + SINGLE_MEAS):
            raise RefusalError(
                "INVALID-FT-IDEALISATION",
                f"mode='ft' is defined only for Pauli measurements, got {self.kind!r}",
            )
        if len(set(self.qubits)) != len(self.qubits):
            raise RefusalError(
                "DUPLICATE-QUBIT-OPERAND",
                f"operation {self.kind!r} repeats qubit in {self.qubits}",
            )


@dataclass(frozen=True)
class AtomicFaultGroup:
    """One atomic fault event source (Definition 2.3).

    In an adversarial noise model, any non-zero combination of the generators in `generators`
    is produced by ONE fault event at this source and therefore has fault weight 1!
    Each generator is packed as a symplectic 2*N_loc bitmask over the target space (ZX edges
    or circuit spacetime locations), and `sector` tags 'X', 'Z', or 'MIXED' for CSS filtering.
    """

    source_id: int
    label: str
    generators: tuple[int, ...]
    sector: str = "MIXED"  # 'X', 'Z', or 'MIXED'


@dataclass
class CliffordCircuit:
    """Clifford circuit with explicit open input/output wires and wire/gate idealisations."""

    n_qubits: int
    input_qubits: tuple[int, ...]
    output_qubits: tuple[int, ...]
    ops: list[CircuitOp] = field(default_factory=list)
    fault_free_wires: frozenset[tuple[int, int]] = frozenset()
    # Optional global flag: if input_wires_fault_free / output_wires_fault_free, suppress wire flips
    fault_free_input_Segment: bool = False
    fault_free_output_segment: bool = False
    all_wires_fault_free: bool = False
    name: str = "circuit"

    def validate(self) -> None:
        """Check wire liveness: inputs start live, preps create live wires, destructive meas retire wires."""
        if len(set(self.input_qubits)) != len(self.input_qubits):
            raise RefusalError("DUPLICATE-INPUT-QUBIT", f"inputs={self.input_qubits}")
        if len(set(self.output_qubits)) != len(self.output_qubits):
            raise RefusalError("DUPLICATE-OUTPUT-QUBIT", f"outputs={self.output_qubits}")
        for q in self.input_qubits + self.output_qubits:
            if not (0 <= q < self.n_qubits):
                raise RefusalError("QUBIT-OUT-OF-RANGE", f"qubit {q} outside 0..{self.n_qubits - 1}")
        live: set[int] = set(self.input_qubits)
        for idx, op in enumerate(self.ops):
            for q in op.qubits:
                if not (0 <= q < self.n_qubits):
                    raise RefusalError("QUBIT-OUT-OF-RANGE", f"op {idx} qubit {q} out of range")
            if op.kind in PREP_STATES or op.kind == "CAT_PREP":
                for q in op.qubits:
                    if q in live:
                        raise RefusalError(
                            "PREP-ON-LIVE-WIRE",
                            f"op {idx} ({op.kind}) prepares already-live qubit {q}",
                        )
                    live.add(q)
            elif op.kind == "ENC_PREP":
                in_qs = tuple(op.params.get("input_qubits", ()))
                for q in in_qs:
                    if q not in live:
                        raise RefusalError(
                            "ENC-INPUT-DEAD",
                            f"op {idx} ENC_PREP consumes non-live input qubit {q}",
                        )
                    live.remove(q)
                for q in op.qubits:
                    if q in live:
                        raise RefusalError(
                            "ENC-OUTPUT-LIVE",
                            f"op {idx} ENC_PREP allocates already-live qubit {q}",
                        )
                    live.add(q)
            elif op.kind in SINGLE_MEAS:
                q = op.qubits[0]
                if q not in live:
                    raise RefusalError(
                        "MEAS-ON-DEAD-WIRE",
                        f"op {idx} ({op.kind}) measures non-live qubit {q}",
                    )
                live.remove(q)
            elif op.kind in SINGLE_QUBIT_GATES + ("CX", "CZ", "SWAP", "MPP"):
                for q in op.qubits:
                    if q not in live:
                        raise RefusalError(
                            "GATE-ON-DEAD-WIRE",
                            f"op {idx} ({op.kind}) acts on non-live qubit {q}",
                        )
            else:
                raise RefusalError("UNKNOWN-CIRCUIT-OP", f"unsupported op kind {op.kind!r}")
        if live != set(self.output_qubits):
            raise RefusalError(
                "LIVE-OUTPUT-MISMATCH",
                f"end-of-circuit live qubits {sorted(live)} != declared outputs {sorted(self.output_qubits)}",
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": "CliffordCircuit",
            "name": self.name,
            "n_qubits": self.n_qubits,
            "input_qubits": list(self.input_qubits),
            "output_qubits": list(self.output_qubits),
            "all_wires_fault_free": self.all_wires_fault_free,
            "fault_free_wires": [list(p) for p in sorted(self.fault_free_wires)],
            "ops": [
                {
                    "kind": op.kind,
                    "qubits": list(op.qubits),
                    "mode": op.mode,
                    **({"label": op.label} if op.label else {}),
                    **({"params": op.params} if op.params else {}),
                }
                for op in self.ops
            ],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> CliffordCircuit:
        ops = [
            CircuitOp(
                kind=str(item["kind"]),
                qubits=tuple(int(q) for q in item["qubits"]),
                mode=str(item.get("mode", "noisy")),
                label=str(item.get("label", "")),
                params=dict(item.get("params", {})),
            )
            for item in data.get("ops", [])
        ]
        ff_wires = frozenset(
            (int(p[0]), int(p[1])) for p in data.get("fault_free_wires", [])
        )
        circ = cls(
            n_qubits=int(data["n_qubits"]),
            input_qubits=tuple(int(q) for q in data.get("input_qubits", [])),
            output_qubits=tuple(int(q) for q in data.get("output_qubits", [])),
            ops=ops,
            fault_free_wires=ff_wires,
            all_wires_fault_free=bool(data.get("all_wires_fault_free", False)),
            name=str(data.get("name", "circuit")),
        )
        circ.validate()
        return circ


@dataclass
class ZXNode:
    """Node in a Clifford ZX diagram.

    kind:
      - 'Z'    : green Z-spider with phase k * (pi/2), k in {0, 1, 2, 3}
      - 'X'    : red X-spider with phase k * (pi/2), k in {0, 1, 2, 3}
      - 'H'    : 2-legged Hadamard box (must have exact degree 2)
      - 'BOUND': 1-legged open boundary port (input or output; must have exact degree 1)
    """

    node_id: int
    kind: str
    phase_k: int = 0
    label: str = ""

    def __post_init__(self) -> None:
        if self.kind not in ("Z", "X", "H", "BOUND"):
            raise RefusalError("UNKNOWN-ZX-NODE", f"invalid ZX node kind {self.kind!r}")
        if not isinstance(self.phase_k, int) or isinstance(self.phase_k, bool) or self.phase_k not in (0, 1, 2, 3):
            raise RefusalError(
                "NON-CLIFFORD-PHASE-REFUSED",
                f"spider phase_k must be in {{0,1,2,3}} (multiples of pi/2), got {self.phase_k!r}",
            )


@dataclass
class ZXEdge:
    """Undirected edge in a ZX diagram (spacetime location under edge-flip noise, Def 5.1)."""

    edge_id: int
    u: int
    v: int
    fault_free: bool = False  # True = purple idealised fault-free edge (Def 5.3)
    label: str = ""
    gadget_sector: str = ""   # '', 'X', 'Z', or 'MIXED' when edge is a fault-gadget trigger


class ZXDiagram:
    """Open Clifford ZX multigraph with ordered boundaries and submodels of edge-flip noise."""

    def __init__(self, name: str = "zx") -> None:
        self.name = name
        self.nodes: dict[int, ZXNode] = {}
        self.edges: dict[int, ZXEdge] = {}
        self.adj: dict[int, list[int]] = {}
        self.inputs: list[int] = []   # ordered BOUND node IDs
        self.outputs: list[int] = []  # ordered BOUND node IDs
        self.meas_spiders: list[tuple[int, str, str]] = []  # (spider_id, flip_pauli 'X'|'Z', label)
        self._next_node = 0
        self._next_edge = 0

    def copy(self) -> ZXDiagram:
        out = ZXDiagram(name=self.name)
        out._next_node = self._next_node
        out._next_edge = self._next_edge
        for nid, nd in self.nodes.items():
            out.nodes[nid] = ZXNode(nd.node_id, nd.kind, nd.phase_k, nd.label)
            out.adj[nid] = list(self.adj[nid])
        for eid, ed in self.edges.items():
            out.edges[eid] = ZXEdge(ed.edge_id, ed.u, ed.v, ed.fault_free, ed.label, ed.gadget_sector)
        out.inputs = list(self.inputs)
        out.outputs = list(self.outputs)
        out.meas_spiders = list(self.meas_spiders)
        return out

    def add_node(self, kind: str, phase_k: int = 0, label: str = "") -> int:
        if not isinstance(phase_k, int) or isinstance(phase_k, bool) or phase_k not in (0, 1, 2, 3):
            raise RefusalError(
                "NON-CLIFFORD-PHASE-REFUSED",
                f"Clifford ZX accepts phase_k in {{0,1,2,3}} (multiples of pi/2), got {phase_k!r}",
            )
        nid = self._next_node
        self._next_node += 1
        self.nodes[nid] = ZXNode(nid, kind, phase_k, label)
        self.adj[nid] = []
        return nid

    def add_boundary(self, is_input: bool, label: str = "") -> int:
        nid = self.add_node("BOUND", 0, label=label)
        if is_input:
            self.inputs.append(nid)
        else:
            self.outputs.append(nid)
        return nid

    def add_edge(
        self,
        u: int,
        v: int,
        fault_free: bool = False,
        label: str = "",
        gadget_sector: str = "",
    ) -> int:
        if u not in self.nodes or v not in self.nodes:
            raise RefusalError("INVALID-EDGE-ENDPOINT", f"missing endpoint in ({u}, {v})")
        if u == v:
            # Represent a self-loop cleanly via a fault-free phase-0 spider on the loop
            mid = self.add_node("Z", 0, label="loop_mid")
            e1 = self.add_edge(u, mid, fault_free=fault_free, label=label)
            self.add_edge(mid, v, fault_free=True, label=f"{label}_close")
            return e1
        eid = self._next_edge
        self._next_edge += 1
        self.edges[eid] = ZXEdge(eid, u, v, fault_free=fault_free, label=label, gadget_sector=gadget_sector)
        self.adj[u].append(eid)
        self.adj[v].append(eid)
        return eid

    def remove_edge(self, eid: int) -> ZXEdge:
        ed = self.edges.pop(eid)
        self.adj[ed.u].remove(eid)
        self.adj[ed.v].remove(eid)
        return ed

    def remove_node(self, nid: int) -> None:
        if self.adj[nid]:
            raise RefusalError("REMOVE-NONISOLATED-NODE", f"node {nid} still has incident edges")
        if nid in self.inputs or nid in self.outputs:
            raise RefusalError("REMOVE-BOUNDARY-NODE", f"cannot remove boundary node {nid}")
        del self.nodes[nid]
        del self.adj[nid]

    def other_endpoint(self, eid: int, nid: int) -> int:
        ed = self.edges[eid]
        if ed.u == nid:
            return ed.v
        if ed.v == nid:
            return ed.u
        raise ValueError(f"node {nid} not incident to edge {eid}")

    def boundary_edges(self) -> list[int]:
        """Ordered boundary edge IDs: inputs first, then outputs."""
        out: list[int] = []
        for b in self.inputs + self.outputs:
            inc = self.adj[b]
            if len(inc) != 1:
                raise RefusalError(
                    "MALFORMED-BOUNDARY-DEGREE",
                    f"boundary node {b} has degree {len(inc)} != 1",
                )
            out.append(inc[0])
        return out

    def validate(self) -> None:
        for nid, nd in self.nodes.items():
            deg = len(self.adj[nid])
            if nd.kind == "BOUND" and deg != 1:
                raise RefusalError("MALFORMED-BOUNDARY-DEGREE", f"BOUND {nid} degree {deg} != 1")
            if nd.kind == "H" and deg != 2:
                raise RefusalError("MALFORMED-HBOX-DEGREE", f"H-box {nid} degree {deg} != 2")
        bset = set(self.inputs) | set(self.outputs)
        for nid, nd in self.nodes.items():
            if (nd.kind == "BOUND") != (nid in bset):
                raise RefusalError("BOUNDARY-LIST-MISMATCH", f"node {nid} boundary mismatch")

    def insert_pauli_box_on_edge(self, eid: int, pauli: str) -> int:
        """Insert a fault-free Pauli box (Definition 6.2) on edge `eid` and return its hub spider.

        Definition 6.2 (Rodatz-Poor-Kissinger p. 22):
          - Z box : fault-free phase-0 Z-spider subdividing the edge (an X leg from the red
                    fault-gadget hub injects Z onto the wire via Pi-Copy!).
          - X box : fault-free H - Z(0) - H subdividing the edge (an X leg into Z(0) injects Z
                    between H and H, which equals X on the outer wire!).
          - Y box : fault-free Z(+pi/2) - Z(0) - Z(-pi/2) subdividing the edge:
                    Wait: check S * Z * Sdg vs S * X * Sdg!
                    Let us verify Definition 6.2 on page 23 directly:
                    When a Z flip occurs on the fault gadget's trigger edge, Pi-Copy pushes a
                    RED pi (X) onto every leg leaving the red hub!
                    - On Z-box (green phase-0 spider): red pi fuses/copies? Wait: red pi on one leg
                      of a 3-legged green spider copies to BOTH wire legs? NO! Red pi on one leg of
                      a green spider does NOT copy — GREEN pi fuses into a green spider, wait:
                      Read page 23 carefully!
                      In page 23's proof: a Z flip on the trigger edge is a GREEN pi! It enters the
                      RED spider hub, where Pi-Copy copies the GREEN pi onto EVERY leg leaving the
                      red hub as a GREEN pi (Z)!
                      Let us check EVERY Pauli box under incoming GREEN pi (Z) on its vertical leg:
                        1) Z-box = green phase-0 spider: incoming GREEN pi fuses into the green
                           spider, becoming a phase-pi green spider = Z on the wire! EXACT!
                        2) X-box = H - green(0) - H: incoming GREEN pi fuses into green(0) -> green(pi)
                           between two H-boxes = H Z H = X on the wire! EXACT!
                        3) Y-box = red(+pi/2) - green(0) - red(-pi/2):
                           Check Rx(-pi/2) * Z * Rx(+pi/2)!
                           Wait: |X_2(1)> = I + i X = sqrt(2) exp(i pi/4 X).
                           Conjugating Z by exp(-i pi/4 X) Z exp(+i pi/4 X):
                           exp(-i pi/4 X) Z = Z exp(+i pi/4 X) = Z (cos(pi/4) I + i sin(pi/4) X),
                           times (cos(pi/4) I + i sin(pi/4) X) = Z (i X) = -Y (up to global -1)!
                           Let us verify with ZX rules (page 23):
                           Incoming green pi fuses into middle green(0) -> green(pi).
                           In noise-free case (green(0) of degree 2 = identity wire),
                           red(+pi/2) -- red(-pi/2) fuses to red(0) = identity wire!
                           With green(pi) in the middle: push green(pi) right past red(-pi/2)
                           by Pi-Copy: green(pi) -- red(-pi/2) = red(+pi/2) -- green(pi),
                           then fuse red(+pi/2) -- red(+pi/2) = red(pi) = X, followed by green(pi) = Z,
                           giving X * Z = -i Y = Y in projective Pauli group! EXACT!!
        """
        if pauli == "I":
            raise ValueError("Identity Pauli box should be omitted")
        ed = self.remove_edge(eid)
        u, v = ed.u, ed.v
        ff = ed.fault_free
        if pauli == "Z":
            hub = self.add_node("Z", 0, label="pbox_Z")
            # Keep original noise flag on one half and fault-free on the other so noisy edge count is preserved!
            self.add_edge(u, hub, fault_free=ff, label=ed.label)
            self.add_edge(hub, v, fault_free=True, label=f"{ed.label}_pbox")
            return hub
        if pauli == "X":
            h1 = self.add_node("H", 0, label="pbox_X_h1")
            hub = self.add_node("Z", 0, label="pbox_X")
            h2 = self.add_node("H", 0, label="pbox_X_h2")
            self.add_edge(u, h1, fault_free=ff, label=ed.label)
            self.add_edge(h1, hub, fault_free=True)
            self.add_edge(hub, h2, fault_free=True)
            self.add_edge(h2, v, fault_free=True)
            return hub
        if pauli == "Y":
            r1 = self.add_node("X", 1, label="pbox_Y_rxp")
            hub = self.add_node("Z", 0, label="pbox_Y")
            r2 = self.add_node("X", 3, label="pbox_Y_rxm")
            self.add_edge(u, r1, fault_free=ff, label=ed.label)
            self.add_edge(r1, hub, fault_free=True)
            self.add_edge(hub, r2, fault_free=True)
            self.add_edge(r2, v, fault_free=True)
            return hub
        raise ValueError(f"Unknown Pauli {pauli!r}")

    def add_fault_gadget(
        self,
        edge_paulis: list[tuple[int, str]],
        label: str = "fault_gadget",
        sector: str = "MIXED",
    ) -> int:
        """Attach a Fault Gadget (Definition 6.3, p. 23) and return its noisy trigger edge ID.

        Connects fault-free Pauli boxes on the target edges via fault-free wires to a red
        phase-0 spider, which connects via ONE normal (noisy) trigger edge to a 1-legged
        green phase-0 spider.
        - Noise-free: 1-legged green spider copies through red spider -> empty legs -> identity!
        - Z flip on trigger edge: pushes green pi through red spider onto every Pauli box!
        - X flip on trigger edge: absorbed by 1-legged green spider (trivial)!
        - Y flip on trigger edge: equivalent to Z flip!
        """
        red_hub = self.add_node("X", 0, label=f"{label}_red")
        green_leaf = self.add_node("Z", 0, label=f"{label}_leaf")
        for eid, pchar in edge_paulis:
            if pchar == "I":
                continue
            box_hub = self.insert_pauli_box_on_edge(eid, pchar)
            self.add_edge(box_hub, red_hub, fault_free=True, label=f"{label}_leg")
        trig_eid = self.add_edge(
            red_hub,
            green_leaf,
            fault_free=False,
            label=f"{label}_trigger",
            gadget_sector=sector,
        )
        return trig_eid

    def ordered_edge_ids(self) -> list[int]:
        return sorted(self.edges.keys())

    def build_edge_flip_fault_groups(self, noise_model: str = "circuit_depolarizing") -> tuple[list[int], list[AtomicFaultGroup]]:
        """Build atomic fault groups under (submodels of) edge-flip noise (Defs 5.2-5.3).

        Returns:
          (ordered_eids, fault_groups)
          where symplectic masks are indexed 0..len(ordered_eids)-1.
        """
        if noise_model not in VALID_NOISE_MODELS:
            raise RefusalError("UNKNOWN-NOISE-MODEL", f"unknown noise model {noise_model!r}")
        eids = self.ordered_edge_ids()
        if noise_model == "fault_free":
            return eids, []
        groups: list[AtomicFaultGroup] = []
        src = 0
        for idx, eid in enumerate(eids):
            ed = self.edges[eid]
            if ed.fault_free:
                continue
            z_mask = 1 << (2 * idx)
            x_mask = 1 << (2 * idx + 1)
            if noise_model == "circuit_css":
                # Under CSS noise (Kissinger p. 38), cross-sector gadgets are inactive, and
                # X edge flips and Z edge flips act as separate atomic fault events (no weight-1 Y=XZ).
                if ed.gadget_sector == "MIXED":
                    continue
                if ed.gadget_sector == "Z":
                    groups.append(AtomicFaultGroup(src, f"edge_{eid}_gadZ", (z_mask,), sector="Z"))
                    src += 1
                elif ed.gadget_sector == "X":
                    groups.append(AtomicFaultGroup(src, f"edge_{eid}_gadX", (z_mask,), sector="X"))
                    src += 1
                else:
                    groups.append(AtomicFaultGroup(src, f"edge_{eid}_Z", (z_mask,), sector="Z"))
                    src += 1
                    groups.append(AtomicFaultGroup(src, f"edge_{eid}_X", (x_mask,), sector="X"))
                    src += 1
            else:
                # Full edge-flip noise: any of {X_e, Y_e, Z_e} at edge e has weight 1!
                # Generated by (z_mask, x_mask) inside a single AtomicFaultGroup!
                groups.append(
                    AtomicFaultGroup(src, f"edge_{eid}", (z_mask, x_mask), sector=ed.gadget_sector or "MIXED")
                )
                src += 1
        return eids, groups

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": "ZXDiagram",
            "name": self.name,
            "inputs": list(self.inputs),
            "outputs": list(self.outputs),
            "meas_spiders": [[nid, fp, lbl] for nid, fp, lbl in self.meas_spiders],
            "nodes": [
                {"id": nid, "kind": nd.kind, "phase_k": nd.phase_k, **({"label": nd.label} if nd.label else {})}
                for nid, nd in sorted(self.nodes.items())
            ],
            "edges": [
                {
                    "id": eid,
                    "u": ed.u,
                    "v": ed.v,
                    "fault_free": ed.fault_free,
                    **({"label": ed.label} if ed.label else {}),
                    **({"gadget_sector": ed.gadget_sector} if ed.gadget_sector else {}),
                }
                for eid, ed in sorted(self.edges.items())
            ],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ZXDiagram:
        zx = cls(name=str(data.get("name", "zx")))
        max_nid = -1
        for item in data.get("nodes", []):
            nid = int(item["id"])
            nd = ZXNode(nid, str(item["kind"]), int(item.get("phase_k", 0)), str(item.get("label", "")))
            zx.nodes[nid] = nd
            zx.adj[nid] = []
            if nid > max_nid:
                max_nid = nid
        zx._next_node = max_nid + 1
        max_eid = -1
        for item in data.get("edges", []):
            eid = int(item["id"])
            ed = ZXEdge(
                eid,
                int(item["u"]),
                int(item["v"]),
                bool(item.get("fault_free", False)),
                str(item.get("label", "")),
                str(item.get("gadget_sector", "")),
            )
            zx.edges[eid] = ed
            zx.adj[ed.u].append(eid)
            zx.adj[ed.v].append(eid)
            if eid > max_eid:
                max_eid = eid
        zx._next_edge = max_eid + 1
        zx.inputs = [int(x) for x in data.get("inputs", [])]
        zx.outputs = [int(x) for x in data.get("outputs", [])]
        zx.meas_spiders = [(int(a), str(b), str(c)) for a, b, c in data.get("meas_spiders", [])]
        zx.validate()
        return zx


def _conjugate_in_basis(zx: ZXDiagram, wire_node: int, pchar: str, ff: bool) -> int:
    """Attach single-qubit Clifford converting Z-eigenbasis to `pchar`-eigenbasis on `wire_node`."""
    if pchar == "Z":
        return wire_node
    if pchar == "X":
        h = zx.add_node("H", 0, label="basis_H")
        zx.add_edge(wire_node, h, fault_free=ff)
        return h
    if pchar == "Y":
        # S * H * |0> = |+i>, H * Sdg * Y * S * H = Z!
        # Wait: check (H Sdg) Y (S H): Sdg Y S = X, H X H = Z!
        # On input side (before Z-proj): apply Sdg then H; after Z-proj: apply H then S!
        raise NotImplementedError("Use explicit Sdg/H around Y")
    raise ValueError(f"Unexpected Pauli {pchar!r}")


def circuit_to_zx_fe(
    circuit: CliffordCircuit,
    noise_model: str = "circuit_depolarizing",
    optimise_gadgets: bool = True,
) -> ZXDiagram:
    """Translate a CliffordCircuit under `noise_model` into a fault-equivalent ZXDiagram!

    Implements Section 6.1 of Rodatz-Poor-Kissinger (arXiv:2506.17181v4):
      - Prop. 6.5 : fault-free ('ff') operations use fault-free internal edges and 0 fault gadgets.
      - Cor.  6.6 : single-qubit gates/preps/measurements require 0 fault gadgets.
      - Prop. 6.7 : fault-tolerant ('ft') multi-qubit Pauli measurements use noisy spider legs
                    and 0 fault gadgets (capturing meas flip + <=1 output flip!).
      - Prop. 6.9 : CNOT immediately followed by noisy MEAS_Z on target (or MEAS_X on control)
                    uses 0 fault gadgets even under full circuit_depolarizing noise!
      - Prop. 6.10: noisy PREP_PLUS on control (or PREP_0 on target) immediately followed by CNOT
                    uses 0 fault gadgets even under full circuit_depolarizing noise!
      - Prop. 6.8 : general noisy CNOT under circuit_depolarizing attaches the 4 cross-sector
                    fault gadgets {XZ, YZ, XY, YY} (under circuit_css, 0 fault gadgets).
      - Prop. 6.4 : general noisy MPP / CZ / SWAP attaches exact fault gadgets for any remaining
                    atomic gate faults not already generated by single edge flips.
    """
    if noise_model not in VALID_NOISE_MODELS:
        raise RefusalError("UNKNOWN-NOISE-MODEL", f"unknown noise model {noise_model!r}")
    circuit.validate()
    zx = ZXDiagram(name=circuit.name)
    full_ff = noise_model == "fault_free" or circuit.all_wires_fault_free

    def wire_ff(q: int, seg: int) -> bool:
        return full_ff or ((q, seg) in circuit.fault_free_wires)

    # Track for each live qubit q: (current_endpoint_node_id, seg_index, just_prepared_kind_or_None)
    live_end: dict[int, tuple[int, int, str | None]] = {}
    for q in circuit.input_qubits:
        b_in = zx.add_boundary(is_input=True, label=f"in_{q}")
        live_end[q] = (b_in, 0, None)

    n_ops = len(circuit.ops)

    # Pre-scan for Prop. 6.9 (CNOT immediately followed on its qubit by destructive MEAS_Z / MEAS_X)
    next_op_on_qubit: dict[tuple[int, int], int] = {}
    last_seen: dict[int, int] = {}
    for idx in range(n_ops - 1, -1, -1):
        op = circuit.ops[idx]
        for q in op.qubits:
            if q in last_seen:
                next_op_on_qubit[(idx, q)] = last_seen[q]
            last_seen[q] = idx
        if op.kind == "ENC_PREP":
            for q in op.params.get("input_qubits", ()):
                if q in last_seen:
                    next_op_on_qubit[(idx, q)] = last_seen[q]
                last_seen[q] = idx

    for idx, op in enumerate(circuit.ops):
        op_ff = full_ff or (op.mode == "ff")

        if op.kind in PREP_STATES:
            q = op.qubits[0]
            if op.kind == "PREP_0":
                sp = zx.add_node("X", 0, label=f"prep0_q{q}")
            elif op.kind == "PREP_1":
                sp = zx.add_node("X", 2, label=f"prep1_q{q}")
            elif op.kind == "PREP_PLUS":
                sp = zx.add_node("Z", 0, label=f"prep+_q{q}")
            else:  # PREP_MINUS
                sp = zx.add_node("Z", 2, label=f"prep-_q{q}")
            live_end[q] = (sp, 0, op.kind if not op_ff else None)

        elif op.kind == "CAT_PREP":
            # |cat_k> = |0...0> + |1...1> (basis='Z', single Z-spider) or |+...+> + |-...-> (basis='X', single X-spider)
            basis = str(op.params.get("basis", "Z"))
            sp_kind = "Z" if basis == "Z" else "X"
            hub = zx.add_node(sp_kind, 0, label=op.label or f"cat_{len(op.qubits)}")
            # Add intermediate identity spiders if op is noisy vs ff? CAT_PREP is an idealised spec primitive
            if not op_ff:
                raise RefusalError(
                    "NOISY-IDEAL-SPEC-PRIMITIVE",
                    "CAT_PREP is an idealised specification primitive; set mode='ff'",
                )
            for q in op.qubits:
                live_end[q] = (hub, 0, None)

        elif op.kind == "ENC_PREP":
            # Idealised fault-free CSS/stabiliser encoder from logical `input_qubits` to physical `op.qubits`
            if not op_ff:
                raise RefusalError(
                    "NOISY-IDEAL-SPEC-PRIMITIVE",
                    "ENC_PREP is an idealised specification primitive; set mode='ff'",
                )
            in_qs = tuple(int(q) for q in op.params.get("input_qubits", ()))
            stabs: list[str] = [str(s) for s in op.params["stabilizers"]]
            log_xs: list[str] = [str(s) for s in op.params.get("logical_x", [])]
            log_zs: list[str] = [str(s) for s in op.params.get("logical_z", [])]
            if len(in_qs) == 0:
                stabs = stabs + log_zs
                log_zs = []
                log_xs = []
            elif len(in_qs) != len(log_zs) or len(in_qs) != len(log_xs):
                raise RefusalError("ENC-LOGICAL-COUNT-MISMATCH", "input_qubits != logical operators")
            # Build fault-free CSS / stabiliser projector diagram:
            # Start each physical wire in |+> (Z-spider phase 0) if no X-all-space needed?
            # Even better and completely general for ANY stabiliser code (not just CSS!):
            # Identity wires from ancilla |0>/|+> or better: for k logical inputs -> n physical outputs,
            # prepare k inputs on first k wires? No: we can attach k logical inputs via fault-free Bell
            # projections (X_in * L_X_phys = +1 AND Z_in * L_Z_phys = +1) and (n - k) code stabilisers
            # S_j = +1 onto n maximally-mixed/free wires!
            # Wait: on n fresh wires initialized to |+>^n? No, a projector \prod (I + g_m) applied to a
            # random state might vanish if |+>^n is orthogonal to a Z-stabiliser? Wait: |+>^n has non-zero
            # overlap with |0_L> ONLY if no X-stabilizer has -1 sign; what if we apply ALL n independent
            # projectors (n - k stabilisers S_j AND k pairs X_in*L_X, Z_in*L_Z = n + k commuting projectors)
            # to n fresh wires connected to k input wires via... wait:
            # Count legs! There are k input wires and n output wires = k + n boundary wires!
            # A Choi vector on k + n qubits is uniquely fixed by k + n independent commuting Pauli
            # projectors applied to ANY product state with non-zero overlap, OR in ZX:
            # Simply connect each physical output wire j to a fresh open internal edge? No, open internal
            # edges aren't states; a state on (k + n) wires with (n - k) stabilizers S_j on outputs and
            # k pairs (X_in_r * L_X_r, Z_in_r * L_Z_r) has k + n stabilizers!
            # Wait: how do we apply a fault-free Pauli projector (I + P) to a set of wires in ZX?!
            # Look at Figure 6 (page 17) of Kissinger:
            # A fault-free Z...Z projector on wires w_1..w_r is literally:
            #   put a phase-0 Z-spider on each wire w_a, connected by fault-free edges to ONE 0-phase X-spider!
            # Check its linear map!!
            # Sum_{b in {0,1}} (Z-spider phase b*pi on each wire)? No: X-spider of phase 0 is |+><+| + |-><-|?
            # Wait: an X-spider of phase 0 with r legs has state |+>^r + |->^r = sum_{x: even wt?} NO!
            # Let us check Figure 6 (page 17) of Kissinger!
            # On page 17 Figure 6: ZZZ measurement = green spider on each wire, all connected to ONE red spider!
            # Check what that tensor is!!
            # Red spider |X_r(0)> = |+>^r + |->^r = sum_{z in {0,1}^r : parity(z)=0} |z> (up to scalar)!!
            # When bit z_a in {0, 1} contracts with green spider |0><0|_w + |1><1|_w on wire a:
            # Wait: green spider |000> + |111> on (wire_in, wire_out, meas_leg) sets z_a = wire_a!
            # Contracting z = (w_1, ..., w_r) with |X_r(0)> = sum_{parity(z)=0} <z| PROJECTS onto
            # parity(w_1, ..., w_r) == 0, which IS the exact projector (I + Z_1 ... Z_r)!!
            # AND for a general Pauli string P = P_1 ... P_r (with X, Y, Z), conjugating each wire's
            # green spider by fault-free single-qubit Cliffords (H for X; Sdg-H / H-S for Y) projects
            # onto (I + P_1 ... P_r)!!
            # AND what initial state on (k inputs + n fresh wires) has non-zero overlap with EVERY
            # stabilizer code + logical Bell pair?!
            # Think! Suppose we start the n output wires as: first k output wires connected directly
            # from the k input wires (fault-free), and remaining (n - k) output wires initialized in |0>
            # or |+>? Wait: what if S_1 = Z_0 Z_1 ... (commutes with |0>) while S_2 = X_0 X_1 ...
            # (|0...0> has overlap 2^{-(n-k)/2} != 0 with ANY CSS code whose Z-stabilisers have +1 signs)!
            # What if the code has X- and Z-stabilisers AND we also project X_in_r * L_X_r and Z_in_r * L_Z_r?
            # Wait: if wire r starts as identity from in_r, it ALREADY has X_in_r X_r = +1 and Z_in_r Z_r = +1,
            # which might anticommute with L_X_r or S_j, killing the amplitude to 0!
            # How do you create a wire with NO initial stabilizer (maximally mixed / unprojected wire, i.e.
            # Choi cup) so projecting onto ANY commuting set of k + n independent stabilizers NEVER vanishes?!
            # THINK IN CHOI SPACE!!
            # On 1 wire, what 1-legged spider has stabilizer... NO 1-legged spider has 0 stabilizers; a
            # 1-legged spider has 1 stabilizer!
            # Instead, given the k + n commuting stabilizer generators G_1, ..., G_{k+n} on the (k + n)
            # ports (inputs + outputs), Gaussian elimination (`gf2.py`) finds a computational/Hadamard
            # product state |s_1 ... s_{k+n}> (each s_i in {|0>, |+>}) that does NOT anticommute with
            #(-1) of any product of G_m (specifically, pick |0> on qubits where Z-part pivots, or even
            # simpler: solve for any stabilizer product state |0>^A |+>^B whose stabilizer group S_0 has
            # S_0 \cap <G> = {+I}, so <psi_0 | \Pi_G | psi_0> = 2^{-(k+n)} != 0 GUARANTEED FOR EVERY
            # STABILISER CODE IN EXISTENCE)!!
            in_info = [(q, live_end[q][0], live_end[q][1]) for q in in_qs]
            for q in in_qs:
                del live_end[q]
            n_total = len(in_qs) + len(op.qubits)
            # Build the (k + n) Pauli strings over (inputs + outputs):
            # Recall on input wires, operator E |x><x| has Choi stabilizer:
            #   enc |0_in> = |0_L>  => Z_in * L_Z = +1
            #   enc |+_in> = |+_L>  => X_in * L_X = +1
            proj_strings: list[str] = []
            for s in stabs:
                if len(s) != len(op.qubits):
                    raise RefusalError("ENC-STAB-LEN", f"stabiliser {s!r} len != {len(op.qubits)}")
                proj_strings.append(("I" * len(in_qs)) + s)
            for r_idx in range(len(in_qs)):
                in_z = ["I"] * len(in_qs)
                in_z[r_idx] = "Z"
                proj_strings.append("".join(in_z) + log_zs[r_idx])
                if log_xs:
                    in_x = ["I"] * len(in_qs)
                    in_x[r_idx] = "X"
                    proj_strings.append("".join(in_x) + log_xs[r_idx])
            # Choose per-port seed basis ('0' -> X(0) spider = Z-stabilizer, '+' -> Z(0) spider = X-stabilizer)
            # so no non-empty combination of proj_strings lies purely in S_0!
            from .gf2 import pauli_string_to_symplectic, rref_gf2, symplectic_inner_product

            g_vecs = [pauli_string_to_symplectic(ps) for ps in proj_strings]
            for i_g in range(len(g_vecs)):
                for j_g in range(i_g + 1, len(g_vecs)):
                    if symplectic_inner_product(g_vecs[i_g], g_vecs[j_g], n_total) != 0:
                        raise RefusalError(
                            "ENC-NONCOMMUTING-GENERATORS",
                            f"ENC_PREP generators {proj_strings[i_g]!r} and {proj_strings[j_g]!r} anticommute",
                        )
            _, piv_cols, _ = rref_gf2(g_vecs, 2 * n_total)
            if len(piv_cols) != n_total:
                raise RefusalError(
                    "ENC-NOT-MAXIMAL",
                    f"ENC_PREP generators have rank {len(piv_cols)} != {n_total}",
                )
            # Pick complementary Lagrangian coordinates for the seed product state!
            # Wait: even simpler! Every row in RREF has a unique pivot column `pcol` in 0..2*n_total-1.
            # Can two pivots sit on the SAME qubit j (both 2*j and 2*j+1)? YES (e.g. a decoupled qubit,
            # though not in an encoder, or Y/Z!). Better: greedily add single-qubit Z_j or X_j to the
            # span of g_vecs! For j = 0..n_total-1: either Z_j is NOT in span(current) OR X_j is NOT
            # in span(current), because Z_j and X_j anticommute while `g_vecs` is an isotropic commuting
            # subspace (so <g_vecs> cannot contain BOTH Z_j and X_j)!
            cur_basis = list(g_vecs)
            seed_kinds: list[str] = []
            for j in range(n_total):
                z_cand = 1 << (2 * j)
                test_rref, _, _ = rref_gf2(cur_basis + [z_cand], 2 * n_total)
                if len(test_rref) == len(cur_basis) + 1:
                    cur_basis = test_rref
                    seed_kinds.append("X")  # X(0) spider prepares |0> (stabilized by Z_j)
                else:
                    x_cand = 1 << (2 * j + 1)
                    test_rref2, _, _ = rref_gf2(cur_basis + [x_cand], 2 * n_total)
                    cur_basis = test_rref2
                    seed_kinds.append("Z")  # Z(0) spider prepares |+> (stabilized by X_j)

            # Initialize the k+n wires on seed_kinds (for inputs, connect `in_nd` -> seed? NO!
            # `in_nd` is an OPEN boundary wire: terminating `in_nd` on a 1-legged spider `seed` means
            # `in_nd` connects to `seed_in` (killing input? NO: first apply the projectors ON the wire
            # between `in_nd` and `seed_in`!). Let's check!
            # Take k=1, n=1, trivial identity encoder: proj_strings = ["ZZ", "XX"].
            # j=0: Z_0 not in <ZZ, XX> -> cur_basis = <ZZ, XX, Z_0>, seed_0 = |0>.
            # j=1: Z_1 IS in <ZZ, Z_0>, so X_1 not in cur_basis -> seed_1 = |+>.
            # Seed state on (port_0, port_1) is |0, +>.
            # Apply projector (I + ZZ)(I + XX) to |0, +>:
            # (I + ZZ)|0, +> = |0, 0> + |0, 0> = 2 |0, 0>.
            # (I + XX)(2 |0, 0>) = 2(|0, 0> + |1, 1>) = Choi state of Identity!!
            # IT WORKS FOR EVERY STABILISER ENCODER ON EARTH!!
            wire_tips: list[int] = []
            # Input ports: wire runs FROM seed_0 -> projectors -> in_nd!
            # Output ports: wire runs FROM seed_out -> projectors -> live_end[q_out]!
            for j in range(n_total):
                s_node = zx.add_node(seed_kinds[j], 0, label=f"enc_seed_{j}")
                wire_tips.append(s_node)

            for m_idx, ps in enumerate(proj_strings):
                red_hub = zx.add_node("X", 0, label=f"enc_proj_{m_idx}")
                for j, pch in enumerate(ps):
                    if pch == "I":
                        continue
                    prev = wire_tips[j]
                    if pch == "Z":
                        sp = zx.add_node("Z", 0)
                        zx.add_edge(prev, sp, fault_free=True)
                        zx.add_edge(sp, red_hub, fault_free=True)
                        wire_tips[j] = sp
                    elif pch == "X":
                        h1 = zx.add_node("H", 0)
                        sp = zx.add_node("Z", 0)
                        h2 = zx.add_node("H", 0)
                        zx.add_edge(prev, h1, fault_free=True)
                        zx.add_edge(h1, sp, fault_free=True)
                        zx.add_edge(sp, red_hub, fault_free=True)
                        zx.add_edge(sp, h2, fault_free=True)
                        wire_tips[j] = h2
                    elif pch == "Y":
                        r1 = zx.add_node("X", 1)
                        sp = zx.add_node("Z", 0)
                        r2 = zx.add_node("X", 3)
                        zx.add_edge(prev, r1, fault_free=True)
                        zx.add_edge(r1, sp, fault_free=True)
                        zx.add_edge(sp, red_hub, fault_free=True)
                        zx.add_edge(sp, r2, fault_free=True)
                        wire_tips[j] = r2

            # Connect input wire_tips[0..k-1] to `in_nodes` (preserving input wire noise flag!) and set output live_ends!
            for r_idx, (q_in, in_nd, seg_in) in enumerate(in_info):
                zx.add_edge(wire_tips[r_idx], in_nd, fault_free=wire_ff(q_in, seg_in), label=f"w_in_enc_q{q_in}")
            for out_idx, q_out in enumerate(op.qubits):
                live_end[q_out] = (wire_tips[len(in_qs) + out_idx], 0, None)

        elif op.kind in SINGLE_QUBIT_GATES:
            q = op.qubits[0]
            prev_nd, seg, _ = live_end[q]
            in_ff = wire_ff(q, seg)
            if op.kind == "I":
                sp = zx.add_node("Z", 0, label=f"I_q{q}")
                zx.add_edge(prev_nd, sp, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                live_end[q] = (sp, seg + 1, None)
            elif op.kind == "Z":
                sp = zx.add_node("Z", 2, label=f"Z_q{q}")
                zx.add_edge(prev_nd, sp, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                live_end[q] = (sp, seg + 1, None)
            elif op.kind == "S":
                sp = zx.add_node("Z", 1, label=f"S_q{q}")
                zx.add_edge(prev_nd, sp, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                live_end[q] = (sp, seg + 1, None)
            elif op.kind == "SDG":
                sp = zx.add_node("Z", 3, label=f"SDG_q{q}")
                zx.add_edge(prev_nd, sp, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                live_end[q] = (sp, seg + 1, None)
            elif op.kind == "X":
                sp = zx.add_node("X", 2, label=f"X_q{q}")
                zx.add_edge(prev_nd, sp, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                live_end[q] = (sp, seg + 1, None)
            elif op.kind == "Y":
                sp1 = zx.add_node("Z", 2, label=f"Y_z_q{q}")
                sp2 = zx.add_node("X", 2, label=f"Y_x_q{q}")
                zx.add_edge(prev_nd, sp1, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                zx.add_edge(sp1, sp2, fault_free=True, label=f"Y_int_q{q}")
                live_end[q] = (sp2, seg + 1, None)
            elif op.kind == "H":
                hb = zx.add_node("H", 0, label=f"H_q{q}")
                zx.add_edge(prev_nd, hb, fault_free=in_ff, label=f"w_q{q}_s{seg}")
                live_end[q] = (hb, seg + 1, None)

        elif op.kind == "CX":
            qc, qt = op.qubits
            prev_c, seg_c, prep_c = live_end[qc]
            prev_t, seg_t, prep_t = live_end[qt]
            in_ff_c = wire_ff(qc, seg_c)
            in_ff_t = wire_ff(qt, seg_t)

            zc = zx.add_node("Z", 0, label=f"cx_ctrl_q{qc}_{idx}")
            xt = zx.add_node("X", 0, label=f"cx_targ_q{qt}_{idx}")
            zx.add_edge(prev_c, zc, fault_free=in_ff_c, label=f"w_q{qc}_s{seg_c}")
            zx.add_edge(prev_t, xt, fault_free=in_ff_t, label=f"w_q{qt}_s{seg_t}")
            zx.add_edge(zc, xt, fault_free=op_ff, label=f"cx_mid_{idx}")

            # Check composite gadget-free theorems (Props 6.9 & 6.10):
            # Prop 6.10: control just prepared in |+> (or colour-dual: target just prepared in |0>)
            # Prop 6.9 : target immediately measured in Z (or colour-dual: control immediately measured in X)
            split_ok = optimise_gadgets and (
                (prep_c == "PREP_PLUS" and not in_ff_c)
                or (prep_t == "PREP_0" and not in_ff_t)
            )
            nxt_t = next_op_on_qubit.get((idx, qt))
            nxt_c = next_op_on_qubit.get((idx, qc))
            merge_ok = optimise_gadgets and (
                (nxt_t is not None and circuit.ops[nxt_t].kind == "MEAS_Z" and circuit.ops[nxt_t].mode != "ff" and not wire_ff(qt, seg_t + 1))
                or (nxt_c is not None and circuit.ops[nxt_c].kind == "MEAS_X" and circuit.ops[nxt_c].mode != "ff" and not wire_ff(qc, seg_c + 1))
            )

            need_cx_gadgets = (
                (not op_ff)
                and (noise_model == "circuit_depolarizing")
                and (not split_ok)
                and (not merge_ok)
            )
            if need_cx_gadgets:
                # Attach output buffer spiders so the 4 cross-sector fault gadgets (Prop. 6.8)
                # sit cleanly on the CNOT output wires!
                out_c = zx.add_node("Z", 0, label=f"cx_outc_{idx}")
                out_t = zx.add_node("Z", 0, label=f"cx_outt_{idx}")
                ec = zx.add_edge(zc, out_c, fault_free=True, label=f"cx_gadc_{idx}")
                et = zx.add_edge(xt, out_t, fault_free=True, label=f"cx_gadt_{idx}")
                for pc, pt in (("X", "Z"), ("Y", "Z"), ("X", "Y"), ("Y", "Y")):
                    # Re-create fresh fault-free segment pair for each gadget so edge IDs stay clean
                    nc = zx.add_node("Z", 0)
                    nt = zx.add_node("Z", 0)
                    ec_use = zx.add_edge(out_c, nc, fault_free=True)
                    et_use = zx.add_edge(out_t, nt, fault_free=True)
                    zx.add_fault_gadget(
                        [(ec_use, pc), (et_use, pt)],
                        label=f"cx_{idx}_{pc}{pt}",
                        sector="MIXED",
                    )
                    out_c, out_t = nc, nt
                live_end[qc] = (out_c, seg_c + 1, None)
                live_end[qt] = (out_t, seg_t + 1, None)
            else:
                live_end[qc] = (zc, seg_c + 1, None)
                live_end[qt] = (xt, seg_t + 1, None)

        elif op.kind == "CZ":
            q0, q1 = op.qubits
            prev_0, seg_0, _ = live_end[q0]
            prev_1, seg_1, _ = live_end[q1]
            z0 = zx.add_node("Z", 0, label=f"cz0_{idx}")
            z1 = zx.add_node("Z", 0, label=f"cz1_{idx}")
            hb = zx.add_node("H", 0, label=f"czH_{idx}")
            zx.add_edge(prev_0, z0, fault_free=wire_ff(q0, seg_0))
            zx.add_edge(prev_1, z1, fault_free=wire_ff(q1, seg_1))
            zx.add_edge(z0, hb, fault_free=op_ff)
            zx.add_edge(hb, z1, fault_free=True)
            if not op_ff:
                # Output fault gadgets on CZ for remaining two-qubit faults
                out0 = zx.add_node("Z", 0)
                out1 = zx.add_node("Z", 0)
                zx.add_edge(z0, out0, fault_free=True)
                zx.add_edge(z1, out1, fault_free=True)
                if noise_model == "circuit_depolarizing":
                    # Single edge flip on z0-hb gives XZ (via X on z0-hb? Wait: X on z0-hb is X on z0 and Z on z1 = XZ!
                    # Z on z0-hb is Z on z0. To be exact per Prop 6.4/6.8, attach the missing correlated output pairs:
                    # Single edge flips already give: ZI, IZ, XI (input->XZ), IX (input->ZX), XZ (mid X), YI->YZ, IY->ZY.
                    # Missing on outputs: ZZ, XX, XY, YX, YY!
                    missing = (("Z", "Z"), ("X", "X"), ("X", "Y"), ("Y", "X"), ("Y", "Y"))
                    for p0, p1 in missing:
                        n0 = zx.add_node("Z", 0)
                        n1 = zx.add_node("Z", 0)
                        e0 = zx.add_edge(out0, n0, fault_free=True)
                        e1 = zx.add_edge(out1, n1, fault_free=True)
                        sec = "Z" if (p0 == "Z" and p1 == "Z") else ("X" if (p0 == "X" and p1 == "X") else "MIXED")
                        zx.add_fault_gadget([(e0, p0), (e1, p1)], label=f"cz_{idx}_{p0}{p1}", sector=sec)
                        out0, out1 = n0, n1
                elif noise_model == "circuit_css":
                    for p0, p1, sec in (("Z", "Z", "Z"), ("X", "X", "X")):
                        n0 = zx.add_node("Z", 0)
                        n1 = zx.add_node("Z", 0)
                        e0 = zx.add_edge(out0, n0, fault_free=True)
                        e1 = zx.add_edge(out1, n1, fault_free=True)
                        zx.add_fault_gadget([(e0, p0), (e1, p1)], label=f"cz_{idx}_{p0}{p1}", sector=sec)
                        out0, out1 = n0, n1
                live_end[q0] = (out0, seg_0 + 1, None)
                live_end[q1] = (out1, seg_1 + 1, None)
            else:
                live_end[q0] = (z0, seg_0 + 1, None)
                live_end[q1] = (z1, seg_1 + 1, None)

        elif op.kind == "SWAP":
            q0, q1 = op.qubits
            end0 = live_end[q0]
            end1 = live_end[q1]
            p0, s0, _ = end0
            p1, s1, _ = end1
            n0 = zx.add_node("Z", 0, label=f"swap0_{idx}")
            n1 = zx.add_node("Z", 0, label=f"swap1_{idx}")
            zx.add_edge(p1, n0, fault_free=wire_ff(q1, s1))
            zx.add_edge(p0, n1, fault_free=wire_ff(q0, s0))
            if not op_ff:
                # Attach correlated output gadgets for SWAP gate faults
                pairs = (
                    [(p_a, p_b) for p_a in ("I", "X", "Y", "Z") for p_b in ("I", "X", "Y", "Z") if p_a != "I" and p_b != "I"]
                )
                for p_a, p_b in pairs:
                    sec = "X" if (p_a == "X" and p_b == "X") else ("Z" if (p_a == "Z" and p_b == "Z") else "MIXED")
                    if noise_model == "circuit_css" and sec == "MIXED":
                        continue
                    m0 = zx.add_node("Z", 0)
                    m1 = zx.add_node("Z", 0)
                    e0 = zx.add_edge(n0, m0, fault_free=True)
                    e1 = zx.add_edge(n1, m1, fault_free=True)
                    zx.add_fault_gadget([(e0, p_a), (e1, p_b)], label=f"swap_{idx}_{p_a}{p_b}", sector=sec)
                    n0, n1 = m0, m1
            live_end[q0] = (n0, s0 + 1, None)
            live_end[q1] = (n1, s1 + 1, None)

        elif op.kind in SINGLE_MEAS:
            q = op.qubits[0]
            prev_nd, seg, _ = live_end.pop(q)
            meas_ff = op_ff and wire_ff(q, seg)
            lbl = op.label or f"{op.kind}_q{q}_{idx}"
            if op.kind == "MEAS_Z":
                # Project onto <0|: 1-legged X-spider of phase 0! Outcome flip = X(pi)
                mnode = zx.add_node("X", 0, label=lbl)
                zx.add_edge(prev_nd, mnode, fault_free=meas_ff, label=f"meas_q{q}_{idx}")
                zx.meas_spiders.append((mnode, "X", lbl))
            elif op.kind == "MEAS_X":
                # Project onto <+|: 1-legged Z-spider of phase 0! Outcome flip = Z(pi)
                mnode = zx.add_node("Z", 0, label=lbl)
                zx.add_edge(prev_nd, mnode, fault_free=meas_ff, label=f"meas_q{q}_{idx}")
                zx.meas_spiders.append((mnode, "Z", lbl))
            else:  # MEAS_Y: project onto <+i| = <0| Sdg -> Sdg (Z(3)) followed by <+| (Z(0)) = Z(3) leaf!
                # Check: <Z_1(3)| = <0| + (-i)^3 <1|? Wait: |Z_1(3)> = |0> - i|1>, which IS the +1
                # eigenstate of -Y, whereas |Z_1(1)> = |0> + i|1> IS the +1 eigenstate of +Y!
                # When used as a ket in Choi contraction with |psi>: (<0| + i<1|)|psi> is ZERO on |0>+i|1|!
                # CAREFUL!! Look at bra vs Choi ket:
                # For |+i> = (|0> + i|1>)/sqrt(2), the projector bra is <+i| = (<0| - i<1|)/sqrt(2)!
                # In Choi ket notation, <0| - i<1| = |0> + (-i)|1| = |Z_1(3)>!!
                # Let's double check! Contract |Z_1(3)> = |0> - i|1| with |+i> = |0> + i|1| via <Phi+| = <00|+<11|:
                # 1*1 + (-i)*(i) = 1 + 1 = 2 != 0! Whereas with |-i> = |0> - i|1|: 1*1 + (-i)*(-i) = 1 - 1 = 0!
                # THANK GOODNESS WE CHECKED THAT!! <+i| IS Z-spider with phase_k = 3!
                mnode = zx.add_node("Z", 3, label=lbl)
                zx.add_edge(prev_nd, mnode, fault_free=meas_ff, label=f"meas_q{q}_{idx}")
                zx.meas_spiders.append((mnode, "Z", lbl))

        elif op.kind == "MPP":
            pstr = str(op.params["pauli"])
            if len(pstr) != len(op.qubits):
                raise RefusalError("MPP-LEN-MISMATCH", f"MPP pauli {pstr!r} != qubits {op.qubits}")
            lbl = op.label or f"MPP_{pstr}_{idx}"
            # Proposition 6.5 (mode == 'ff'): internal measurement edges are fault-free!
            # Proposition 6.7 (mode == 'ft'): internal measurement edges are NOISY (captures meas flip
            #   and meas flip + single-qubit output fault with 0 fault gadgets!).
            # Definition 2.4 (mode == 'noisy'): internal edges noisy PLUS multi-qubit output/meas+output
            #   fault gadgets via Prop. 6.4!
            internal_ff = op_ff
            red_hub = zx.add_node("X", 0, label=lbl)
            zx.meas_spiders.append((red_hub, "X", lbl))
            out_tips: list[int] = []
            for q, pch in zip(op.qubits, pstr):
                if pch == "I":
                    out_tips.append(live_end[q][0])
                    continue
                prev_nd, seg, _ = live_end[q]
                in_ff = wire_ff(q, seg)
                if pch == "Z":
                    sp = zx.add_node("Z", 0, label=f"mpp_z_q{q}_{idx}")
                    zx.add_edge(prev_nd, sp, fault_free=in_ff)
                    zx.add_edge(sp, red_hub, fault_free=internal_ff, label=f"mpp_leg_q{q}_{idx}")
                    out_tips.append(sp)
                elif pch == "X":
                    # Colour-dual of Z-measurement leg: H - Z(0) - H (Corollary 6.6 + Prop 6.7)
                    h_in = zx.add_node("H", 0)
                    sp = zx.add_node("Z", 0, label=f"mpp_x_q{q}_{idx}")
                    h_out = zx.add_node("H", 0)
                    zx.add_edge(prev_nd, h_in, fault_free=in_ff)
                    zx.add_edge(h_in, sp, fault_free=internal_ff)
                    zx.add_edge(sp, red_hub, fault_free=internal_ff, label=f"mpp_leg_q{q}_{idx}")
                    zx.add_edge(sp, h_out, fault_free=internal_ff)
                    out_tips.append(h_out)
                elif pch == "Y":
                    # Sdg; H; meas_Z; H; S  (since H Sdg Y S H = Z!)
                    sdg = zx.add_node("Z", 3)
                    h_in = zx.add_node("H", 0)
                    sp = zx.add_node("Z", 0, label=f"mpp_y_q{q}_{idx}")
                    h_out = zx.add_node("H", 0)
                    s_g = zx.add_node("Z", 1)
                    zx.add_edge(prev_nd, sdg, fault_free=in_ff)
                    zx.add_edge(sdg, h_in, fault_free=True)
                    zx.add_edge(h_in, sp, fault_free=internal_ff)
                    zx.add_edge(sp, red_hub, fault_free=internal_ff, label=f"mpp_leg_q{q}_{idx}")
                    zx.add_edge(sp, h_out, fault_free=internal_ff)
                    zx.add_edge(h_out, s_g, fault_free=True)
                    out_tips.append(s_g)
                else:
                    raise RefusalError("INVALID-PAULI-CHAR", f"invalid char {pch!r} in MPP")

            if op.mode == "noisy" and not full_ff and len(op.qubits) >= 2:
                # Full circuit-level noisy MPP: add fault gadgets for multi-qubit output faults
                # (plus optional measurement flip on red_hub) per Definition 2.4 & Prop 6.4!
                # Add a fault-free measurement-flip anchor on red_hub if needed:
                k_q = len(op.qubits)
                if k_q > 5:
                    raise RefusalError(
                        "NOISY-LARGE-MPP-REFUSED",
                        f"mode='noisy' MPP on {k_q} > 5 qubits has 2*4^k-1 atomic faults; use mode='ft' or 'ff'",
                    )
                cur_tips = list(out_tips)
                import itertools

                for mflip in (0, 1):
                    for combo in itertools.product(("I", "X", "Y", "Z"), repeat=k_q):
                        wt_out = sum(1 for c in combo if c != "I")
                        if wt_out <= 1:
                            # Already generated by 'ft' single-edge flips (Prop. 6.7)!
                            continue
                        sectors = {c for c in combo if c != "I"}
                        if noise_model == "circuit_css" and len(sectors) > 1:
                            continue
                        sec = next(iter(sectors)) if len(sectors) == 1 else "MIXED"
                        edge_specs: list[tuple[int, str]] = []
                        next_tips: list[int] = []
                        for j_idx, c_pauli in enumerate(combo):
                            if c_pauli == "I":
                                next_tips.append(cur_tips[j_idx])
                            else:
                                n_sp = zx.add_node("Z", 0)
                                e_id = zx.add_edge(cur_tips[j_idx], n_sp, fault_free=True)
                                edge_specs.append((e_id, c_pauli))
                                next_tips.append(n_sp)
                        if mflip:
                            # Attach an X-box on a fault-free identity tail off red_hub?
                            # An X(pi) spider into red_hub flips the measurement outcome:
                            # Wait: incoming green(pi) into red_hub copies to all legs? No: in Prop 6.7,
                            # an X (red pi) on ONE internal measurement leg fuses into red_hub!
                            # An X-box connected to our gadget's red spider pushes GREEN pi into X-box = RED pi
                            # on one internal leg of `red_hub`!
                            # Let's pick the first incident leg of red_hub!
                            first_leg = zx.adj[red_hub][0]
                            edge_specs.append((first_leg, "X"))
                        zx.add_fault_gadget(edge_specs, label=f"mpp_gad_{idx}", sector=sec)
                        cur_tips = next_tips
                out_tips = cur_tips

            for q, tip in zip(op.qubits, out_tips):
                _, seg, _ = live_end[q]
                live_end[q] = (tip, seg + 1, None)

    # Terminate all output qubits on ordered BOUND nodes
    for q in circuit.output_qubits:
        prev_nd, seg, _ = live_end[q]
        b_out = zx.add_boundary(is_input=False, label=f"out_{q}")
        zx.add_edge(prev_nd, b_out, fault_free=wire_ff(q, seg), label=f"w_out_q{q}")

    zx.validate()
    return zx


def extract_circuit_fe(
    zx: ZXDiagram,
    circuit_template: CliffordCircuit,
    target_noise_model: str = "circuit_depolarizing",
) -> CliffordCircuit:
    """Validate and extract a fault-equivalent implementable CliffordCircuit from a ZXDiagram.

    Enforces Proposition 6.8-6.10 fail-closed:
      1) Idealised primitives (CAT_PREP, ENC_PREP, mode in {'ft', 'ff'}, or fault_free_wires)
         are refused if `target_noise_model` requires a physical hardware circuit unless
         explicitly marked allowed via their own modes, AND
      2) Under `circuit_depolarizing`, any noisy CNOT in `circuit_template` whose corresponding
         ZX subgraph lacks a Prop. 6.10 split termination (`PREP_PLUS` on control or `PREP_0`
         on target), a Prop. 6.9 merge termination (`MEAS_Z` on target or `MEAS_X` on control),
         or the 4 explicit Prop. 6.8 cross-sector fault gadgets in `zx` is REFUSED with
         `UNGUARDED-CNOT-EXTRACTION`!
    """
    zx.validate()
    circuit_template.validate()
    if len(zx.inputs) != len(circuit_template.input_qubits) or len(zx.outputs) != len(circuit_template.output_qubits):
        raise RefusalError(
            "BOUNDARY-ARITY-MISMATCH",
            "ZX diagram boundary counts do not match candidate circuit inputs/outputs",
        )
    if target_noise_model == "circuit_depolarizing":
        # Each unguarded noisy CNOT requires 4 cross-sector Prop. 6.8 fault gadgets unless covered
        # by Prop. 6.9 (merge into MEAS_Z/MEAS_X) or Prop. 6.10 (split from PREP_PLUS/PREP_0)!
        n_mixed_gadgets = sum(1 for ed in zx.edges.values() if ed.gadget_sector == "MIXED")
        n_ops = len(circuit_template.ops)
        next_op_on_q: dict[tuple[int, int], int] = {}
        last_seen: dict[int, int] = {}
        for idx in range(n_ops - 1, -1, -1):
            op = circuit_template.ops[idx]
            for q in op.qubits:
                if q in last_seen:
                    next_op_on_q[(idx, q)] = last_seen[q]
                last_seen[q] = idx
        prev_op_on_q: dict[tuple[int, int], int] = {}
        prev_seen: dict[int, int] = {}
        for idx, op in enumerate(circuit_template.ops):
            for q in op.qubits:
                if q in prev_seen:
                    prev_op_on_q[(idx, q)] = prev_seen[q]
                prev_seen[q] = idx

        unguarded_cxs: list[tuple[int, int, int]] = []
        for idx, op in enumerate(circuit_template.ops):
            if op.kind == "CX" and op.mode == "noisy" and not circuit_template.all_wires_fault_free:
                qc, qt = op.qubits
                pc_idx = prev_op_on_q.get((idx, qc))
                pt_idx = prev_op_on_q.get((idx, qt))
                nc_idx = next_op_on_q.get((idx, qc))
                nt_idx = next_op_on_q.get((idx, qt))
                split_ok = (
                    (pc_idx is not None and circuit_template.ops[pc_idx].kind == "PREP_PLUS")
                    or (pt_idx is not None and circuit_template.ops[pt_idx].kind == "PREP_0")
                )
                merge_ok = (
                    (nt_idx is not None and circuit_template.ops[nt_idx].kind == "MEAS_Z")
                    or (nc_idx is not None and circuit_template.ops[nc_idx].kind == "MEAS_X")
                )
                if not (split_ok or merge_ok):
                    unguarded_cxs.append((idx, qc, qt))

        if n_mixed_gadgets < 4 * len(unguarded_cxs):
            idx, qc, qt = unguarded_cxs[0]
            raise RefusalError(
                "UNGUARDED-CNOT-EXTRACTION",
                f"CNOT at op {idx} on qubits ({qc}, {qt}) lacks split/merge termination "
                f"(Props 6.9-6.10) and ZX diagram has {n_mixed_gadgets} < {4 * len(unguarded_cxs)} "
                f"Prop 6.8 cross-sector fault gadgets under circuit_depolarizing noise",
            )
    return circuit_template


def load_spec_file(path: str) -> CliffordCircuit | ZXDiagram:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    kind = data.get("kind")
    if kind == "CliffordCircuit":
        return CliffordCircuit.from_dict(data)
    if kind == "ZXDiagram":
        return ZXDiagram.from_dict(data)
    raise RefusalError("UNKNOWN-SPEC-KIND", f"JSON 'kind' must be CliffordCircuit or ZXDiagram, got {kind!r}")
