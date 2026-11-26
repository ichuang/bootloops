"""Leg B, part 2: an LRAT proof checker (stdlib only, fail-closed).

Checks a refutation in the LRAT format (Cruz-Filipe, Heule, Hunt, Kaufmann,
Schneider-Kamp 2017), text or binary, against a DIMACS CNF.  Only the RUP
fragment is accepted: every added clause must follow by unit propagation
along its listed hints, each hint being unit or falsified at its turn.  A
negative (RAT) hint is refused, because accepting RAT would need a second
argument about satisfiability preservation that this checker does not carry;
CaDiCaL's LRAT output does not use it.

A proof is accepted only if it derives the empty clause.
"""


class ProofError(Exception):
    pass


def parse_dimacs(text):
    nvars, clauses, cur, declared = 0, [], [], None
    for line in text.splitlines():
        line = line.strip()
        if not line or line[0] == "c":
            continue
        if line[0] == "p":
            parts = line.split()
            if len(parts) != 4 or parts[1] != "cnf":
                raise ProofError("bad DIMACS header")
            nvars, declared = int(parts[2]), int(parts[3])
            continue
        for tok in line.split():
            v = int(tok)
            if v == 0:
                clauses.append(tuple(cur))
                cur = []
            else:
                if abs(v) > nvars:
                    raise ProofError("literal exceeds declared variable count")
                cur.append(v)
    if cur:
        raise ProofError("unterminated clause")
    if declared is None or declared != len(clauses):
        raise ProofError("clause count does not match DIMACS header")
    return nvars, clauses


def _steps_text(data):
    for ln, line in enumerate(data.decode("ascii").splitlines(), 1):
        parts = line.split()
        if not parts:
            continue
        try:
            if len(parts) > 1 and parts[1] == "d":
                ids = [int(x) for x in parts[2:]]
                if not ids or ids[-1] != 0:
                    raise ValueError
                yield ("d", ids[:-1])
                continue
            nums = [int(x) for x in parts]
            z = nums.index(0, 1)
            if nums[-1] != 0:
                raise ValueError
            yield ("a", nums[0], nums[1:z], nums[z + 1:-1])
        except ValueError:
            raise ProofError(f"malformed LRAT line {ln}") from None


def _steps_binary(data):
    pos, end = 0, len(data)

    def varint():
        nonlocal pos
        x, shift = 0, 0
        while True:
            if pos >= end:
                raise ProofError("truncated binary LRAT")
            b = data[pos]
            pos += 1
            x |= (b & 0x7F) << shift
            if not b & 0x80:
                return x
            shift += 7

    def signed_list():
        out = []
        while True:
            u = varint()
            if u == 0:
                return out
            out.append(-(u >> 1) if u & 1 else (u >> 1))

    while pos < end:
        tag = data[pos]
        pos += 1
        if tag == 0x61:      # 'a'
            cid = varint()
            lits = signed_list()
            hints = signed_list()
            yield ("a", cid, lits, hints)
        elif tag == 0x64:    # 'd'
            yield ("d", signed_list())
        else:
            raise ProofError(f"unknown binary LRAT tag {tag:#x} at byte {pos - 1}")


def steps(data):
    if not data:
        raise ProofError("empty proof")
    return _steps_binary(data) if data[0] in (0x61, 0x64) else _steps_text(data)


def check(cnf_text, proof_bytes):
    """Return statistics if the proof refutes the CNF; raise ProofError otherwise."""
    nvars, clauses = parse_dimacs(cnf_text)
    db = {i + 1: c for i, c in enumerate(clauses)}
    if any(len(c) == 0 for c in clauses):
        return dict(steps=0, hints=0, trivially_unsat=True)
    val = {}          # literal -> True when assigned true
    nsteps = nhints = 0
    last_id = len(clauses)
    for st in steps(proof_bytes):
        if st[0] == "d":
            for i in st[1]:
                db.pop(i, None)
            continue
        _, cid, lits, hints = st
        if cid <= last_id:
            raise ProofError(f"clause id {cid} is not increasing")
        last_id = cid
        val.clear()
        conflict = False
        for l in lits:
            if abs(l) > nvars:
                # extension variables would need a definitional argument
                raise ProofError(f"clause {cid} uses a variable outside the CNF")
            if val.get(l):
                conflict = True      # tautology: l and -l both in the clause
            val[-l] = True
        if not conflict:
            for h in hints:
                nhints += 1
                if h <= 0:
                    raise ProofError(f"clause {cid}: RAT hint {h} is not accepted")
                c = db.get(h)
                if c is None:
                    raise ProofError(f"clause {cid}: hint {h} is unknown or deleted")
                unit = 0
                for l in c:
                    if val.get(l):
                        raise ProofError(f"clause {cid}: hint {h} is satisfied, not unit")
                    if not val.get(-l):
                        if unit and unit != l:
                            raise ProofError(f"clause {cid}: hint {h} is not unit")
                        unit = l
                if unit == 0:
                    conflict = True
                    break
                val[unit] = True
        if not conflict:
            raise ProofError(f"clause {cid}: hints do not reach a conflict")
        db[cid] = tuple(lits)
        nsteps += 1
        if not lits:
            return dict(steps=nsteps, hints=nhints, trivially_unsat=False)
    raise ProofError("proof ends without deriving the empty clause")
