"""Leg B, part 3: drive CaDiCaL and check its LRAT proof with lrat.check.

The solver is not trusted.  UNSAT is accepted only when the in-package checker
accepts the proof against the CNF this package wrote; SAT is accepted only
when the returned model decodes to a word that passes the direct test.
"""
import hashlib
import os
import shutil
import subprocess
import tempfile
import time

from . import certify, cnf, gf2, lrat

QUANTITIES = ("d_x", "d_z", "wmin_x", "wmin_z")


def find_cadical():
    cand = os.environ.get("QECCERT_CADICAL")
    if cand:
        return cand if os.path.isfile(cand) and os.access(cand, os.X_OK) else None
    return shutil.which("cadical")


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


_FLAGS = {}


def _model_flags(solver):
    """Recent CaDiCaL omits some variables from the printed model unless asked;
    a partial model would be misread, so ask for all literals where supported."""
    if solver not in _FLAGS:
        try:
            h = subprocess.run([solver, "--help"], capture_output=True, text=True, timeout=10).stdout
        except (OSError, subprocess.TimeoutExpired):
            h = ""
        _FLAGS[solver] = ["--modelalllits=true"] if "modelalllits" in h else []
    return _FLAGS[solver]


def solve(formula, workdir, tag, timeout, solver=None):
    """Returns ("unsat", info) | ("sat", word, info) | ("unknown", reason)."""
    solver = solver or find_cadical()
    if not solver:
        return ("unknown", "cadical not found (put it on PATH or set QECCERT_CADICAL)")
    cnf_path = os.path.join(workdir, f"{tag}.cnf")
    proof_path = os.path.join(workdir, f"{tag}.lrat")
    text = formula.dimacs()
    with open(cnf_path, "w") as f:
        f.write(text)
    t0 = time.monotonic()
    try:
        p = subprocess.run([solver, "--lrat", "--no-binary", "-q", *_model_flags(solver),
                            cnf_path, proof_path],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return ("unknown", f"cadical exceeded {timeout} s")
    secs = round(time.monotonic() - t0, 2)
    info = dict(cnf_sha256=hashlib.sha256(text.encode()).hexdigest(),
                vars=formula.nvars, clauses=len(formula.clauses), solve_seconds=secs)
    if p.returncode == 20:
        with open(proof_path, "rb") as f:
            proof = f.read()
        t1 = time.monotonic()
        try:
            st = lrat.check(text, proof)
        except lrat.ProofError as e:
            return ("unknown", f"cadical said UNSAT but its LRAT proof was rejected: {e}")
        info.update(lrat_sha256=_sha(proof_path), lrat_bytes=len(proof),
                    lrat_steps=st["steps"], lrat_hints=st["hints"],
                    check_seconds=round(time.monotonic() - t1, 2))
        return ("unsat", info)
    if p.returncode == 10:
        word, seen = 0, set()
        for line in p.stdout.splitlines():
            if line.startswith("v "):
                for tok in line[2:].split():
                    v = int(tok)
                    seen.add(abs(v))
                    if 0 < v <= formula.n:
                        word |= 1 << (v - 1)
        if not all(j in seen for j in range(1, formula.n + 1)):
            return ("unknown", "cadical said SAT but printed a partial model")
        if not formula.holds(word):
            return ("unknown", "cadical said SAT but its model fails the direct test")
        return ("sat", word, info)
    return ("unknown", f"cadical exit code {p.returncode}: {(p.stderr or p.stdout)[-200:]}")


def solver_version(solver=None):
    solver = solver or find_cadical()
    if not solver:
        return None
    try:
        return subprocess.run([solver, "--version"], capture_output=True, text=True,
                              timeout=10).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        return None


def run_leg(receipt, code, targets=None, workdir=None, timeout=600, quantities=QUANTITIES):
    """Prove q >= b for each quantity, where b is the caller's target or else
    the witnessed upper bound already in the receipt."""
    targets = targets or {}
    tmp = None
    if workdir is None:
        tmp = tempfile.mkdtemp(prefix="qeccert-sat-")
        workdir = tmp
    os.makedirs(workdir, exist_ok=True)
    ver = solver_version()
    try:
        for q in quantities:
            entry = receipt["facts"][q]
            ev = dict(leg="B", quantity=q, method="SAT (sequential-counter CNF), LRAT checked "
                                                  "by qeccert.lrat", solver=ver)
            b = targets.get(q, entry["upper"])
            if b is None:
                ev.update(status="skipped", reason="no target and no witnessed upper bound")
                receipt["evidence"].append(ev)
                continue
            if b < 1:
                ev.update(status="skipped", reason="target below 1")
                receipt["evidence"].append(ev)
                continue
            res = solve(cnf.for_quantity(code, q, b - 1), workdir,
                        f"{code.sha256()[:12]}_{q}_lt{b}", timeout)
            if res[0] == "unsat":
                ev.update(status="proved", claim=f"{q} >= {b}", **res[1])
                certify.add_lower_bound(receipt, q, b, ev)
            elif res[0] == "sat":
                word = res[1]
                wgt = gf2.wt(word)
                ev.update(status="refuted", reason=f"found a word of weight {wgt} < {b}",
                          **res[2])
                receipt["evidence"].append(ev)
                _set_upper(receipt, code, q, word)
            else:
                ev.update(status="skipped", reason=res[1])
                receipt["evidence"].append(ev)
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return receipt


def _set_upper(receipt, code, q, word):
    e = receipt["facts"][q]
    w = gf2.wt(word)
    if e["lower"] is not None and w < e["lower"]:
        raise certify.Refusal(f"{q}: a verified word of weight {w} contradicts the proven "
                              f"lower bound {e['lower']}: the legs disagree")
    if e["upper"] is None or w < e["upper"]:
        n = code.n
        e["upper"], e["witness"] = w, certify._bits(word, n)
        if q in ("d_x", "d_z"):
            dual = gf2.nullspace(code.hx if q == "d_x" else code.hz, n)
            p = certify._partner(word, dual)
            e["partner"] = certify._bits(p, n) if p is not None else None
        e["exact"] = e["lower"] == e["upper"]
        certify._combine(receipt["facts"])
