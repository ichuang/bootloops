"""External CSS code interchange, sparse/Matrix-Market loaders, and coordinate maps.

Schema `ftqc_codes-import-1` records producer provenance, field and code type,
source SHA-256 digests, index base (0 or 1), explicit bijective coordinate and
check-row permutations, optional paired candidate logical operators, and
optional finite distance witnesses while preserving redundant physical check
rows and exact check order.
"""
import hashlib
import json
import os

from . import certify, codes, gf2

SCHEMA = "ftqc_codes-import-1"


def _is_strict_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def _sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _bits(v, n):
    return "".join("1" if (v >> j) & 1 else "0" for j in range(n))


def _unbits(s, n):
    if not isinstance(s, str) or len(s) != n or set(s) - {"0", "1"}:
        raise codes.SpecError(f"expected 0/1 string of length {n}")
    return int(s[::-1], 2)


def validate_permutation(perm, size, label="coordinate"):
    """Validate that `perm` is an explicit bijection on `0..size-1`."""
    if not isinstance(perm, (list, tuple)):
        raise codes.SpecError(f"{label} map must be a list of integers")
    if len(perm) != size or not all(_is_strict_int(x) for x in perm):
        raise certify.Refusal(f"{label} map must have {size} integer entries")
    if set(perm) != set(range(size)):
        raise certify.Refusal(f"{label} map is not a bijection on 0..{size - 1}")
    return [int(x) for x in perm]


def transport_word(word, perm_q):
    """Push bitmask `word` along external-to-native qubit map `j -> perm_q[j]`."""
    return gf2.from_supp(perm_q[j] for j in gf2.supp(word))


def apply_coordinate_map(n, hx, hz, cmap):
    """Apply external-to-native `(perm_q, perm_hx, perm_hz)` preserving row count."""
    if not isinstance(cmap, dict):
        raise codes.SpecError("coordinate_map must be an object")
    unknown = set(cmap) - {"perm_q", "perm_hx", "perm_hz", "convention"}
    if unknown:
        raise codes.SpecError(f"unknown coordinate_map fields: {sorted(unknown)}")
    perm_q = validate_permutation(cmap.get("perm_q", list(range(n))), n, "perm_q")
    perm_hx = validate_permutation(
        cmap.get("perm_hx", list(range(len(hx)))), len(hx), "perm_hx"
    )
    perm_hz = validate_permutation(
        cmap.get("perm_hz", list(range(len(hz)))), len(hz), "perm_hz"
    )
    hx_out = [0] * len(hx)
    for i, r in enumerate(hx):
        hx_out[perm_hx[i]] = transport_word(r, perm_q)
    hz_out = [0] * len(hz)
    for i, r in enumerate(hz):
        hz_out[perm_hz[i]] = transport_word(r, perm_q)
    return hx_out, hz_out, {
        "convention": "push_external_to_native",
        "perm_q": perm_q,
        "perm_hx": perm_hx,
        "perm_hz": perm_hz,
    }


def parse_matrix_market(text, *, index_base=1, allow_duplicates=False):
    """Parse a binary Matrix Market coordinate matrix into `(nrows, ncols, supports)`."""
    if not _is_strict_int(index_base) or index_base not in (0, 1):
        raise codes.SpecError(f"index_base must be 0 or 1, got {index_base!r}")
    raw_lines = text.splitlines()
    if not raw_lines or not raw_lines[0].startswith("%%MatrixMarket"):
        raise codes.SpecError("missing %%MatrixMarket header")
    hdr = raw_lines[0].strip().split()
    if len(hdr) != 5 or hdr[1] != "matrix" or hdr[2] != "coordinate":
        raise codes.SpecError(f"unsupported Matrix Market header: {raw_lines[0]!r}")
    field_kind, sym_kind = hdr[3].lower(), hdr[4].lower()
    if field_kind not in ("integer", "pattern"):
        raise certify.Refusal(
            f"binary CSS Matrix Market import requires integer or pattern field, got {field_kind!r}"
        )
    if sym_kind != "general":
        raise codes.SpecError(f"unsupported Matrix Market symmetry {sym_kind!r}; expected 'general'")

    body = [ln.strip() for ln in raw_lines[1:] if ln.strip() and not ln.lstrip().startswith("%")]
    if not body:
        raise codes.SpecError("Matrix Market file lacks dimension line")
    dims = body[0].split()
    if len(dims) != 3:
        raise codes.SpecError(f"invalid Matrix Market dimension line: {body[0]!r}")
    try:
        nrows, ncols, nnz = int(dims[0]), int(dims[1]), int(dims[2])
    except ValueError:
        raise codes.SpecError(f"non-integer Matrix Market dimensions: {body[0]!r}") from None
    if nrows <= 0 or ncols <= 0 or nnz < 0:
        raise codes.SpecError(f"invalid Matrix Market shape {nrows}x{ncols} nnz={nnz}")
    entries = body[1:]
    if len(entries) != nnz:
        raise codes.SpecError(f"Matrix Market declares nnz={nnz} but has {len(entries)} entries")

    row_sets = [set() for _ in range(nrows)]
    for line_idx, entry in enumerate(entries, start=1):
        parts = entry.split()
        expected_cols = 2 if field_kind == "pattern" else 3
        if len(parts) != expected_cols:
            raise codes.SpecError(f"malformed Matrix Market entry #{line_idx}: {entry!r}")
        try:
            r_raw, c_raw = int(parts[0]), int(parts[1])
            val = 1 if field_kind == "pattern" else int(parts[2])
        except ValueError:
            raise codes.SpecError(f"non-integer Matrix Market entry #{line_idx}: {entry!r}") from None
        if val != 1:
            raise certify.Refusal(
                f"binary Matrix Market entry #{line_idx} has value {val}; expected 1 "
                "(non-binary or higher-field matrix requires explicit binary expansion)"
            )
        r = r_raw - index_base
        c = c_raw - index_base
        if not (0 <= r < nrows) or not (0 <= c < ncols):
            raise certify.Refusal(
                f"Matrix Market coordinate ({r_raw}, {c_raw}) out of range for "
                f"{nrows}x{ncols} with index_base={index_base}"
            )
        if c in row_sets[r]:
            if not allow_duplicates:
                raise certify.Refusal(
                    f"duplicate Matrix Market coordinate ({r_raw}, {c_raw}) in row {r_raw}"
                )
            row_sets[r].remove(c)
        else:
            row_sets[r].add(c)
    return nrows, ncols, [sorted(s) for s in row_sets]


def parse_sparse_supports(supports, nrows, ncols, *, index_base=0, allow_duplicates=False, label="matrix"):
    """Validate sparse row supports and convert to ordered bitmask rows."""
    if not _is_strict_int(index_base) or index_base not in (0, 1):
        raise codes.SpecError(f"{label}: index_base must be 0 or 1")
    if not _is_strict_int(nrows) or nrows <= 0 or not _is_strict_int(ncols) or ncols <= 0:
        raise codes.SpecError(f"{label}: nrows and ncols must be positive integers")
    if not isinstance(supports, list) or len(supports) != nrows:
        raise codes.SpecError(f"{label}: expected {nrows} sparse row supports")
    out = []
    for r_idx, row in enumerate(supports):
        if not isinstance(row, list):
            raise codes.SpecError(f"{label} row {r_idx}: support must be a list")
        seen = set()
        for c_raw in row:
            if not _is_strict_int(c_raw):
                raise codes.SpecError(f"{label} row {r_idx}: column index must be an integer")
            c = c_raw - index_base
            if not (0 <= c < ncols):
                raise certify.Refusal(
                    f"{label} row {r_idx}: coordinate {c_raw} out of bounds for width {ncols} "
                    f"(index_base={index_base})"
                )
            if c in seen:
                if not allow_duplicates:
                    raise certify.Refusal(
                        f"{label} row {r_idx}: duplicate coordinate {c_raw}"
                    )
                seen.remove(c)
            else:
                seen.add(c)
        out.append(gf2.from_supp(seen))
    return out


def load_matrix_block(block, *, expected_n=None, default_index_base=0, base_dir=None, label="matrix"):
    """Load a dense, sparse-support, or Matrix Market block into `(nrows, ncols, bit_rows, meta)`."""
    if not isinstance(block, dict) or "format" not in block:
        raise codes.SpecError(f"{label} block requires 'format'")
    fmt = block["format"]
    if fmt == "dense_strings":
        rows_raw = block.get("rows")
        if not isinstance(rows_raw, list) or not rows_raw:
            raise codes.SpecError(f"{label}: dense_strings requires non-empty 'rows'")
        ncols = len(rows_raw[0]) if isinstance(rows_raw[0], str) else -1
        if ncols <= 0:
            raise codes.SpecError(f"{label}: invalid dense row")
        bit_rows = [_unbits(s, ncols) for s in rows_raw]
        raw_bytes = ("\n".join(rows_raw) + "\n").encode("utf-8")
        actual_sha = _sha256_bytes(raw_bytes)
        if "sha256" in block and block["sha256"] != actual_sha:
            raise certify.Refusal(f"{label} dense sha256 mismatch")
        meta = {"format": fmt, "rows": len(bit_rows), "n": ncols, "sha256": actual_sha}
    elif fmt == "sparse_rows":
        if "rows" not in block or "n" not in block or "supports" not in block:
            raise codes.SpecError(f"{label}: sparse_rows requires 'rows', 'n', and 'supports'")
        nrows, ncols = block["rows"], block["n"]
        ibase = block.get("index_base", default_index_base)
        allow_dup = bool(block.get("allow_duplicates", False))
        bit_rows = parse_sparse_supports(
            block["supports"], nrows, ncols, index_base=ibase, allow_duplicates=allow_dup, label=label
        )
        canonical_bytes = json.dumps(
            [[j for j in gf2.supp(r)] for r in bit_rows], separators=(",", ":")
        ).encode("utf-8")
        actual_sha = _sha256_bytes(canonical_bytes)
        if "sha256" in block and block["sha256"] != actual_sha:
            raise certify.Refusal(f"{label} sparse_rows sha256 mismatch")
        meta = {
            "format": fmt,
            "rows": nrows,
            "n": ncols,
            "index_base": ibase,
            "sha256": actual_sha,
        }
    elif fmt == "matrix_market":
        ibase = block.get("index_base", 1)
        allow_dup = bool(block.get("allow_duplicates", False))
        if "path" in block:
            rel_path = block["path"]
            if not isinstance(rel_path, str) or not rel_path:
                raise codes.SpecError(f"{label}: invalid Matrix Market path")
            full_path = os.path.join(base_dir or ".", rel_path)
            with open(full_path, "rb") as f:
                raw_bytes = f.read()
            src_ref = rel_path
        elif "text" in block:
            if not isinstance(block["text"], str):
                raise codes.SpecError(f"{label}: Matrix Market 'text' must be a string")
            raw_bytes = block["text"].encode("utf-8")
            src_ref = "inline"
        else:
            raise codes.SpecError(f"{label}: matrix_market requires 'path' or 'text'")
        actual_sha = _sha256_bytes(raw_bytes)
        if "sha256" in block and block["sha256"] != actual_sha:
            raise certify.Refusal(
                f"{label} file sha256 mismatch: manifest expects {block['sha256']}, got {actual_sha}"
            )
        nrows, ncols, supports = parse_matrix_market(
            raw_bytes.decode("utf-8"), index_base=ibase, allow_duplicates=allow_dup
        )
        bit_rows = [gf2.from_supp(s) for s in supports]
        meta = {
            "format": fmt,
            "source": src_ref,
            "rows": nrows,
            "n": ncols,
            "index_base": ibase,
            "sha256": actual_sha,
        }
    else:
        raise codes.SpecError(f"{label}: unsupported matrix format {fmt!r}")

    if expected_n is not None and ncols != expected_n:
        raise certify.Refusal(f"{label} width {ncols} does not match code n={expected_n}")
    return len(bit_rows), ncols, bit_rows, meta


def verify_candidate_logicals(
    code,
    lx,
    lz,
    *,
    signs_x=None,
    signs_z=None,
    slot_names=None,
    convention="positive_css",
):
    """Verify paired X/Z logical representatives against `code`."""
    if not code.commutes():
        raise certify.Refusal("cannot verify candidate logicals on non-commuting checks")
    if convention not in ("positive_css", "signed_css"):
        raise codes.SpecError(f"unknown candidate_logicals convention {convention!r}")
    ex, ez = gf2.Echelon(code.hx), gf2.Echelon(code.hz)
    k = code.n - ex.rank - ez.rank
    if len(lx) != k or len(lz) != k:
        raise certify.Refusal(
            f"candidate logical count ({len(lx)}, {len(lz)}) does not match k={k}"
        )
    for idx, xw in enumerate(lx):
        if xw <= 0 or xw >> code.n:
            raise certify.Refusal(f"logical X[{idx}] support out of range")
        if any(gf2.dot(xw, h) for h in code.hz):
            raise certify.Refusal(f"logical X[{idx}] anticommutes with a Z check")
    for idx, zw in enumerate(lz):
        if zw <= 0 or zw >> code.n:
            raise certify.Refusal(f"logical Z[{idx}] support out of range")
        if any(gf2.dot(zw, h) for h in code.hx):
            raise certify.Refusal(f"logical Z[{idx}] anticommutes with an X check")

    ex_ext = gf2.Echelon(code.hx)
    for idx, xw in enumerate(lx):
        if not ex_ext.add(xw):
            raise certify.Refusal(f"logical X[{idx}] is dependent modulo X stabilizers")
    ez_ext = gf2.Echelon(code.hz)
    for idx, zw in enumerate(lz):
        if not ez_ext.add(zw):
            raise certify.Refusal(f"logical Z[{idx}] is dependent modulo Z stabilizers")

    for i in range(k):
        for j in range(k):
            want = 1 if i == j else 0
            if gf2.dot(lx[i], lz[j]) != want:
                raise certify.Refusal(
                    f"candidate logical pairing L_X[{i}] . L_Z[{j}] = {gf2.dot(lx[i], lz[j])}, expected {want}"
                )

    if signs_x is None:
        if convention == "signed_css":
            raise codes.SpecError("signed_css convention requires explicit signs_x and signs_z")
        sx = [0] * k
    else:
        if not isinstance(signs_x, list) or len(signs_x) != k or any(s not in (0, 1) for s in signs_x):
            raise codes.SpecError("signs_x must be a list of k bits (0 or 1)")
        sx = [int(s) for s in signs_x]

    if signs_z is None:
        if convention == "signed_css":
            raise codes.SpecError("signed_css convention requires explicit signs_x and signs_z")
        sz = [0] * k
    else:
        if not isinstance(signs_z, list) or len(signs_z) != k or any(s not in (0, 1) for s in signs_z):
            raise codes.SpecError("signs_z must be a list of k bits (0 or 1)")
        sz = [int(s) for s in signs_z]

    if slot_names is None:
        names = [f"q{i}" for i in range(k)]
    else:
        if (
            not isinstance(slot_names, list)
            or len(slot_names) != k
            or any(not isinstance(nm, str) or not nm for nm in slot_names)
            or len(set(slot_names)) != k
        ):
            raise certify.Refusal("slot_names must be k unique non-empty strings")
        names = list(slot_names)

    return {
        "convention": convention,
        "slot_names": names,
        "signs_x": sx,
        "signs_z": sz,
        "lx_supports": [gf2.supp(w) for w in lx],
        "lz_supports": [gf2.supp(w) for w in lz],
    }


def _parse_witness_entry(entry, code, quantity, perm_q, default_index_base=0):
    if not isinstance(entry, dict) or "support" not in entry:
        raise codes.SpecError(f"witness {quantity} requires 'support'")
    ibase = entry.get("index_base", default_index_base)
    coords = entry.get("coordinates", "native")
    if coords not in ("external", "native"):
        raise codes.SpecError(f"witness {quantity}: coordinates must be 'external' or 'native'")
    w = parse_sparse_supports(
        [entry["support"]], 1, code.n, index_base=ibase, allow_duplicates=False, label=f"witness {quantity}"
    )[0]
    if w == 0:
        raise certify.Refusal(f"witness {quantity} cannot be empty")
    if coords == "external":
        w = transport_word(w, perm_q)
    checks = code.hz if quantity == "d_x" else code.hx
    other = code.hx if quantity == "d_x" else code.hz
    if any(gf2.dot(w, h) for h in checks):
        raise certify.Refusal(f"witness {quantity} fails opposite parity check")
    partner = None
    if "partner_support" in entry and entry["partner_support"] is not None:
        p = parse_sparse_supports(
            [entry["partner_support"]], 1, code.n, index_base=ibase, allow_duplicates=False,
            label=f"witness {quantity} partner",
        )[0]
        if coords == "external":
            p = transport_word(p, perm_q)
        partner = p
    else:
        for z in gf2.nullspace(other, code.n):
            if gf2.dot(z, w):
                partner = z
                break
    if partner is None or any(gf2.dot(partner, h) for h in other) or gf2.dot(partner, w) != 1:
        raise certify.Refusal(f"witness {quantity} is a stabilizer or has invalid partner")
    return {
        "weight": gf2.wt(w),
        "support": gf2.supp(w),
        "partner_support": gf2.supp(partner),
    }


def load_import_manifest(manifest_or_path, *, base_dir=None):
    """Load and validate an `ftqc_codes-import-1` JSON manifest into a `CSSCode`."""
    if isinstance(manifest_or_path, str):
        full_path = os.path.abspath(manifest_or_path)
        if base_dir is None:
            base_dir = os.path.dirname(full_path)
        with open(full_path) as f:
            d = json.load(f)
    elif isinstance(manifest_or_path, dict):
        d = manifest_or_path
    else:
        raise codes.SpecError("manifest must be a dict or file path")

    if d.get("schema") != SCHEMA:
        raise codes.SpecError(f"expected schema {SCHEMA!r}, got {d.get('schema')!r}")

    code_type = str(d.get("code_type", "")).lower()
    if code_type in ("non_css", "non-css", "symplectic", "stabilizer_non_css"):
        raise certify.Refusal(
            "non-CSS stabilizer code requires a separate symplectic-code extension; "
            "ftqc_codes CSS interchange only admits CSS checks"
        )
    if code_type != "css":
        raise codes.SpecError(f"unsupported code_type {d.get('code_type')!r}; expected 'css'")

    field = d.get("field")
    if field not in ("GF(2)", 2):
        raise certify.Refusal(
            f"unsupported field {field!r}; binary CSS import requires field 'GF(2)' "
            "(expand extension-field checks to binary before import)"
        )

    producer = d.get("producer")
    if (
        not isinstance(producer, dict)
        or not isinstance(producer.get("name"), str)
        or not producer["name"].strip()
        or not (producer.get("version") or producer.get("revision"))
    ):
        raise codes.SpecError("import manifest requires producer.name and producer.version or revision")

    params = d.get("parameters")
    if not isinstance(params, dict) or not _is_strict_int(params.get("n")) or params["n"] <= 0:
        raise codes.SpecError("import manifest requires positive integer parameters.n")
    n = params["n"]
    default_ibase = d.get("index_base", 0)
    if not _is_strict_int(default_ibase) or default_ibase not in (0, 1):
        raise codes.SpecError("index_base must be 0 or 1")

    sources = d.get("sources")
    if not isinstance(sources, dict) or "hx" not in sources or "hz" not in sources:
        raise codes.SpecError("import manifest requires sources.hx and sources.hz")

    _, _, hx_raw, hx_meta = load_matrix_block(
        sources["hx"], expected_n=n, default_index_base=default_ibase, base_dir=base_dir, label="hx"
    )
    _, _, hz_raw, hz_meta = load_matrix_block(
        sources["hz"], expected_n=n, default_index_base=default_ibase, base_dir=base_dir, label="hz"
    )

    name = d.get("name") or "imported_css"
    raw_code = codes.CSSCode(name, n, hx_raw, hz_raw)
    raw_sha = raw_code.sha256()
    if "raw_sha256" in d and d["raw_sha256"] is not None and d["raw_sha256"] != raw_sha:
        raise certify.Refusal(
            f"raw_sha256 mismatch: manifest expects {d['raw_sha256']}, got {raw_sha}"
        )

    if "coordinate_map" in d and d["coordinate_map"] is not None:
        hx, hz, cmap_meta = apply_coordinate_map(n, hx_raw, hz_raw, d["coordinate_map"])
        perm_q = cmap_meta["perm_q"]
    else:
        hx, hz = hx_raw, hz_raw
        perm_q = list(range(n))
        cmap_meta = None

    code = codes.CSSCode(name, n, hx, hz)
    if "sha256" in d and d["sha256"] is not None and d["sha256"] != code.sha256():
        raise certify.Refusal(
            f"imported matrix sha256 {code.sha256()} does not match declared {d['sha256']}"
        )
    if not code.commutes():
        raise certify.Refusal("imported H_X and H_Z do not commute")

    ex, ez = gf2.Echelon(hx), gf2.Echelon(hz)
    k = n - ex.rank - ez.rank
    if "k" in params and params["k"] is not None:
        if not _is_strict_int(params["k"]) or params["k"] != k:
            raise certify.Refusal(f"declared parameters.k={params['k']} disagrees with computed k={k}")

    source_hashes = {"hx": hx_meta["sha256"], "hz": hz_meta["sha256"]}
    cand_meta = None
    if d.get("candidate_logicals") is not None:
        cl = d["candidate_logicals"]
        if not isinstance(cl, dict) or "lx" not in cl or "lz" not in cl:
            raise codes.SpecError("candidate_logicals requires 'lx' and 'lz' matrix blocks")
        _, _, lx_raw, lx_m = load_matrix_block(
            cl["lx"], expected_n=n, default_index_base=default_ibase, base_dir=base_dir, label="lx"
        )
        _, _, lz_raw, lz_m = load_matrix_block(
            cl["lz"], expected_n=n, default_index_base=default_ibase, base_dir=base_dir, label="lz"
        )
        source_hashes["lx"] = lx_m["sha256"]
        source_hashes["lz"] = lz_m["sha256"]
        coords = cl.get("coordinates", "external")
        if coords not in ("external", "native"):
            raise codes.SpecError("candidate_logicals.coordinates must be 'external' or 'native'")
        if coords == "external" and cmap_meta is not None:
            lx = [transport_word(w, perm_q) for w in lx_raw]
            lz = [transport_word(w, perm_q) for w in lz_raw]
        else:
            lx, lz = lx_raw, lz_raw
        cand_meta = verify_candidate_logicals(
            code,
            lx,
            lz,
            signs_x=cl.get("signs_x"),
            signs_z=cl.get("signs_z"),
            slot_names=cl.get("slot_names"),
            convention=cl.get("convention", "positive_css"),
        )

    wit_meta = {}
    if d.get("witnesses") is not None:
        if not isinstance(d["witnesses"], dict):
            raise codes.SpecError("witnesses must be a dict keyed by d_x / d_z")
        for q, wentry in d["witnesses"].items():
            if q not in ("d_x", "d_z"):
                raise codes.SpecError(f"unsupported witness key {q!r}; expected d_x or d_z")
            wit_meta[q] = _parse_witness_entry(wentry, code, q, perm_q, default_index_base=default_ibase)

    import_record = {
        "schema": SCHEMA,
        "ancestry": "supplied",
        "producer": dict(producer),
        "code_type": "css",
        "field": "GF(2)",
        "index_base": default_ibase,
        "raw_sha256": raw_sha,
        "sha256": code.sha256(),
        "source_hashes": source_hashes,
        "coordinate_map": cmap_meta,
        "provenance": dict(d.get("provenance") or {}),
        "paper_claims": dict(d.get("paper_claims") or {}),
        "simulation_estimates": list(d.get("simulation_estimates") or []),
    }
    code.meta = {
        "family": "import",
        "import": import_record,
    }
    if cand_meta is not None:
        code.meta["candidate_logicals"] = cand_meta
    if wit_meta:
        code.meta["witnesses"] = wit_meta
    return code
