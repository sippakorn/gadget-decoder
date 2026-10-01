"""Session files, schema "gadgetlab.session/1": load, save, validation and inputs_hash (§5, E8)."""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from gadgetlab.code import Code, build_code
from gadgetlab.errors import InputError
from gadgetlab.geometry import L_MAX, L_MIN, Geometry
from gadgetlab.graphs import CUSTOM, PRESETS, Graph, GraphError, canonical_edges, default_layout, validate_cubic
from gadgetlab.spec import PRESET_NAMES, SpecError, bit_rows, parse_rows

SCHEMA = "gadgetlab.session/1"
SCOPES = ("global", "local_closed")
RESULT_KINDS = ("decode", "table", "search")
TOP_KEYS = ("schema", "graph", "L", "syndrome", "terminals", "scope", "patch", "region", "spec", "results", "notes")
MAX_TERMINALS = 6


class SessionError(InputError):
    pass


@dataclass
class Patch:
    vertices: list[int]  # ascending, distinct
    layers: tuple[int, int]

    def to_json(self) -> dict:
        return {"vertices": list(self.vertices), "layers": list(self.layers)}


@dataclass
class Session:
    graph: Graph
    L: int
    syndrome: list[int] = field(default_factory=list)  # Z-check indices with sigma = 1, ascending
    terminals: list[int] = field(default_factory=list)  # qubit indices, in terminal order
    scope: str = "global"
    patch: Patch | None = None
    region: list[int] = field(default_factory=list)  # Z-check indices, in region order
    spec: dict | None = None  # {"preset": name, "rows": {bit string: entry}}
    results: dict = field(default_factory=lambda: {k: None for k in RESULT_KINDS})
    notes: str = ""

    @property
    def code(self) -> Code:
        return build_code(self.graph.n, self.graph.edges, self.L)

    @property
    def geometry(self) -> Geometry:
        return self.code.geometry

    def sigma(self, checks=None) -> np.ndarray:
        s = np.zeros(self.geometry.m, dtype=np.uint8)
        s[list(self.syndrome if checks is None else checks)] = 1
        return s

    def patch_qubits(self) -> tuple[int, ...] | None:
        """Q for the current scope: None in global scope, the patch qubits in local_closed scope."""
        if self.scope == "global":
            return None
        if self.patch is None:
            raise InputError("local_closed scope needs a patch: choose patch vertices and a layer range")
        return tuple(self.geometry.patch(self.patch.vertices, *self.patch.layers))

    def spec_rows(self) -> dict | None:
        return None if self.spec is None else self.spec["rows"]

    def to_json(self) -> dict:
        geo = self.geometry
        g = self.graph
        return {
            "schema": SCHEMA,
            "graph": {
                "preset": g.preset,
                "seed": g.seed,
                "n": g.n,
                "edges": [[u, w] for u, w in g.edges],
                "layout": {str(v): [g.layout[v][0], g.layout[v][1]] for v in range(g.n)},
            },
            "L": self.L,
            "syndrome": [geo.check_json(i) for i in self.syndrome],
            "terminals": [geo.qubit_json(j) for j in self.terminals],
            "scope": self.scope,
            "patch": None if self.patch is None else self.patch.to_json(),
            "region": [geo.check_json(i) for i in self.region],
            "spec": copy.deepcopy(self.spec),
            "results": copy.deepcopy(self.results),
            "notes": self.notes,
        }


def new_session(graph: Graph, L: int) -> Session:
    build_code(graph.n, graph.edges, L)  # validates L and the construction
    return Session(graph=graph, L=L)


# ---------------------------------------------------------------- loading and validation


def load(path) -> Session:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError as exc:
        raise SessionError(f"cannot read {path}: {exc.strerror}") from None
    try:
        obj = json.loads(text)
    except json.JSONDecodeError as exc:
        raise SessionError(f"{path} is not valid JSON: {exc}") from None
    return from_json(obj)


def save(session: Session, path) -> None:
    Path(path).write_text(dumps(session.to_json()), encoding="utf-8")


def dumps(obj) -> str:
    return json.dumps(obj, indent=2, ensure_ascii=False) + "\n"


def _int(value, field: str, lo: int | None = None, hi: int | None = None) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise SessionError(f"{field}: expected an integer, got {value!r}")
    if hi is None and lo is not None and value < lo:
        raise SessionError(f"{field}: {value} is below the minimum {lo}")
    if (lo is not None and value < lo) or (hi is not None and value > hi):
        raise SessionError(f"{field}: {value} is out of range {lo}..{hi}")
    return value


def _list(value, field: str) -> list:
    if not isinstance(value, list):
        raise SessionError(f"{field}: expected a list, got {type(value).__name__}")
    return value


def _keys(obj, required, field: str, optional=()) -> None:
    if not isinstance(obj, dict):
        raise SessionError(f"{field}: expected an object, got {type(obj).__name__}")
    for k in required:
        if k not in obj:
            raise SessionError(f"{field}: missing key {k!r}")
    extra = sorted(set(obj) - set(required) - set(optional))
    if extra:
        raise SessionError(f"{field}: unexpected key {extra[0]!r}")


def parse_graph(obj, field: str = "graph") -> Graph:
    """The graph object of a session file: edges must already be canonical."""
    _keys(obj, ("n", "edges"), field, optional=("preset", "seed", "layout"))
    preset = obj.get("preset") or CUSTOM
    if preset not in PRESETS and preset != CUSTOM:
        raise SessionError(f"{field}.preset: unknown preset {preset!r}")
    seed = obj.get("seed")
    if seed is not None:
        _int(seed, f"{field}.seed", 0)
    n = _int(obj["n"], f"{field}.n", 4)
    edges_in = _list(obj["edges"], f"{field}.edges")
    try:
        edges = canonical_edges(edges_in)
        if [list(e) for e in edges] != [list(e) if isinstance(e, (list, tuple)) else e for e in edges_in]:
            raise GraphError("edges must be in canonical order: each edge [u, w] with u < w, the list sorted")
        validate_cubic(n, edges)
    except GraphError as exc:
        raise SessionError(f"{field}.edges: {exc}") from None
    layout_in = obj.get("layout")
    if layout_in is None:
        layout = default_layout(preset, n, edges)
    else:
        _keys(layout_in, [str(v) for v in range(n)], f"{field}.layout")
        layout = {}
        for v in range(n):
            point = layout_in[str(v)]
            if (not isinstance(point, list) or len(point) != 2
                    or not all(isinstance(c, (int, float)) and not isinstance(c, bool) and 0 <= c <= 1 for c in point)):
                raise SessionError(f"{field}.layout.{v}: expected [x, y] with 0 <= x, y <= 1, got {point!r}")
            layout[v] = (float(point[0]), float(point[1]))
    return Graph(n=n, edges=edges, preset=preset, seed=seed, layout=layout)


def from_json(obj) -> Session:
    """Validate a session object (§5); every problem raises SessionError naming the field."""
    try:
        return _from_json(obj)
    except SessionError:
        raise
    except InputError as exc:
        raise SessionError(str(exc)) from None


def _from_json(obj) -> Session:
    if not isinstance(obj, dict):
        raise SessionError("a session must be a JSON object")
    if obj.get("schema") != SCHEMA:
        raise SessionError(f'schema: expected "{SCHEMA}", got {obj.get("schema")!r}')
    _keys(obj, ("schema", "graph", "L"), "session", optional=TOP_KEYS)
    graph = parse_graph(obj["graph"])
    L = _int(obj["L"], "L", L_MIN, L_MAX)
    session = Session(graph=graph, L=L)
    geo = session.geometry

    def checks(name: str) -> list[int]:
        out = []
        for k, item in enumerate(_list(obj.get(name, []), name)):
            i = geo.parse_check(item, f"{name}[{k}]")
            if i in out:
                raise SessionError(f"{name}[{k}]: check (e = {item['e']}, t = {item['t']}) appears twice")
            out.append(i)
        return out

    session.syndrome = sorted(checks("syndrome"))
    session.region = checks("region")

    terminals = []
    for k, item in enumerate(_list(obj.get("terminals", []), "terminals")):
        j = geo.parse_qubit(item, f"terminals[{k}]")
        if j in terminals:
            raise SessionError(f"terminals[{k}]: this qubit is already terminal {terminals.index(j) + 1}")
        terminals.append(j)
    if len(terminals) > MAX_TERMINALS:
        raise SessionError(f"terminals: at most {MAX_TERMINALS} terminals, got {len(terminals)}")
    session.terminals = terminals

    scope = obj.get("scope", "global")
    if scope not in SCOPES:
        raise SessionError(f'scope: expected "global" or "local_closed", got {scope!r}')
    session.scope = scope

    patch = obj.get("patch")
    if patch is not None:
        _keys(patch, ("vertices", "layers"), "patch")
        vertices = []
        for k, v in enumerate(_list(patch["vertices"], "patch.vertices")):
            _int(v, f"patch.vertices[{k}]", 0, graph.n - 1)
            if v in vertices:
                raise SessionError(f"patch.vertices[{k}]: vertex {v} appears twice")
            vertices.append(v)
        layers = _list(patch["layers"], "patch.layers")
        if len(layers) != 2:
            raise SessionError("patch.layers: expected [l0, l1]")
        l0 = _int(layers[0], "patch.layers[0]", 0, L - 1)
        l1 = _int(layers[1], "patch.layers[1]", 0, L - 1)
        if l0 > l1:
            raise SessionError(f"patch.layers: need l0 <= l1, got [{l0}, {l1}]")
        session.patch = Patch(sorted(vertices), (l0, l1))

    spec = obj.get("spec")
    if spec is not None:
        _keys(spec, ("rows",), "spec", optional=("preset",))
        preset = spec.get("preset", "custom")
        if preset not in PRESET_NAMES:
            raise SessionError(f"spec.preset: expected one of {', '.join(PRESET_NAMES)}, got {preset!r}")
        try:
            rows = parse_rows(spec["rows"], len(terminals), "spec.rows")
        except SpecError as exc:
            raise SessionError(f"{exc} (the spec needs one row per assignment of the {len(terminals)} terminals)") from None
        session.spec = {"preset": preset, "rows": rows}

    results = obj.get("results")
    session.results = {k: None for k in RESULT_KINDS}
    if results is not None:
        _keys(results, (), "results", optional=RESULT_KINDS)
        for kind in RESULT_KINDS:
            value = results.get(kind)
            if value is not None:
                VALIDATORS[kind](geo, value, f"results.{kind}")
                session.results[kind] = copy.deepcopy(value)

    notes = obj.get("notes", "")
    if not isinstance(notes, str):
        raise SessionError("notes: expected a string")
    session.notes = notes
    return session


# ---------------------------------------------------------------- result validation

STATUSES = ("optimal", "infeasible")


def _cost(value, field: str) -> None:
    if value is not None:
        _int(value, field, 0)


def _error_list(geo: Geometry, value, field: str) -> None:
    if value is None:
        return
    for k, item in enumerate(_list(value, field)):
        geo.parse_qubit(item, f"{field}[{k}]")


def _validate_decode(geo: Geometry, obj, field: str) -> None:
    _keys(obj, ("status", "cost", "error", "verified", "solver", "ms", "inputs_hash"), field)
    if obj["status"] not in STATUSES:
        raise SessionError(f"{field}.status: expected one of {STATUSES}, got {obj['status']!r}")
    _cost(obj["cost"], f"{field}.cost")
    _error_list(geo, obj["error"], f"{field}.error")


def _validate_table(geo: Geometry, obj, field: str) -> None:
    _keys(obj, ("terminals", "rows", "unconstrained_cost", "match", "inputs_hash"), field)
    terminals = _list(obj["terminals"], f"{field}.terminals")
    for k, item in enumerate(terminals):
        geo.parse_qubit(item, f"{field}.terminals[{k}]")
    rows = _list(obj["rows"], f"{field}.rows")
    expected = bit_rows(len(terminals))
    if [r.get("a") if isinstance(r, dict) else None for r in rows] != expected:
        raise SessionError(f"{field}.rows: expected rows {expected} in order")
    for k, row in enumerate(rows):
        _keys(row, ("a", "status", "cost", "error"), f"{field}.rows[{k}]")
        if row["status"] not in STATUSES:
            raise SessionError(f"{field}.rows[{k}].status: expected one of {STATUSES}")
        _cost(row["cost"], f"{field}.rows[{k}].cost")
        _error_list(geo, row["error"], f"{field}.rows[{k}].error")
    _cost(obj["unconstrained_cost"], f"{field}.unconstrained_cost")
    if obj["match"] is not None:
        _keys(obj["match"], ("matches", "c", "margin"), f"{field}.match")


def _validate_search(geo: Geometry, obj, field: str) -> None:
    _keys(obj, ("params", "summary", "matches", "inputs_hash"), field)
    if not isinstance(obj["params"], dict):
        raise SessionError(f"{field}.params: expected an object")
    _keys(obj["summary"], ("candidates", "filtered", "infeasible_syndrome", "rejected", "matched",
                           "solver_calls", "cache_hits", "elapsed_s", "exhaustive"), f"{field}.summary")
    for k, item in enumerate(_list(obj["matches"], f"{field}.matches")):
        _keys(item, ("rank", "pattern", "c", "margin", "table"), f"{field}.matches[{k}]")
        _validate_table(geo, item["table"], f"{field}.matches[{k}].table")


VALIDATORS = {"decode": _validate_decode, "table": _validate_table, "search": _validate_search}


# ---------------------------------------------------------------- inputs_hash (V11)


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def hash_inputs(obj) -> str:
    """First 16 hex digits of SHA-256 over the canonical JSON; the viewer computes the same."""
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()[:16]


def _base_inputs(session: Session) -> dict:
    patch = None
    if session.scope == "local_closed" and session.patch is not None:
        patch = session.patch.to_json()
    return {
        "graph": {"n": session.graph.n, "edges": [[u, w] for u, w in session.graph.edges]},
        "L": session.L,
        "scope": session.scope,
        "patch": patch,
    }


def hash_params(params: dict) -> dict:
    return {k: params.get(k) for k in ("budget", "seed", "top_k", "max_weight")}


def inputs(kind: str, session: Session, params: dict | None = None, syndrome=None) -> dict:
    """Everything a result of this kind depends on, as plain JSON."""
    syndrome = sorted(session.syndrome if syndrome is None else syndrome)
    out = _base_inputs(session)
    out["kind"] = kind
    if kind == "decode":
        out["syndrome"] = syndrome
    elif kind == "table":
        out["syndrome"] = syndrome
        out["terminals"] = list(session.terminals)
        out["spec"] = session.spec_rows()
    elif kind == "search":
        region = set(session.region)
        out["syndrome_outside_region"] = [i for i in syndrome if i not in region]
        out["region"] = list(session.region)
        out["terminals"] = list(session.terminals)
        out["spec"] = session.spec_rows()
        out["params"] = hash_params(params or {})
    else:
        raise ValueError(f"unknown result kind {kind!r}")
    return out


def inputs_hash(kind: str, session: Session, params: dict | None = None, syndrome=None) -> str:
    return hash_inputs(inputs(kind, session, params, syndrome))
