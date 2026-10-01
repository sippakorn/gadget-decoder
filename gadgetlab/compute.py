"""Session-level computations shared by the CLI and the server, so both produce identical JSON.

Translates between typed addresses (through the session's geometry) and the index-only
engine modules solve, table, spec and search.
"""

from __future__ import annotations

import random

import numpy as np

from gadgetlab import distance
from gadgetlab.code import code_info
from gadgetlab.errors import InputError
from gadgetlab.graphs import CUSTOM, NAMED_PRESETS, Graph, from_edges, from_preset
from gadgetlab.search import SearchOutcome, SearchParams, SearchProblem, run_search
from gadgetlab.session import Session, inputs_hash, new_session, parse_graph
from gadgetlab.solve import DEFAULT_BACKEND, Decoder, Result
from gadgetlab.spec import bit_rows, match as spec_match, problems as spec_problems
from gadgetlab.table import CostTable, compute_table

TIMING_KEYS = ("ms", "elapsed_s")


def info(session: Session) -> dict:
    return code_info(session.graph, session.L, session.code)


def build(graph_obj, L) -> tuple[Session, dict]:
    """A new session from {"preset", "n", "seed"} or an edge list (or a session's graph object)."""
    graph = graph_from_request(graph_obj)
    session = new_session(graph, L)
    return session, info(session)


def graph_from_request(obj) -> Graph:
    if not isinstance(obj, dict):
        raise InputError('graph: expected an object such as {"preset": "petersen"}')
    if obj.get("edges") is not None:
        if "layout" in obj or obj.get("preset") not in (None, CUSTOM):
            return parse_graph({k: obj[k] for k in ("preset", "seed", "n", "edges", "layout") if k in obj})
        return from_edges(obj["edges"], obj.get("n"))
    preset = obj.get("preset")
    if preset is None:
        raise InputError('graph: give a preset, e.g. {"preset": "petersen"}, or an edge list')
    if preset in NAMED_PRESETS:
        return from_preset(preset)
    return from_preset(preset, n=obj.get("n"), seed=obj.get("seed"))


def decoder_for(session: Session, backend: str = DEFAULT_BACKEND, debug: bool = True) -> Decoder:
    return Decoder(session.code.HZ, session.patch_qubits(), backend, check=debug)


def _error_json(session: Session, r: Result):
    geo = session.geometry
    return [geo.qubit_json(j) for j in r.ones()] if r.status == "optimal" else None


def decode(session: Session, backend: str = DEFAULT_BACKEND, debug: bool = True) -> dict:
    r = decoder_for(session, backend, debug).solve(session.sigma(), ())
    return {
        "status": r.status,
        "cost": r.cost,
        "error": _error_json(session, r),
        "verified": True if (debug and r.status == "optimal") else None,
        "solver": r.solver,
        "ms": round(r.ms, 1),
        "inputs_hash": inputs_hash("decode", session),
    }


def usable_spec_rows(session: Session) -> dict | None:
    rows = session.spec_rows()
    if rows is None or set(rows) != set(bit_rows(len(session.terminals))) or spec_problems(rows):
        return None
    return rows


def table_json(session: Session, table: CostTable, syndrome=None) -> dict:
    geo = session.geometry
    rows = usable_spec_rows(session)
    return {
        "terminals": [geo.qubit_json(t) for t in table.terminals],
        "rows": [{"a": a, "status": r.status, "cost": r.cost, "error": _error_json(session, r)}
                 for a, r in table.rows.items()],
        "unconstrained_cost": table.unconstrained.cost,
        "match": None if rows is None else spec_match(rows, table.costs()).to_json(),
        "inputs_hash": inputs_hash("table", session, syndrome=syndrome),
    }


def table(session: Session, backend: str = DEFAULT_BACKEND, debug: bool = True) -> dict:
    if not session.terminals:
        raise InputError("choose at least one terminal before computing a cost table")
    t = compute_table(decoder_for(session, backend, debug), session.sigma(), session.terminals)
    return table_json(session, t)


def distances(session: Session) -> dict:
    code = session.code
    return distance.distances(code.HX, code.HZ)


def random_syndrome(session: Session, weight, seed) -> dict:
    """{"syndrome"}: H_Z x for a uniformly random x of the given weight over all qubits."""
    geo = session.geometry
    if not isinstance(weight, int) or isinstance(weight, bool) or not 0 <= weight <= geo.N:
        raise InputError(f"weight must be an integer between 0 and {geo.N}, got {weight!r}")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise InputError(f"seed must be a non-negative integer, got {seed!r}")
    x = np.zeros(geo.N, dtype=np.int64)
    x[random.Random(seed).sample(range(geo.N), weight)] = 1
    sigma = (session.code.HZ.astype(np.int64) @ x) % 2
    return {"syndrome": [geo.check_json(int(i)) for i in np.flatnonzero(sigma)]}


# ---------------------------------------------------------------- search


def search_params(obj: dict | None) -> SearchParams:
    """SearchParams from {budget, seed, top_k, max_weight, workers, solver}; missing keys take defaults."""
    obj = obj or {}
    defaults = SearchParams()

    def pick(key, default):
        value = obj.get(key)
        return default if value is None else value

    return SearchParams(
        budget=pick("budget", defaults.budget),
        seed=pick("seed", defaults.seed),
        top_k=pick("top_k", defaults.top_k),
        max_weight=obj.get("max_weight"),
        workers=obj.get("workers"),
        backend=pick("solver", defaults.backend),
    )


def search_problem(session: Session) -> SearchProblem:
    if not session.terminals:
        raise InputError("choose at least one terminal before searching")
    if session.spec is None:
        raise InputError("write a spec before searching")
    if not session.region:
        raise InputError("the region is empty: add at least one check")
    return SearchProblem(
        H=session.code.HZ,
        sigma0=session.sigma(),
        region=tuple(session.region),
        terminals=tuple(session.terminals),
        rows=session.spec_rows(),
        Q=session.patch_qubits(),
    )


def search(session: Session, params: SearchParams, debug: bool = True, progress=None, should_cancel=None) -> dict:
    problem = search_problem(session)
    outcome = run_search(problem, params, check=debug, progress=progress, should_cancel=should_cancel)
    return search_json(session, outcome)


def search_json(session: Session, outcome: SearchOutcome) -> dict:
    geo = session.geometry
    p = outcome.params
    params = {
        "budget": p.budget,
        "seed": p.seed,
        "top_k": p.top_k,
        "max_weight": p.max_weight,
        "workers": p.workers,
        "solver": p.backend,
        "dropped_region": [geo.check_json(session.region[pos]) for pos in outcome.dropped],
    }
    matches = []
    for rank, rec in enumerate(outcome.matches, start=1):
        syndrome = pattern_syndrome(session, rec.pattern)
        matches.append({
            "rank": rank,
            "pattern": rec.pattern,
            "c": rec.c,
            "margin": rec.margin,
            "table": table_json(session, rec.table, syndrome=syndrome),
        })
    return {"params": params, "summary": outcome.summary, "matches": matches,
            "inputs_hash": inputs_hash("search", session, params)}


def pattern_syndrome(session: Session, pattern: str) -> list[int]:
    """The base syndrome with its entries on the region replaced by the pattern, as check indices."""
    if len(pattern) != len(session.region):
        raise InputError(f"pattern {pattern!r} does not match a region of {len(session.region)} checks")
    on = set(session.syndrome)
    for i, bit in zip(session.region, pattern):
        (on.add if bit == "1" else on.discard)(i)
    return sorted(on)


# ---------------------------------------------------------------- comparison


def strip_timing(obj):
    """A copy of a result with the timing fields (ms, elapsed_s) removed, for comparisons."""
    if isinstance(obj, dict):
        return {k: strip_timing(v) for k, v in obj.items() if k not in TIMING_KEYS}
    if isinstance(obj, list):
        return [strip_timing(v) for v in obj]
    return obj


def differences(a, b, path: str = "") -> list[str]:
    """Paths where two JSON values differ, ignoring timing fields."""
    a, b = strip_timing(a), strip_timing(b)
    out: list[str] = []
    _diff(a, b, path or "result", out)
    return out


def _diff(a, b, path: str, out: list[str]) -> None:
    if len(out) >= 20:
        return
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{path}.{k}: present in only one")
            else:
                _diff(a[k], b[k], f"{path}.{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: lengths {len(a)} and {len(b)}")
        for k, (x, y) in enumerate(zip(a, b)):
            _diff(x, y, f"{path}[{k}]", out)
    elif a != b or type(a) is not type(b):
        out.append(f"{path}: {a!r} != {b!r}")
