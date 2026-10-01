"""Search over syndrome patterns on a region (§3.12): candidates, lazy evaluation, cache,
process pool, progress, cancel and jobs.

Works only with a check matrix and integer indices of qubits and checks.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import random
import signal
import threading
import time
import traceback
import uuid
from dataclasses import dataclass, replace

import numpy as np

from gadgetlab import f2
from gadgetlab.errors import InputError
from gadgetlab.solve import BACKENDS, DEFAULT_BACKEND, Decoder
from gadgetlab.spec import bit_rows, match as spec_match, parse_rows, problems as spec_problems
from gadgetlab.table import CostTable, assignment, check_invariants, validate_terminals

DEFAULT_BUDGET = 4096
DEFAULT_SEED = 0
DEFAULT_TOP_K = 50
MAX_BUDGET = 1_000_000
MAX_WORKERS = 64


def default_workers() -> int:
    return max(1, (os.cpu_count() or 1) - 1)


class SearchCancelled(Exception):
    def __init__(self, summary: dict):
        super().__init__("search cancelled")
        self.summary = summary


def _check_int(name: str, value, lo: int, hi: int | None = None) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < lo or (hi is not None and value > hi):
        bound = f"between {lo} and {hi}" if hi is not None else f"at least {lo}"
        raise InputError(f"{name} must be an integer {bound}, got {value!r}")
    return value


@dataclass(frozen=True)
class SearchParams:
    budget: int = DEFAULT_BUDGET
    seed: int = DEFAULT_SEED
    top_k: int = DEFAULT_TOP_K
    max_weight: int | None = None
    workers: int | None = None  # None means default_workers()
    backend: str = DEFAULT_BACKEND

    def resolved(self) -> SearchParams:
        _check_int("budget", self.budget, 1, MAX_BUDGET)
        _check_int("seed", self.seed, 0)
        _check_int("top_k", self.top_k, 1)
        if self.max_weight is not None:
            _check_int("max_weight", self.max_weight, 0)
        workers = default_workers() if self.workers is None else _check_int("workers", self.workers, 1, MAX_WORKERS)
        if self.backend not in BACKENDS:
            raise InputError(f"unknown solver {self.backend!r}; expected one of {', '.join(BACKENDS)}")
        return replace(self, workers=workers)


@dataclass(frozen=True)
class SearchProblem:
    H: np.ndarray
    sigma0: np.ndarray
    region: tuple[int, ...]  # ordered Z-check indices
    terminals: tuple[int, ...]
    rows: dict  # spec rows, bit string -> entry
    Q: tuple[int, ...] | None  # None means global scope


@dataclass
class MatchRecord:
    pattern: str  # aligned with the full region
    weight: int  # number of ones among the varied region positions
    c: int
    margin: int | None  # None means infinity
    table: CostTable


def rank_key(rec: MatchRecord):
    """Margin descending (infinity first), then c ascending, then |p| ascending, then p lexicographically."""
    return (0 if rec.margin is None else 1, 0 if rec.margin is None else -rec.margin, rec.c, rec.weight, rec.pattern)


@dataclass
class SearchOutcome:
    params: SearchParams
    summary: dict
    matches: list[MatchRecord]  # ranked, at most top_k
    dropped: list[int]  # region positions dropped because they lie outside C_Q


class _Stop(Exception):
    pass


def candidate_patterns(R: int, budget: int, seed: int, should_cancel=None) -> tuple[list[int], bool]:
    """Patterns as integers whose bits, most significant first, are p_1 ... p_R.

    Every pattern in lexicographic order when 2^R <= budget (exhaustive); otherwise budget
    distinct patterns drawn uniformly with random.Random(seed), returned in lexicographic order.
    When the budget exceeds half of the 2^R patterns, the excluded patterns are drawn instead.
    """
    total = 1 << R
    if total <= budget:
        return list(range(total)), True
    rng = random.Random(seed)
    k = budget if 2 * budget <= total else total - budget
    chosen: set[int] = set()
    draws = 0
    while len(chosen) < k:
        chosen.add(rng.getrandbits(R))
        draws += 1
        if should_cancel is not None and draws % 65536 == 0 and should_cancel():
            raise _Stop
    if k == budget:
        return sorted(chosen), False
    return [p for p in range(total) if p not in chosen], False


class Evaluator:
    """Lazy evaluation of one candidate pattern against the spec (§3.12)."""

    def __init__(self, problem: SearchProblem, backend: str, check: bool):
        self.decoder = Decoder(problem.H, problem.Q, backend, check)
        self.terminals = validate_terminals(self.decoder, problem.terminals)
        k = len(self.terminals)
        self.rows = parse_rows(problem.rows, k, "spec.rows")
        issues = spec_problems(self.rows)
        if issues:
            raise InputError(f"the spec cannot be matched: {issues[0]}")
        order = bit_rows(k)
        self.order = order
        self.eq_rows = [a for a in order if self.rows[a][0] == "eq"]
        self.ge_rows = [a for a in order if self.rows[a][0] == "ge"]
        self.any_rows = [a for a in order if self.rows[a][0] == "any"]
        self.sigma0 = f2.as_f2(problem.sigma0).reshape(-1).copy()
        if self.sigma0.shape[0] != self.decoder.m:
            raise InputError(f"the base syndrome has length {self.sigma0.shape[0]}, expected {self.decoder.m}")
        region = [int(i) for i in problem.region]
        if not region:
            raise InputError("the region is empty: add at least one check")
        if len(set(region)) != len(region):
            raise InputError("the region contains a check twice")
        for i in region:
            if not 0 <= i < self.decoder.m:
                raise InputError(f"region check {i} is out of range 0..{self.decoder.m - 1}")
        self.region = region
        self.varied = [pos for pos, i in enumerate(region) if i in self.decoder.enforced_set]
        self.dropped = [pos for pos, i in enumerate(region) if i not in self.decoder.enforced_set]
        self.R = len(self.varied)
        self.colspace = f2.ColumnSpaceTest(problem.H) if problem.Q is None else None

    def sigma_for(self, p: int) -> np.ndarray:
        sigma = self.sigma0.copy()
        for k, pos in enumerate(self.varied):
            sigma[self.region[pos]] = (p >> (self.R - 1 - k)) & 1
        return sigma

    def pattern_string(self, p: int) -> str:
        sigma = self.sigma_for(p)
        return "".join(str(int(sigma[i])) for i in self.region)

    def evaluate(self, p: int) -> tuple:
        """(p, "rejected" | "matched", MatchRecord or None, solver calls, cache hits)."""
        calls, hits = self.decoder.solver_calls, self.decoder.cache_hits
        record = self._evaluate(p)
        kind = "rejected" if record is None else "matched"
        return (p, kind, record, self.decoder.solver_calls - calls, self.decoder.cache_hits - hits)

    def _evaluate(self, p: int) -> MatchRecord | None:
        sigma = self.sigma_for(p)
        results = {}

        def solve(a):
            results[a] = self.decoder.solve(sigma, assignment(self.terminals, a))
            return results[a]

        c = None
        for a in self.eq_rows:  # fix c, rejecting at the first inconsistency
            r = solve(a)
            if r.status != "optimal":
                return None
            if c is None:
                c = r.cost - self.rows[a][1]
            elif r.cost - self.rows[a][1] != c:
                return None
        for a in self.ge_rows:  # then reject at the first violation
            r = solve(a)
            if r.status == "optimal" and r.cost < c + self.rows[a][1]:
                return None
        for a in self.any_rows:  # complete the table of a match
            solve(a)
        table = CostTable(self.terminals, {a: results[a] for a in self.order}, self.decoder.solve(sigma, ()))
        check_invariants(table)
        m = spec_match(self.rows, table.costs())
        if not m.matches or m.c != c:
            raise RuntimeError("lazy evaluation and full matching disagree")
        return MatchRecord(self.pattern_string(p), bin(p).count("1"), m.c, m.margin, table)


def run_search(problem: SearchProblem, params: SearchParams, *, check: bool = True,
               progress=None, should_cancel=None) -> SearchOutcome:
    """Run a search. progress(done, total, summary) is called as candidates finish; should_cancel()
    is polled between candidates, and a cancelled search raises SearchCancelled."""
    params = params.resolved()
    start = time.perf_counter()
    ev = Evaluator(problem, params.backend, check)
    counts = {"candidates": 0, "filtered": 0, "infeasible_syndrome": 0, "rejected": 0,
              "matched": 0, "solver_calls": 0, "cache_hits": 0}
    exhaustive = (1 << ev.R) <= params.budget

    def summary() -> dict:
        return {**counts, "elapsed_s": round(time.perf_counter() - start, 3), "exhaustive": exhaustive}

    try:
        patterns, exhaustive = candidate_patterns(ev.R, params.budget, params.seed, should_cancel)
    except _Stop:
        raise SearchCancelled(summary()) from None
    counts["candidates"] = len(patterns)
    test_colspace = ev.colspace is not None and not ev.colspace.full_row_rank
    todo = []
    for k, p in enumerate(patterns):
        if k % 16384 == 0 and should_cancel is not None and should_cancel():
            raise SearchCancelled(summary())
        if params.max_weight is not None and bin(p).count("1") > params.max_weight:
            counts["filtered"] += 1
        elif test_colspace and not ev.colspace.contains(ev.sigma_for(p)):
            counts["infeasible_syndrome"] += 1
        else:
            todo.append(p)
    matches: list[MatchRecord] = []

    def report() -> None:
        if progress is not None:
            done = sum(counts[k] for k in ("filtered", "infeasible_syndrome", "rejected", "matched"))
            progress(done, counts["candidates"], summary())

    def absorb(outcome) -> None:
        _p, kind, record, calls, hits = outcome
        counts[kind] += 1
        counts["solver_calls"] += calls
        counts["cache_hits"] += hits
        if record is not None:
            matches.append(record)
            if len(matches) > 2 * params.top_k + 64:
                matches.sort(key=rank_key)
                del matches[params.top_k:]
        report()

    report()
    cancelled = False
    if params.workers <= 1 or len(todo) <= 1:
        for p in todo:
            if should_cancel is not None and should_cancel():
                cancelled = True
                break
            try:
                outcome = ev.evaluate(p)
            except Exception:
                if should_cancel is not None and should_cancel():  # e.g. Ctrl+C also reached the solver
                    cancelled = True
                    break
                raise
            absorb(outcome)
    else:
        cancelled = _run_pool(problem, params, check, todo, absorb, should_cancel)
    if cancelled:
        raise SearchCancelled(summary())
    matches.sort(key=rank_key)
    return SearchOutcome(params, summary(), matches[: params.top_k], ev.dropped)


# ---------------------------------------------------------------- process pool

_EVALUATOR: Evaluator | None = None
_CANCEL = None


def _worker_init(problem: SearchProblem, backend: str, check: bool, cancel_event) -> None:
    global _EVALUATOR, _CANCEL
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # the parent handles Ctrl+C
    _EVALUATOR = Evaluator(problem, backend, check)
    _CANCEL = cancel_event


def _worker_evaluate(chunk: list[int]) -> list[tuple]:
    out = []
    for p in chunk:  # the cancel flag is checked between candidates
        if _CANCEL.is_set():
            break
        out.append(_EVALUATOR.evaluate(p))
    return out


def _run_pool(problem, params, check, todo, absorb, should_cancel) -> bool:
    """Evaluate todo in a spawn-context process pool; returns True when cancelled."""
    ctx = mp.get_context("spawn")
    cancel_event = ctx.Event()
    workers = min(params.workers, len(todo))
    size = max(1, min(8, len(todo) // (workers * 16)))
    chunks = [todo[k:k + size] for k in range(0, len(todo), size)]
    pool = ctx.Pool(workers, initializer=_worker_init, initargs=(problem, params.backend, check, cancel_event))
    try:
        it = pool.imap_unordered(_worker_evaluate, chunks, chunksize=1)
        remaining = len(chunks)
        while remaining:
            if should_cancel is not None and should_cancel():
                cancel_event.set()
                return True
            try:
                outcomes = it.next(timeout=0.1)
            except mp.TimeoutError:
                continue
            except Exception:
                if should_cancel is not None and should_cancel():  # e.g. Ctrl+C also reached a solver
                    cancel_event.set()
                    return True
                raise
            remaining -= 1
            for outcome in outcomes:
                absorb(outcome)
        return False
    finally:
        pool.terminate()
        pool.join()


# ---------------------------------------------------------------- jobs


class Job:
    """A search running in a background thread, with progress, cancel and a JSON snapshot."""

    def __init__(self, target, debug: bool = False):
        self.id = uuid.uuid4().hex[:12]
        self.state = "running"
        self.progress = {"done": 0, "total": 0}
        self.summary = None
        self.result = None
        self.error = None
        self.traceback = None
        self.debug = debug
        self.created = time.time()
        self._target = target
        self._cancel = threading.Event()
        self._lock = threading.Lock()
        self._thread = threading.Thread(target=self._run, name=f"search-{self.id}", daemon=True)

    def start(self) -> Job:
        self._thread.start()
        return self

    def _on_progress(self, done: int, total: int, summary: dict) -> None:
        with self._lock:
            self.progress = {"done": done, "total": total}
            self.summary = summary

    def _run(self) -> None:
        try:
            result = self._target(self._on_progress, self._cancel.is_set)
        except SearchCancelled as exc:
            with self._lock:
                self.state = "cancelled"
                self.summary = exc.summary
        except Exception as exc:  # reported through the job, never raised into the server
            with self._lock:
                if self.state != "cancelled":
                    self.state = "error"
                    self.error = str(exc) or type(exc).__name__
                    self.traceback = traceback.format_exc()
        else:
            with self._lock:
                if self.state == "running":
                    self.state = "done"
                    self.result = result
                    self.summary = result.get("summary", self.summary)

    def cancel(self) -> str:
        self._cancel.set()
        with self._lock:
            if self.state == "running":
                self.state = "cancelled"
            return self.state

    def finished(self, timeout: float | None = None) -> bool:
        self._thread.join(timeout)
        return not self._thread.is_alive()

    def snapshot(self) -> dict:
        with self._lock:
            out = {"state": self.state, "progress": dict(self.progress), "summary": self.summary, "result": self.result}
            if self.state == "error":
                out["error"] = self.error
                if self.debug:
                    out["traceback"] = self.traceback
            return out
