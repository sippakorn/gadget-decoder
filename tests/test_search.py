"""T8 search: candidate counts, sampling, ranking, cancel, lazy = full evaluation; cache and filters."""

import time

import numpy as np
import pytest

from gadgetlab import code as codemod
from gadgetlab.graphs import from_preset
from gadgetlab.search import (Job, MatchRecord, SearchCancelled, SearchParams, SearchProblem, candidate_patterns,
                              rank_key, run_search)
from gadgetlab.solve import Decoder
from gadgetlab.spec import bit_rows, match, preset_rows
from gadgetlab.table import compute_table


def _code(preset, L):
    g = from_preset(preset)
    return codemod.build_code(g.n, g.edges, L)


def _problem(c, region, terminals, rows, Q=None, sigma_on=()):
    sigma0 = np.zeros(c.m, dtype=np.uint8)
    sigma0[list(sigma_on)] = 1
    return SearchProblem(c.HZ, sigma0, tuple(region), tuple(terminals), rows, None if Q is None else tuple(Q))


def test_exhaustive_candidate_count_is_two_to_the_region_size():
    c = _code("k4", 2)
    geo = c.geometry
    problem = _problem(c, [0, 2, 4], [geo.edge_qubit(0, 0), geo.edge_qubit(5, 1)], preset_rows("wire", 2))
    out = run_search(problem, SearchParams(workers=1, backend="brute"))
    s = out.summary
    assert s["candidates"] == 8 and s["exhaustive"] is True
    assert s["filtered"] + s["infeasible_syndrome"] + s["rejected"] + s["matched"] == 8
    assert s["infeasible_syndrome"] == 0
    assert s["solver_calls"] >= 2 * 8


def test_sampled_patterns_are_distinct_and_reproducible():
    a, exhaustive = candidate_patterns(12, 100, seed=4)
    b, _ = candidate_patterns(12, 100, seed=4)
    c, _ = candidate_patterns(12, 100, seed=5)
    assert not exhaustive and len(a) == 100 and len(set(a)) == 100
    assert all(0 <= p < 2 ** 12 for p in a)
    assert a == b and a != c
    big, _ = candidate_patterns(80, 50, seed=1)
    assert len(set(big)) == 50 and max(big) < 2 ** 80
    assert candidate_patterns(3, 8, 0) == (list(range(8)), True)


def test_sampled_search_summary():
    c = _code("petersen", 2)
    geo = c.geometry
    problem = _problem(c, list(range(12)), [geo.edge_qubit(3, 0), geo.vertex_qubit(2, 0)], preset_rows("wire", 2))
    out = run_search(problem, SearchParams(budget=40, seed=2, workers=1))
    s = out.summary
    assert s["candidates"] == 40 and s["exhaustive"] is False
    assert s["filtered"] + s["infeasible_syndrome"] + s["rejected"] + s["matched"] == 40


def test_ranking_follows_section_3_12():
    def rec(pattern, c, margin):
        return MatchRecord(pattern, pattern.count("1"), c, margin, None)

    records = [rec("0110", 2, 1), rec("0001", 3, None), rec("1000", 2, 1), rec("0011", 2, 3),
               rec("0100", 1, 1), rec("1111", 0, None), rec("0010", 2, 1)]
    order = [r.pattern for r in sorted(records, key=rank_key)]
    # infinite margin first (c ascending), then margin descending, then c, then |p|, then p
    assert order == ["1111", "0001", "0011", "0100", "0010", "1000", "0110"]


def _full_matches(c, problem):
    """Reference: compute every row of every candidate, then match."""
    decoder = Decoder(c.HZ, problem.Q, "brute")
    out = set()
    R = len(problem.region)
    for p in range(2 ** R):
        sigma = problem.sigma0.copy()
        for k, i in enumerate(problem.region):
            sigma[i] = (p >> (R - 1 - k)) & 1
        table = compute_table(decoder, sigma, problem.terminals)
        m = match(problem.rows, table.costs())
        if m.matches:
            out.add((format(p, f"0{R}b"), m.c, m.margin))
    return out


@pytest.mark.parametrize("rows", [
    preset_rows("wire", 2),
    {"00": ["eq", 0], "01": ["any"], "10": ["ge", 1], "11": ["ge", 0]},
    {"00": ["ge", 0], "01": ["eq", 1], "10": ["eq", 1], "11": ["any"]},
])
def test_lazy_evaluation_finds_the_same_matches_as_full_evaluation(rows):
    c = _code("k4", 3)
    geo = c.geometry
    region = [geo.z_check(e, t) for e, t in ((0, 0), (1, 0), (2, 1), (3, 1), (4, 0), (5, 1), (3, 0))]
    terminals = [geo.edge_qubit(0, 1), geo.vertex_qubit(3, 1)]
    problem = _problem(c, region, terminals, rows, sigma_on=[geo.z_check(5, 0)])
    out = run_search(problem, SearchParams(top_k=1000, workers=1, backend="brute"))
    lazy = {(m.pattern, m.c, m.margin) for m in out.matches}
    assert lazy == _full_matches(c, problem)
    assert out.summary["matched"] == len(lazy)
    keys = [rank_key(m) for m in out.matches]
    assert keys == sorted(keys)
    for m in out.matches:
        assert list(m.table.rows) == bit_rows(2)  # complete tables, one witness per feasible row
        assert all(r.witness is not None for r in m.table.rows.values() if r.status == "optimal")


def test_pool_gives_the_same_matches_as_one_worker():
    c = _code("k4", 3)
    geo = c.geometry
    region = [geo.z_check(e, t) for e, t in ((0, 0), (1, 0), (2, 1), (3, 1), (4, 0), (5, 1))]
    rows = {"00": ["eq", 0], "01": ["any"], "10": ["ge", 1], "11": ["ge", 0]}
    problem = _problem(c, region, [geo.edge_qubit(0, 1), geo.vertex_qubit(3, 1)], rows)
    one = run_search(problem, SearchParams(workers=1))
    two = run_search(problem, SearchParams(workers=3))
    assert [(m.pattern, m.c, m.margin) for m in one.matches] == [(m.pattern, m.c, m.margin) for m in two.matches]
    assert one.summary["solver_calls"] == two.summary["solver_calls"]
    for a, b in zip(one.matches, two.matches):
        assert all(np.array_equal(a.table.rows[k].witness, b.table.rows[k].witness)
                   for k in a.table.rows if a.table.rows[k].witness is not None)


def test_filters_max_weight():
    c = _code("k4", 2)
    geo = c.geometry
    problem = _problem(c, [0, 1, 2, 3], [geo.edge_qubit(1, 0)], {"0": ["eq", 0], "1": ["ge", 0]})
    out = run_search(problem, SearchParams(max_weight=1, workers=1, backend="brute"))
    assert out.summary["candidates"] == 16 and out.summary["filtered"] == 11
    assert all(m.weight <= 1 for m in out.matches)


def test_local_closed_drops_region_checks_outside_the_patch():
    c = _code("k4", 3)
    geo = c.geometry
    Q = geo.patch([0], 0, 0)  # edge-qubits (E, e, 0) at vertex 0; C_Q = the checks (e, 0) at vertex 0
    at0 = geo.edges_at[0]
    region = [geo.z_check(at0[0], 0), geo.z_check(at0[0], 1), geo.z_check(at0[1], 0)]
    sigma_on = [geo.z_check(at0[0], 1)]
    problem = _problem(c, region, [geo.edge_qubit(at0[0], 0)], {"0": ["eq", 0], "1": ["any"]}, Q=Q, sigma_on=sigma_on)
    out = run_search(problem, SearchParams(workers=1, backend="brute", top_k=10))
    assert out.dropped == [1]
    assert out.summary["candidates"] == 4
    assert all(m.pattern[1] == "1" for m in out.matches)  # the dropped check keeps its base value


def test_cache_key_ignores_checks_that_are_not_enforced():
    c = _code("k4", 3)
    geo = c.geometry
    Q = geo.patch([0], 0, 0)
    decoder = Decoder(c.HZ, Q, "ilp_cbc")
    sigma = np.zeros(c.m, dtype=np.uint8)
    other = sigma.copy()
    outside = next(i for i in range(c.m) if i not in decoder.enforced_set)
    other[outside] = 1
    decoder.solve(sigma)
    decoder.solve(other)
    decoder.solve(sigma, {Q[0]: 1})
    assert (decoder.solver_calls, decoder.cache_hits) == (2, 1)


def _slow_problem():
    c = _code("petersen", 3)
    geo = c.geometry
    region = list(range(0, 30, 2))
    return _problem(c, region, [geo.edge_qubit(3, 1), geo.vertex_qubit(2, 0)], preset_rows("wire", 2))


@pytest.mark.parametrize("workers", [1, 2])
def test_cancel_stops_a_job_within_two_seconds(workers):
    problem = _slow_problem()

    def target(progress, should_cancel):
        run_search(problem, SearchParams(budget=4096, workers=workers), progress=progress, should_cancel=should_cancel)
        return {"summary": {}}

    job = Job(target).start()
    deadline = time.monotonic() + 30
    while job.progress["done"] < 3 and time.monotonic() < deadline:
        time.sleep(0.05)
    assert job.progress["done"] >= 3, "the search made no progress"
    t0 = time.monotonic()
    assert job.cancel() == "cancelled"
    assert job.finished(timeout=2.0), "the job did not stop within 2 seconds"
    assert time.monotonic() - t0 < 2.0
    snap = job.snapshot()
    assert snap["state"] == "cancelled" and snap["result"] is None
    assert 0 < snap["progress"]["done"] < snap["progress"]["total"] == 4096


def test_cancel_flag_raises_search_cancelled():
    with pytest.raises(SearchCancelled) as info:
        run_search(_slow_problem(), SearchParams(budget=64, workers=1), should_cancel=lambda: True)
    assert info.value.summary["candidates"] == 64
