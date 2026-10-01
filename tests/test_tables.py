"""T6 table invariants on random cases, and T11 infeasibility."""

import math

import numpy as np
import pytest

from gadgetlab import code as codemod
from gadgetlab import f2
from gadgetlab.graphs import NAMED_PRESETS, ball, from_preset
from gadgetlab.solve import BACKENDS, Decoder, Result, highs_available
from gadgetlab.table import CostTable, InvariantError, check_invariants, compute_table
from gadgetlab.errors import InputError


def _instance(preset, L):
    g = from_preset(preset)
    return g, codemod.build_code(g.n, g.edges, L)


@pytest.mark.parametrize("preset,L,backend", [("k4", 2, "brute"), ("k4", 3, "ilp_cbc"), ("k33", 2, "brute"),
                                               ("petersen", 2, "ilp_cbc"), ("prism", 3, "ilp_cbc")])
def test_t6_invariants_hold_on_random_tables(preset, L, backend):
    g, c = _instance(preset, L)
    geo = c.geometry
    rng = np.random.default_rng(len(preset) * 7 + L)
    for case in range(12):
        sigma = rng.integers(0, 2, c.m).astype(np.uint8)
        if case % 2:
            U = ball(g.n, g.edges, int(rng.integers(g.n)), 1)
            l0 = int(rng.integers(0, L))
            Q = geo.patch(U, l0, int(rng.integers(l0, L)))
            pool = np.array(Q)
        else:
            Q = None
            pool = np.arange(c.N)
        k = int(rng.integers(1, 4))
        terminals = [int(t) for t in rng.choice(pool, k, replace=False)]
        decoder = Decoder(c.HZ, Q, backend)
        table = compute_table(decoder, sigma, terminals)  # raises InvariantError on a violation
        costs = [math.inf if r.cost is None else r.cost for r in table.rows.values()]
        free = math.inf if table.unconstrained.cost is None else table.unconstrained.cost
        assert min(costs) == free
        assert list(table.rows) == [format(i, f"0{k}b") for i in range(2 ** k)]
        if table.unconstrained.status == "optimal":
            a = "".join(str(int(table.unconstrained.witness[t])) for t in terminals)
            assert table.rows[a].cost == free


def _fake(cost, witness=None):
    if cost is None:
        return Result("infeasible", None, None, "fake", 0.0)
    return Result("optimal", cost, witness, "fake", 0.0)


def test_t6_invariant_violations_are_errors():
    w = np.zeros(4, dtype=np.uint8)
    w[0] = 1
    good = CostTable((0,), {"0": _fake(3, w), "1": _fake(1, w)}, _fake(1, w))
    check_invariants(good)
    with pytest.raises(InvariantError, match=r"\(i\)"):
        check_invariants(CostTable((0,), {"0": _fake(3, w), "1": _fake(2, w)}, _fake(1, w)))
    with pytest.raises(InvariantError, match=r"\(i\)"):
        check_invariants(CostTable((0,), {"0": _fake(None), "1": _fake(None)}, _fake(1, w)))
    w0 = np.zeros(4, dtype=np.uint8)  # the unconstrained witness has terminal value 0
    with pytest.raises(InvariantError, match=r"\(ii\)"):
        check_invariants(CostTable((0,), {"0": _fake(3, w0), "1": _fake(1, w)}, _fake(1, w0)))


def test_table_input_errors():
    g, c = _instance("k4", 2)
    geo = c.geometry
    sigma = np.zeros(c.m, dtype=np.uint8)
    with pytest.raises(InputError, match="between 1 and 6"):
        compute_table(Decoder(c.HZ), sigma, [])
    with pytest.raises(InputError, match="between 1 and 6"):
        compute_table(Decoder(c.HZ), sigma, list(range(7)))
    with pytest.raises(InputError, match="distinct"):
        compute_table(Decoder(c.HZ), sigma, [1, 1])
    Q = geo.patch([0], 0, 0)
    outside = next(j for j in range(c.N) if j not in Q)
    with pytest.raises(InputError, match="outside the patch"):
        compute_table(Decoder(c.HZ, Q), sigma, [outside])


@pytest.mark.parametrize("backend", [b for b in BACKENDS if b != "ilp_highs" or highs_available()])
def test_t11_check_with_all_four_qubits_fixed_to_zero(backend):
    g, c = _instance("petersen", 3)
    geo = c.geometry
    sigma = np.zeros(c.m, dtype=np.uint8)
    i = geo.z_check(4, 1)
    sigma[i] = 1
    phi = {j: 0 for j in geo.z_support(i)}
    if backend == "brute":
        Q = geo.patch(list(geo.edges[4]), 0, 2)  # small enough for brute
        r = Decoder(c.HZ, Q, backend).solve(sigma, phi)
    else:
        r = Decoder(c.HZ, None, backend).solve(sigma, phi)
    assert (r.status, r.cost, r.witness) == ("infeasible", None, None)


@pytest.mark.parametrize("backend", [b for b in BACKENDS if b != "ilp_highs" or highs_available()])
def test_t11_local_closed_with_inconsistent_parities(backend):
    g, c = _instance("cube", 3)
    geo = c.geometry
    v = 2
    Q = geo.patch([v], 1, 1)  # only the edge-qubits (E, e, 1) at v
    sigma = np.zeros(c.m, dtype=np.uint8)
    e = geo.edges_at[v][0]
    sigma[geo.z_check(e, 0)] = 1  # (e, 0) and (e, 1) both see only (E, e, 1) inside Q
    decoder = Decoder(c.HZ, Q, backend)
    r = decoder.solve(sigma)
    assert (r.status, r.cost, r.witness) == ("infeasible", None, None)
    table = compute_table(decoder, sigma, [geo.edge_qubit(e, 1)])
    assert all(row.status == "infeasible" and row.cost is None for row in table.rows.values())


def test_t11_global_scope_without_terminals_is_always_feasible():
    for name in NAMED_PRESETS:
        g = from_preset(name)
        for L in (2, 3, 4):
            c = codemod.build_code(g.n, g.edges, L)
            assert f2.rank(c.HZ) == c.m  # full row rank, so every syndrome is in the column space
    g, c = _instance("heawood", 2)
    rng = np.random.default_rng(11)
    r = Decoder(c.HZ).solve(rng.integers(0, 2, c.m))
    assert r.status == "optimal"
