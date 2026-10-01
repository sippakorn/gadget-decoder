"""T3 solver cross-check: ilp_cbc and brute agree exactly on status and cost, and verify() passes
on every witness. Also ilp_highs when highspy is installed."""

import numpy as np
import pytest

from gadgetlab import code as codemod
from gadgetlab.graphs import ball, from_preset
from gadgetlab.solve import Decoder, Problem, highs_available, verify

INSTANCES = [("k4", 2), ("k4", 3), ("k33", 2)]


def _syndromes(c, rng):
    """150 syndromes H_Z x for random x of weight 1..6, then 50 uniformly random ones."""
    out = []
    for _ in range(150):
        x = np.zeros(c.N, dtype=np.int64)
        x[rng.choice(c.N, int(rng.integers(1, 7)), replace=False)] = 1
        out.append(((c.HZ.astype(np.int64) @ x) % 2).astype(np.uint8))
    for _ in range(50):
        out.append(rng.integers(0, 2, c.m).astype(np.uint8))
    return out


def _random_patch(g, geo, rng):
    if rng.random() < 0.5:
        U = ball(g.n, g.edges, int(rng.integers(g.n)), int(rng.integers(0, 2)))
    else:
        U = sorted(rng.choice(g.n, int(rng.integers(1, g.n)), replace=False).tolist())
    l0 = int(rng.integers(0, geo.L))
    l1 = int(rng.integers(l0, geo.L))
    return tuple(geo.patch(U, l0, l1))


def _terminals(pool, rng):
    k = int(rng.integers(0, 3))
    chosen = rng.choice(pool, k, replace=False) if k else []
    return {int(j): int(rng.integers(0, 2)) for j in chosen}


@pytest.mark.parametrize("preset,L", INSTANCES)
def test_t3_ilp_cbc_and_brute_agree(preset, L, record):
    g = from_preset(preset)
    c = codemod.build_code(g.n, g.edges, L)
    geo = c.geometry
    rng = np.random.default_rng(3 + 10 * L + len(preset))
    backends = ["ilp_cbc", "brute"] + (["ilp_highs"] if highs_available() else [])
    global_decoders = {b: Decoder(c.HZ, None, b) for b in backends}
    stats = {"global": [0, 0], "local_closed": [0, 0]}
    for sigma in _syndromes(c, rng):
        phi = _terminals(np.arange(c.N), rng)
        Q = _random_patch(g, geo, rng)
        local_decoders = {b: Decoder(c.HZ, Q, b) for b in backends}
        local_phi = _terminals(np.array(Q), rng) if Q else {}
        for scope, decoders, patch, terms in (("global", global_decoders, None, phi),
                                              ("local_closed", local_decoders, Q, local_phi)):
            results = {b: d.solve(sigma, terms) for b, d in decoders.items()}
            problem = Problem(c.HZ, sigma, patch, tuple(sorted(terms.items())))
            for r in results.values():
                verify(problem, r)
            outcomes = {b: (r.status, r.cost) for b, r in results.items()}
            assert len(set(outcomes.values())) == 1, (scope, outcomes)
            stats[scope][0] += 1
            stats[scope][1] += results["brute"].status == "infeasible"
    record("T3 solver cross-check", f"{preset} L={L}: backends {', '.join(backends)}; "
           f"global {stats['global'][0]} problems ({stats['global'][1]} infeasible), "
           f"local_closed {stats['local_closed'][0]} problems ({stats['local_closed'][1]} infeasible); all agree")
    assert stats["global"][0] == 200 and stats["local_closed"][0] == 200
    assert stats["local_closed"][1] > 0  # the infeasible path was exercised


def test_brute_refuses_large_solution_spaces():
    g = from_preset("petersen")
    c = codemod.build_code(g.n, g.edges, 3)
    from gadgetlab.solve import BackendRefused
    with pytest.raises(BackendRefused, match="dimension 35"):
        Decoder(c.HZ, None, "brute").solve(np.zeros(c.m, dtype=np.uint8))


def test_verify_rejects_bad_witnesses():
    g = from_preset("k4")
    c = codemod.build_code(g.n, g.edges, 2)
    sigma = np.zeros(c.m, dtype=np.uint8)
    sigma[0] = 1
    r = Decoder(c.HZ, None, "ilp_cbc").solve(sigma)
    problem = Problem(c.HZ, sigma, None, ())
    from dataclasses import replace

    from gadgetlab.solve import VerificationError
    with pytest.raises(VerificationError, match="weight"):
        verify(problem, replace(r, cost=r.cost + 1))
    flipped = r.witness.copy()
    flipped[np.flatnonzero(c.HZ[0])[0]] ^= 1
    with pytest.raises(VerificationError):
        verify(problem, replace(r, witness=flipped))
    Q = tuple(range(4))
    with pytest.raises(VerificationError, match="outside the patch"):
        verify(Problem(c.HZ, sigma, Q, ()), replace(r, witness=np.eye(c.N, dtype=np.uint8)[10], cost=0))
