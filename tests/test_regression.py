"""T5 regression of a recorded claim: sigma = 1 on every Z-check (e, 0), 0 elsewhere, gives
MW(sigma) = 3n/2 - 2 alpha(G). Global scope, ilp_cbc. A failure is reported, never "fixed"."""

import networkx as nx
import numpy as np
import pytest

from gadgetlab import code as codemod
from gadgetlab.graphs import NAMED_PRESETS, from_preset, to_networkx
from gadgetlab.solve import Decoder, Problem, verify


def _cases():
    out = []
    for name in NAMED_PRESETS:
        out += [pytest.param(name, None, L, id=f"{name}-L{L}") for L in (2, 3, 4)]
    for n in (6, 8, 10, 12):
        for seed in (0, 1):
            out += [pytest.param("random_cubic", (n, seed), L, id=f"random-n{n}-s{seed}-L{L}") for L in (2, 3, 4)]
    return out


@pytest.mark.parametrize("name,rnd,L", _cases())
def test_t5_gap0_all_ones(name, rnd, L, record):
    g = from_preset(name) if rnd is None else from_preset(name, n=rnd[0], seed=rnd[1])
    c = codemod.build_code(g.n, g.edges, L)
    geo = c.geometry
    sigma = np.zeros(geo.m, dtype=np.uint8)
    sigma[[geo.z_check(e, 0) for e in range(geo.nE)]] = 1
    result = Decoder(c.HZ, None, "ilp_cbc").solve(sigma)
    verify(Problem(c.HZ, sigma, None, ()), result)
    G = to_networkx(g.n, g.edges)
    alpha = max(len(q) for q in nx.find_cliques(nx.complement(G)))
    predicted = 3 * g.n // 2 - 2 * alpha
    label = f"{name}{'' if rnd is None else f' n={rnd[0]} seed={g.seed}'} L={L}"
    record("T5 MW(gap 0 all ones) vs 3n/2 - 2 alpha", f"{label:28s} MW={result.cost:3d} 3n/2-2a={predicted:3d} "
           f"{'match' if result.cost == predicted else 'MISMATCH'}")
    assert result.status == "optimal"
    assert result.cost == predicted
