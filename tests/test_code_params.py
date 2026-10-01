"""T2 code parameters: k = n/2 + 1, girth as networkx, dX and dZ as warmup's min_logical_weight.

A failure here contradicts a recorded research claim: report it, never change the expectation.
"""

import networkx as nx
import pytest

from gadgetlab import code as codemod
from gadgetlab import distance
from gadgetlab.graphs import NAMED_PRESETS, from_preset, to_networkx


def _graphs():
    out = [(name, None) for name in NAMED_PRESETS]
    out += [("random_cubic", (n, seed)) for n in (6, 8, 10, 12) for seed in (0, 1)]
    return out


def _cases():
    cases = []
    for name, rnd in _graphs():
        g = from_preset(name) if rnd is None else from_preset(name, n=rnd[0], seed=rnd[1])
        for L in (2, 3, 4):
            N = 3 * g.n // 2 * L + g.n * (L - 1)
            label = f"{name}{'' if rnd is None else f'-n{rnd[0]}-s{rnd[1]}'}-L{L}"
            marks = [pytest.mark.slow] if N > 40 else []
            cases.append(pytest.param(name, rnd, L, id=label, marks=marks))
    return cases


@pytest.mark.parametrize("name,rnd,L", _cases())
def test_t2_code_parameters(name, rnd, L, warmup, record):
    g = from_preset(name) if rnd is None else from_preset(name, n=rnd[0], seed=rnd[1])
    c = codemod.build_code(g.n, g.edges, L)
    info = codemod.code_info(g, L, c)
    girth_nx = nx.girth(to_networkx(g.n, g.edges))
    assert info["k"] == g.n // 2 + 1
    assert info["girth"] == girth_nx

    d = distance.distances(c.HX, c.HZ)
    label = f"{name}{'' if rnd is None else f' n={rnd[0]} seed={g.seed}'} L={L}"
    record("T2 distances (dX, dZ)", f"{label:28s} N={c.N:3d} k={info['k']} girth={girth_nx} dX={d['dX']} dZ={d['dZ']}")
    if c.N > 60:
        assert d == {"dX": "not computed", "dZ": "not computed"}
        return
    ref_dX = warmup["min_logical_weight"](c.HZ.astype(int), c.HX.astype(int), 4)
    ref_dZ = warmup["min_logical_weight"](c.HX.astype(int), c.HZ.astype(int), 4)
    assert d["dX"] == ref_dX
    assert d["dZ"] == ref_dZ
    target = min(L, girth_nx)
    if target <= 4:
        as_int = [5 if v == ">4" else v for v in (d["dX"], d["dZ"])]
        assert min(as_int) == target
