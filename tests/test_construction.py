"""T1 construction and T4 single-qubit syndromes; index map and panel data (E2)."""

import numpy as np
import pytest

from gadgetlab import code as codemod
from gadgetlab import geometry
from gadgetlab.graphs import NAMED_PRESETS, from_preset

CASES = [(p, L) for p in NAMED_PRESETS for L in (2, 3, 4)]


@pytest.mark.parametrize("preset,L", CASES)
def test_t1_kronecker_equals_supports(preset, L):
    g = from_preset(preset)
    c = codemod.build_code(g.n, g.edges, L)
    geo = c.geometry
    HX_kron, HZ_kron = codemod.hgp(codemod.incidence_matrix(g.n, g.edges), codemod.rep(L))
    assert np.array_equal(HZ_kron, geo.HZ())
    assert np.array_equal(HX_kron, geo.HX())
    assert not ((HX_kron.astype(int) @ HZ_kron.T.astype(int)) % 2).any()
    nE = 3 * g.n // 2
    assert len(g.edges) == nE
    assert HZ_kron.shape == (nE * (L - 1), nE * L + g.n * (L - 1))
    assert HX_kron.shape == (g.n * L, nE * L + g.n * (L - 1))
    assert (geo.N, geo.m, geo.mx) == (nE * L + g.n * (L - 1), nE * (L - 1), g.n * L)


def test_t1_construction_guard_catches_a_wrong_index_formula(monkeypatch):
    g = from_preset("k4")
    codemod._build_code.cache_clear()
    monkeypatch.setattr(geometry.Geometry, "vertex_qubit",
                        lambda self, v, t: self.nE * self.L + t * self.n + v)  # transposed on purpose
    with pytest.raises(codemod.ConstructionError):
        codemod.build_code(g.n, g.edges, 3)
    codemod._build_code.cache_clear()


def test_index_map_round_trip_and_formulas():
    g = from_preset("petersen")
    geo = geometry.Geometry(g.n, g.edges, 3)
    for j in range(geo.N):
        kind, a, b = geo.qubit(j)
        assert (geo.edge_qubit(a, b) if kind == "E" else geo.vertex_qubit(a, b)) == j
        assert geo.parse_qubit(geo.qubit_json(j), "q") == j
    for i in range(geo.m):
        e, t = geo.check(i)
        assert i == e * (geo.L - 1) + t
        assert geo.parse_check(geo.check_json(i), "c") == i
    assert geo.edge_qubit(3, 1) == 3 * 3 + 1
    assert geo.vertex_qubit(2, 1) == 15 * 3 + 2 * 2 + 1
    assert geo.x_check(4, 2) == 4 * 3 + 2


def test_supports_follow_section_3_5():
    g = from_preset("k33")
    L = 4
    geo = geometry.Geometry(g.n, g.edges, L)
    for e, (u, w) in enumerate(g.edges):
        for t in range(L - 1):
            expected = {("E", e, t), ("E", e, t + 1), ("V", u, t), ("V", w, t)}
            assert {geo.qubit(j) for j in geo.z_support(geo.z_check(e, t))} == expected
    for v in range(g.n):
        for l in range(L):
            expected = {("E", e, l) for e, (a, b) in enumerate(g.edges) if v in (a, b)}
            expected |= {("V", v, t) for t in (l - 1, l) if 0 <= t <= L - 2}
            assert {geo.qubit(j) for j in geo.x_support(geo.x_check(v, l))} == expected


def test_t4_single_qubit_syndromes_k4_L3():
    g = from_preset("k4")
    L = 3
    c = codemod.build_code(g.n, g.edges, L)
    geo = c.geometry
    for j in range(geo.N):
        x = np.zeros(geo.N, dtype=int)
        x[j] = 1
        fired = {geo.check(int(i)) for i in np.flatnonzero((c.HZ.astype(int) @ x) % 2)}
        kind, a, b = geo.qubit(j)
        if kind == "E":
            expected = {(a, t) for t in (b - 1, b) if 0 <= t <= L - 2}
        else:
            expected = {(e, b) for e, (u, w) in enumerate(g.edges) if a in (u, w)}
            assert len(expected) == 3
        assert fired == expected, f"qubit {geo.qubit(j)}"


def test_panels_cover_every_qubit_and_check_once():
    g = from_preset("cube")
    geo = geometry.Geometry(g.n, g.edges, 4)
    panels = geo.panels()
    assert [(p["kind"], p["index"]) for p in panels] == [
        ("layer", 0), ("gap", 0), ("layer", 1), ("gap", 1), ("layer", 2), ("gap", 2), ("layer", 3)]
    assert panels[1]["title"] == "Gap 0: edges are checks, vertices are qubits"
    qubits = [m["qubit"] for p in panels for m in p["edges"] + p["vertices"] if "qubit" in m]
    checks = [m["check"] for p in panels for m in p["edges"] if "check" in m]
    assert sorted(qubits) == list(range(geo.N))
    assert sorted(checks) == list(range(geo.m))


def test_patch_follows_section_3_9():
    g = from_preset("petersen")
    geo = geometry.Geometry(g.n, g.edges, 4)
    U = {0, 1}
    Q = set(geo.patch(U, 1, 2))
    expected = {geo.edge_qubit(e, l) for e, (u, w) in enumerate(g.edges) if u in U or w in U for l in (1, 2)}
    expected |= {geo.vertex_qubit(v, 1) for v in U}
    assert Q == expected
    assert geo.patch(U, 2, 2) == sorted(geo.edge_qubit(e, 2) for e, (u, w) in enumerate(g.edges) if u in U or w in U)
    with pytest.raises(ValueError):
        geo.patch(U, 2, 1)


@pytest.mark.parametrize("bad,field", [
    ({"q": "E", "e": 15, "l": 0}, "terminals[0].e"),
    ({"q": "E", "e": 1, "l": 3}, "terminals[0].l"),
    ({"q": "V", "v": 10, "t": 0}, "terminals[0].v"),
    ({"q": "V", "v": 1, "t": 2}, "terminals[0].t"),
    ({"q": "X", "e": 1, "l": 0}, "terminals[0].q"),
    ({"q": "E", "e": 1}, "terminals[0]: missing key 'l'"),
    ({"q": "E", "e": 1, "l": 0, "t": 0}, "terminals[0]: unexpected key"),
])
def test_qubit_addresses_are_validated_with_field_names(bad, field):
    g = from_preset("petersen")
    geo = geometry.Geometry(g.n, g.edges, 3)
    with pytest.raises(geometry.AddressError, match=field.replace("[", r"\[").replace("]", r"\]")):
        geo.parse_qubit(bad, "terminals[0]")
