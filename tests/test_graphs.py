"""E1: presets, random_cubic, user edge lists, rejections, layouts, ball(), girth and alpha."""

import networkx as nx
import pytest

from gadgetlab import graphs
from gadgetlab.graphs import GraphError


@pytest.mark.parametrize("name,n", [("k4", 4), ("k33", 6), ("prism", 6), ("cube", 8), ("petersen", 10), ("heawood", 14)])
def test_presets_are_canonical_connected_cubic(name, n):
    g = graphs.from_preset(name)
    assert g.n == n and g.preset == name and g.seed is None
    assert g.edges == sorted(g.edges) and all(u < w for u, w in g.edges)
    graphs.validate_cubic(g.n, g.edges)
    assert set(g.layout) == set(range(n))
    assert all(0.0 <= x <= 1.0 and 0.0 <= y <= 1.0 for x, y in g.layout.values())
    assert graphs.default_layout(name, g.n, g.edges) == g.layout  # deterministic


def test_petersen_is_outer_pentagon_and_inner_pentagram():
    g = graphs.from_preset("petersen")
    assert g.edges[:3] == [(0, 1), (0, 4), (0, 5)]
    assert g.layout[0] == (0.5, 0.0)
    assert g.layout[1] == pytest.approx((0.9755, 0.3455), abs=1e-4)
    inner = [g.layout[v] for v in range(5, 10)]
    assert all(abs(x - 0.5) < 0.3 and abs(y - 0.5) < 0.3 for x, y in inner)


def test_presets_match_networkx_constructions():
    assert graphs.from_preset("k4").edges == graphs.canonical_edges(nx.complete_graph(4).edges())
    cube = nx.convert_node_labels_to_integers(nx.hypercube_graph(3), ordering="sorted")
    assert graphs.from_preset("cube").edges == graphs.canonical_edges(cube.edges())
    assert graphs.from_preset("heawood").edges == graphs.canonical_edges(nx.heawood_graph().edges())


def test_random_cubic_is_deterministic_and_records_the_seed():
    a = graphs.from_preset("random_cubic", n=12, seed=5)
    b = graphs.from_preset("random_cubic", n=12, seed=5)
    assert a.edges == b.edges and a.layout == b.layout and a.seed == 5
    assert a.display_name() == "Random cubic (n = 12, seed 5)"


def test_random_cubic_retries_disconnected_seeds(monkeypatch):
    calls = []
    real = nx.random_regular_graph

    def fake(d, n, seed):
        calls.append(seed)
        if len(calls) < 3:
            return nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
        return real(d, n, seed=seed)

    monkeypatch.setattr(graphs.nx, "random_regular_graph", fake)
    g = graphs.from_preset("random_cubic", n=8, seed=7)
    assert calls == [7, 8, 9] and g.seed == 9


@pytest.mark.parametrize("edges,message", [
    ([(0, 1), (1, 2), (2, 3), (3, 0)], "not cubic"),
    ([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7)], "not connected"),
    ([(0, 0)], "loop"),
    ([(0, 1), (1, 0)], "more than once"),
])
def test_user_edge_lists_are_rejected_with_clear_messages(edges, message):
    with pytest.raises(GraphError, match=message):
        graphs.from_edges(edges)


def test_user_edge_list_is_canonicalised():
    g = graphs.from_edges([(3, 2), (1, 0), (2, 0), (3, 1), (0, 3), (2, 1)])
    assert g.edges == [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] and g.preset == "custom"
    assert all(0.0 <= c <= 1.0 for xy in g.layout.values() for c in xy)


def test_ball_and_girth_and_alpha():
    g = graphs.from_preset("petersen")
    assert graphs.ball(g.n, g.edges, 0, 0) == [0]
    assert graphs.ball(g.n, g.edges, 0, 1) == [0, 1, 4, 5]
    assert graphs.ball(g.n, g.edges, 0, 2) == list(range(10))
    for name in graphs.NAMED_PRESETS:
        h = graphs.from_preset(name)
        G = graphs.to_networkx(h.n, h.edges)
        assert graphs.girth(h.n, h.edges) == nx.girth(G)
        reference = max(len(c) for c in nx.find_cliques(nx.complement(G)))
        assert graphs.independence_number(h.n, h.edges) == reference
