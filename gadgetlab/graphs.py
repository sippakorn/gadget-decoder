"""Connected cubic graphs: presets, canonical edges, layouts, ball(), girth and alpha."""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field

import networkx as nx
import numpy as np

from gadgetlab.errors import InputError

NAMED_PRESETS = ("k4", "k33", "prism", "cube", "petersen", "heawood")
PRESETS = NAMED_PRESETS + ("random_cubic",)
CUSTOM = "custom"
RANDOM_RETRIES = 100

DISPLAY_NAMES = {
    "k4": "K4",
    "k33": "K3,3",
    "prism": "Prism",
    "cube": "Cube",
    "petersen": "Petersen",
    "heawood": "Heawood",
    "random_cubic": "Random cubic",
    CUSTOM: "Custom graph",
}


class GraphError(InputError):
    pass


@dataclass
class Graph:
    n: int
    edges: list[tuple[int, int]]  # canonical: u < w, sorted lexicographically
    preset: str  # a name from PRESETS, or "custom" for a user edge list
    seed: int | None = None  # the seed actually used, for random_cubic
    layout: dict[int, tuple[float, float]] = field(default_factory=dict)

    @property
    def nE(self) -> int:
        return len(self.edges)

    def display_name(self) -> str:
        name = DISPLAY_NAMES.get(self.preset, self.preset)
        if self.preset == "random_cubic":
            return f"{name} (n = {self.n}, seed {self.seed})"
        return name


# ---------------------------------------------------------------- edges


def canonical_edges(edge_list) -> list[tuple[int, int]]:
    """Each edge as (u, w) with u < w, the list sorted lexicographically (§3.2)."""
    out = []
    seen = set()
    for k, pair in enumerate(edge_list):
        try:
            u, w = pair
        except (TypeError, ValueError):
            raise GraphError(f"edge {k} must be a pair of vertices, got {pair!r}") from None
        if not (_is_int(u) and _is_int(w)):
            raise GraphError(f"edge {k} must be a pair of integers, got {pair!r}")
        u, w = int(u), int(w)
        if u == w:
            raise GraphError(f"edge {k} = ({u}, {w}) is a loop")
        e = (min(u, w), max(u, w))
        if e in seen:
            raise GraphError(f"edge ({e[0]}, {e[1]}) appears more than once")
        seen.add(e)
        out.append(e)
    return sorted(out)


def _is_int(x) -> bool:
    return isinstance(x, (int, np.integer)) and not isinstance(x, bool)


def validate_cubic(n: int, edges: list[tuple[int, int]]) -> None:
    """Raise GraphError unless (n, edges) is a connected cubic graph on 0..n-1 with canonical edges."""
    if not isinstance(n, int) or isinstance(n, bool):
        raise GraphError(f"n must be an integer, got {n!r}")
    if n < 4 or n % 2:
        raise GraphError(f"a cubic graph needs an even number of vertices, at least 4; got n = {n}")
    if list(edges) != canonical_edges(edges):
        raise GraphError("edges must be in canonical order: each edge (u, w) with u < w, sorted")
    for u, w in edges:
        if u < 0 or w >= n:
            raise GraphError(f"edge ({u}, {w}) uses a vertex outside 0..{n - 1}")
    degree = [0] * n
    for u, w in edges:
        degree[u] += 1
        degree[w] += 1
    for v, d in enumerate(degree):
        if d != 3:
            raise GraphError(f"the graph is not cubic: vertex {v} has degree {d}")
    components = nx.number_connected_components(to_networkx(n, edges))
    if components != 1:
        raise GraphError(f"the graph is not connected: it has {components} components")


def to_networkx(n: int, edges) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    return G


def incident_edges(n: int, edges) -> list[list[int]]:
    """For each vertex, the indices of its edges, ascending."""
    at: list[list[int]] = [[] for _ in range(n)]
    for e, (u, w) in enumerate(edges):
        at[u].append(e)
        at[w].append(e)
    return at


def neighbours(n: int, edges) -> list[list[int]]:
    nb: list[list[int]] = [[] for _ in range(n)]
    for u, w in edges:
        nb[u].append(w)
        nb[w].append(u)
    return nb


# ---------------------------------------------------------------- construction


def from_preset(name: str, n: int | None = None, seed: int | None = None) -> Graph:
    if name not in PRESETS:
        raise GraphError(f"unknown preset {name!r}; expected one of {', '.join(PRESETS)}")
    if name == "random_cubic":
        return random_cubic(n, 0 if seed is None else seed)
    if name == "k4":
        G = nx.complete_graph(4)
    elif name == "k33":
        G = nx.complete_bipartite_graph(3, 3)
    elif name == "prism":
        G = nx.circular_ladder_graph(3)
    elif name == "cube":
        G = nx.convert_node_labels_to_integers(nx.hypercube_graph(3), ordering="sorted")
    elif name == "petersen":
        G = nx.petersen_graph()
    else:
        G = nx.heawood_graph()
    n_g = G.number_of_nodes()
    edges = canonical_edges(G.edges())
    validate_cubic(n_g, edges)
    return Graph(n=n_g, edges=edges, preset=name, seed=None, layout=default_layout(name, n_g, edges))


def random_cubic(n, seed: int = 0) -> Graph:
    if n is None:
        raise GraphError("random_cubic needs n")
    if not isinstance(n, int) or isinstance(n, bool) or n < 4 or n % 2:
        raise GraphError(f"random_cubic needs an even n of at least 4; got {n!r}")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise GraphError(f"random_cubic needs a non-negative integer seed; got {seed!r}")
    s = seed
    for _ in range(RANDOM_RETRIES + 1):
        G = nx.random_regular_graph(3, n, seed=s)
        if nx.is_connected(G):
            edges = canonical_edges(G.edges())
            validate_cubic(n, edges)
            return Graph(n=n, edges=edges, preset="random_cubic", seed=s, layout=spring_layout(n, edges))
        s += 1
    raise GraphError(f"no connected cubic graph found for n = {n} with seeds {seed}..{s - 1}")


def from_edges(edge_list, n: int | None = None) -> Graph:
    """A user edge list; vertices must be 0..n-1."""
    edges = canonical_edges(edge_list)
    if not edges:
        raise GraphError("the edge list is empty")
    if n is None:
        n = max(w for _, w in edges) + 1
    validate_cubic(n, edges)
    return Graph(n=n, edges=edges, preset=CUSTOM, seed=None, layout=spring_layout(n, edges))


# ---------------------------------------------------------------- layouts


def _circle(count: int, radius: float, cx: float = 0.5, cy: float = 0.5) -> list[tuple[float, float]]:
    """Points on a circle, the first at the top, then clockwise on screen (y grows downwards)."""
    pts = []
    for i in range(count):
        a = math.radians(-90.0 + 360.0 * i / count)
        pts.append((cx + radius * math.cos(a), cy + radius * math.sin(a)))
    return pts


def _finish(points) -> dict[int, tuple[float, float]]:
    return {v: (round(min(1.0, max(0.0, x)), 4), round(min(1.0, max(0.0, y)), 4)) for v, (x, y) in enumerate(points)}


def default_layout(preset: str, n: int, edges) -> dict[int, tuple[float, float]]:
    """Deterministic layout in [0, 1]^2: standard drawings for named presets, else a spring layout."""
    if preset == "k4":
        outer = _circle(3, 0.5, cy=0.58)
        return _finish(outer + [(0.5, 0.58)])
    if preset == "k33":
        return _finish([(0.1, 0.1), (0.5, 0.1), (0.9, 0.1), (0.1, 0.9), (0.5, 0.9), (0.9, 0.9)])
    if preset == "prism":
        return _finish(_circle(3, 0.5, cy=0.58) + _circle(3, 0.2, cy=0.58))
    if preset == "cube":
        pts = {0: (0.0, 0.0), 1: (1.0, 0.0), 3: (1.0, 1.0), 2: (0.0, 1.0),
               4: (0.3, 0.3), 5: (0.7, 0.3), 7: (0.7, 0.7), 6: (0.3, 0.7)}
        return _finish([pts[v] for v in range(8)])
    if preset == "petersen":
        return _finish(_circle(5, 0.5) + _circle(5, 0.24))
    if preset == "heawood":
        return _finish(_circle(14, 0.5))
    return spring_layout(n, edges)


def spring_layout(n: int, edges) -> dict[int, tuple[float, float]]:
    """networkx spring layout with seed 0, scaled uniformly into [0, 1]^2 and centred."""
    pos = nx.spring_layout(to_networkx(n, edges), seed=0)
    xs = [float(pos[v][0]) for v in range(n)]
    ys = [float(pos[v][1]) for v in range(n)]
    span = max(max(xs) - min(xs), max(ys) - min(ys)) or 1.0
    ox = (1.0 - (max(xs) - min(xs)) / span) / 2
    oy = (1.0 - (max(ys) - min(ys)) / span) / 2
    return _finish([(ox + (xs[v] - min(xs)) / span, oy + (ys[v] - min(ys)) / span) for v in range(n)])


# ---------------------------------------------------------------- graph quantities


def ball(n: int, edges, s: int, r: int) -> list[int]:
    """Vertices at graph distance at most r from s, ascending."""
    if not 0 <= s < n:
        raise GraphError(f"ball: vertex {s} is out of range 0..{n - 1}")
    if r < 0:
        raise GraphError(f"ball: radius must be non-negative, got {r}")
    nb = neighbours(n, edges)
    dist = {s: 0}
    queue = deque([s])
    while queue:
        u = queue.popleft()
        if dist[u] == r:
            continue
        for w in nb[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                queue.append(w)
    return sorted(dist)


def girth(n: int, edges) -> int | None:
    """Length of a shortest cycle (None for a forest), by a BFS from every vertex."""
    nb = neighbours(n, edges)
    best = math.inf
    for s in range(n):
        dist = [-1] * n
        parent = [-1] * n
        dist[s] = 0
        queue = deque([s])
        while queue:
            u = queue.popleft()
            if 2 * dist[u] + 1 >= best:
                break
            for w in nb[u]:
                if dist[w] < 0:
                    dist[w] = dist[u] + 1
                    parent[w] = u
                    queue.append(w)
                elif w != parent[u]:
                    best = min(best, dist[u] + dist[w] + 1)
    return None if best is math.inf else int(best)


def independence_number(n: int, edges) -> int:
    """alpha(G), the size of a largest independent set, by branching on bitmasks."""
    nb_mask = [0] * n
    for u, w in edges:
        nb_mask[u] |= 1 << w
        nb_mask[w] |= 1 << u
    memo: dict[int, int] = {}

    def mis(S: int) -> int:
        if S == 0:
            return 0
        if S in memo:
            return memo[S]
        # pick the vertex of S with the fewest neighbours inside S
        best_v, best_d = -1, n + 1
        T = S
        while T:
            v = (T & -T).bit_length() - 1
            T &= T - 1
            d = bin(nb_mask[v] & S).count("1")
            if d < best_d:
                best_v, best_d = v, d
        v = best_v
        if best_d <= 1:
            value = 1 + mis(S & ~(1 << v) & ~nb_mask[v])
        else:
            value = max(1 + mis(S & ~(1 << v) & ~nb_mask[v]), mis(S & ~(1 << v)))
        memo[S] = value
        return value

    return mis((1 << n) - 1)
