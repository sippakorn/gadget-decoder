"""Typed addresses <-> indices, supports, patches and panel data for HGP(cycle code of G, R_L).

This is the adapter between the picture of G □ P_L (§3.6) and the index world of the
engine (§3.5). Only the session, server, CLI and code construction use it.
"""

from __future__ import annotations

import numpy as np

from gadgetlab.errors import InputError
from gadgetlab.graphs import incident_edges

L_MIN, L_MAX = 2, 8


class AddressError(InputError):
    pass


def _int_field(obj: dict, key: str, field: str, lo: int, hi: int) -> int:
    if key not in obj:
        raise AddressError(f"{field}: missing key {key!r}")
    value = obj[key]
    if not isinstance(value, int) or isinstance(value, bool):
        raise AddressError(f"{field}.{key}: expected an integer, got {value!r}")
    if not lo <= value <= hi:
        raise AddressError(f"{field}.{key}: {value} is out of range {lo}..{hi}")
    return value


class Geometry:
    """Index map and supports of §3.5 for a graph (n, canonical edges) and L layers."""

    def __init__(self, n: int, edges, L: int):
        if not isinstance(L, int) or isinstance(L, bool) or not L_MIN <= L <= L_MAX:
            raise InputError(f"L must be an integer between {L_MIN} and {L_MAX}, got {L!r}")
        self.n = n
        self.edges = [tuple(e) for e in edges]
        self.L = L
        self.nE = len(self.edges)
        self.N = self.nE * L + n * (L - 1)
        self.m = self.nE * (L - 1)
        self.mx = n * L
        self.edges_at = incident_edges(n, self.edges)

    # ------------------------------------------------------------ index formulas

    def edge_qubit(self, e: int, l: int) -> int:
        return e * self.L + l

    def vertex_qubit(self, v: int, t: int) -> int:
        return self.nE * self.L + v * (self.L - 1) + t

    def z_check(self, e: int, t: int) -> int:
        return e * (self.L - 1) + t

    def x_check(self, v: int, l: int) -> int:
        return v * self.L + l

    def qubit(self, j: int) -> tuple:
        """("E", e, l) or ("V", v, t) for qubit index j."""
        if not 0 <= j < self.N:
            raise AddressError(f"qubit index {j} is out of range 0..{self.N - 1}")
        if j < self.nE * self.L:
            e, l = divmod(j, self.L)
            return ("E", e, l)
        v, t = divmod(j - self.nE * self.L, self.L - 1)
        return ("V", v, t)

    def check(self, i: int) -> tuple[int, int]:
        """(e, t) for Z-check index i."""
        if not 0 <= i < self.m:
            raise AddressError(f"check index {i} is out of range 0..{self.m - 1}")
        return divmod(i, self.L - 1)

    # ------------------------------------------------------------ supports

    def z_support(self, i: int) -> list[int]:
        e, t = self.check(i)
        u, w = self.edges[e]
        return [self.edge_qubit(e, t), self.edge_qubit(e, t + 1), self.vertex_qubit(u, t), self.vertex_qubit(w, t)]

    def x_support(self, k: int) -> list[int]:
        v, l = divmod(k, self.L)
        qubits = [self.edge_qubit(e, l) for e in self.edges_at[v]]
        qubits += [self.vertex_qubit(v, t) for t in (l - 1, l) if 0 <= t <= self.L - 2]
        return qubits

    def HZ(self) -> np.ndarray:
        H = np.zeros((self.m, self.N), dtype=np.uint8)
        for i in range(self.m):
            H[i, self.z_support(i)] = 1
        return H

    def HX(self) -> np.ndarray:
        H = np.zeros((self.mx, self.N), dtype=np.uint8)
        for k in range(self.mx):
            H[k, self.x_support(k)] = 1
        return H

    # ------------------------------------------------------------ patches (§3.9)

    def patch(self, vertices, l0: int, l1: int) -> list[int]:
        """Qubits of Q(U, l0, l1), ascending."""
        U = set(vertices)
        if not 0 <= l0 <= l1 <= self.L - 1:
            raise InputError(f"patch layers must satisfy 0 <= l0 <= l1 <= {self.L - 1}, got [{l0}, {l1}]")
        out = []
        for e, (u, w) in enumerate(self.edges):
            if u in U or w in U:
                out += [self.edge_qubit(e, l) for l in range(l0, l1 + 1)]
        for v in sorted(U):
            out += [self.vertex_qubit(v, t) for t in range(l0, l1)]
        return sorted(out)

    # ------------------------------------------------------------ JSON addresses (§5)

    def qubit_json(self, j: int) -> dict:
        kind, a, b = self.qubit(j)
        return {"q": "E", "e": a, "l": b} if kind == "E" else {"q": "V", "v": a, "t": b}

    def check_json(self, i: int) -> dict:
        e, t = self.check(i)
        return {"e": e, "t": t}

    def parse_qubit(self, obj, field: str) -> int:
        if not isinstance(obj, dict):
            raise AddressError(f'{field}: expected a qubit address like {{"q": "E", "e": 3, "l": 0}}, got {obj!r}')
        kind = obj.get("q")
        if kind == "E":
            _only_keys(obj, ("q", "e", "l"), field)
            e = _int_field(obj, "e", field, 0, self.nE - 1)
            l = _int_field(obj, "l", field, 0, self.L - 1)
            return self.edge_qubit(e, l)
        if kind == "V":
            _only_keys(obj, ("q", "v", "t"), field)
            v = _int_field(obj, "v", field, 0, self.n - 1)
            t = _int_field(obj, "t", field, 0, self.L - 2)
            return self.vertex_qubit(v, t)
        raise AddressError(f'{field}.q: expected "E" or "V", got {kind!r}')

    def parse_check(self, obj, field: str) -> int:
        if not isinstance(obj, dict):
            raise AddressError(f'{field}: expected a Z-check address like {{"e": 3, "t": 0}}, got {obj!r}')
        _only_keys(obj, ("e", "t"), field)
        e = _int_field(obj, "e", field, 0, self.nE - 1)
        t = _int_field(obj, "t", field, 0, self.L - 2)
        return self.z_check(e, t)

    # ------------------------------------------------------------ text

    def describe_qubit(self, j: int) -> str:
        kind, a, b = self.qubit(j)
        if kind == "E":
            u, w = self.edges[a]
            return f"edge-qubit e = {a} ({u}–{w}), l = {b}"
        return f"vertex-qubit v = {a}, t = {b}"

    def describe_check(self, i: int) -> str:
        e, t = self.check(i)
        u, w = self.edges[e]
        return f"Z-check e = {e} ({u}–{w}), t = {t}"

    # ------------------------------------------------------------ panels (§3.6, V2)

    def panels(self) -> list[dict]:
        """What the viewer draws: Layer 0, Gap 0, Layer 1, ..., Layer L-1, each with the same layout."""
        out = []
        for l in range(self.L):
            out.append({
                "kind": "layer", "index": l,
                "title": f"Layer {l}: edges are qubits",
                "edges": [{"edge": e, "qubit": self.edge_qubit(e, l)} for e in range(self.nE)],
                "vertices": [{"vertex": v} for v in range(self.n)],
            })
            if l < self.L - 1:
                t = l
                out.append({
                    "kind": "gap", "index": t,
                    "title": f"Gap {t}: edges are checks, vertices are qubits",
                    "edges": [{"edge": e, "check": self.z_check(e, t)} for e in range(self.nE)],
                    "vertices": [{"vertex": v, "qubit": self.vertex_qubit(v, t)} for v in range(self.n)],
                })
        return out


def _only_keys(obj: dict, keys: tuple, field: str) -> None:
    extra = sorted(set(obj) - set(keys))
    if extra:
        raise AddressError(f"{field}: unexpected key {extra[0]!r}")
