"""H1, R_L, hgp() and code info (§3.3, §3.4, E2, E3)."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from gadgetlab import f2
from gadgetlab.geometry import Geometry
from gadgetlab.graphs import girth


class ConstructionError(AssertionError):
    """The Kronecker and explicit-support constructions disagree, or H_X H_Z^T != 0."""


def incidence_matrix(n: int, edges) -> np.ndarray:
    """H1: n × nE, H1[v, e] = 1 exactly when v is an endpoint of edge e."""
    H1 = np.zeros((n, len(edges)), dtype=np.uint8)
    for e, (u, w) in enumerate(edges):
        H1[u, e] = H1[w, e] = 1
    return H1


def rep(L: int) -> np.ndarray:
    """R_L: (L-1) × L with R_L[i, i] = R_L[i, i+1] = 1."""
    R = np.zeros((L - 1, L), dtype=np.uint8)
    for i in range(L - 1):
        R[i, i] = R[i, i + 1] = 1
    return R


def hgp(H1, H2) -> tuple[np.ndarray, np.ndarray]:
    """(H_X, H_Z) by the Kronecker formula, exactly as hgp() in warmup_check.py."""
    H1 = np.asarray(H1, dtype=np.int64)
    H2 = np.asarray(H2, dtype=np.int64)
    r1, n1 = H1.shape
    r2, n2 = H2.shape
    HX = np.hstack([np.kron(H1, np.eye(n2, dtype=np.int64)), np.kron(np.eye(r1, dtype=np.int64), H2.T)]) % 2
    HZ = np.hstack([np.kron(np.eye(n1, dtype=np.int64), H2), np.kron(H1.T, np.eye(r2, dtype=np.int64))]) % 2
    return HX.astype(np.uint8), HZ.astype(np.uint8)


@dataclass(frozen=True)
class Code:
    geometry: Geometry
    HX: np.ndarray
    HZ: np.ndarray

    @property
    def N(self) -> int:
        return self.geometry.N

    @property
    def m(self) -> int:
        return self.geometry.m


def build_code(n: int, edges, L: int) -> Code:
    """Build H_X and H_Z both ways (Kronecker and explicit supports) and require them to agree."""
    return _build_code(n, tuple(tuple(e) for e in edges), L)


@lru_cache(maxsize=32)
def _build_code(n: int, edges: tuple, L: int) -> Code:
    geo = Geometry(n, edges, L)
    HX, HZ = hgp(incidence_matrix(n, edges), rep(L))
    HX_explicit, HZ_explicit = geo.HX(), geo.HZ()
    if HZ.shape != (geo.m, geo.N) or HX.shape != (geo.mx, geo.N):
        raise ConstructionError(f"unexpected sizes: H_Z {HZ.shape}, H_X {HX.shape}, expected m = {geo.m}, n·L = {geo.mx}, N = {geo.N}")
    if not np.array_equal(HZ, HZ_explicit):
        raise ConstructionError("H_Z from the Kronecker formula differs from H_Z built from the supports of §3.5")
    if not np.array_equal(HX, HX_explicit):
        raise ConstructionError("H_X from the Kronecker formula differs from H_X built from the supports of §3.5")
    if ((HX.astype(np.int64) @ HZ.T.astype(np.int64)) % 2).any():
        raise ConstructionError("H_X H_Z^T is not 0 mod 2")
    HX.setflags(write=False)
    HZ.setflags(write=False)
    return Code(geo, HX, HZ)


def code_info(graph, L: int, code: Code | None = None) -> dict:
    """n, nE, L, N, m, n·L, k = N - rank H_X - rank H_Z, and girth (E3)."""
    if code is None:
        code = build_code(graph.n, graph.edges, L)
    geo = code.geometry
    k = geo.N - f2.rank(code.HX) - f2.rank(code.HZ)
    return {
        "name": graph.display_name(),
        "preset": graph.preset,
        "n": graph.n,
        "nE": geo.nE,
        "L": L,
        "N": geo.N,
        "m": geo.m,
        "x_checks": geo.mx,
        "k": int(k),
        "girth": girth(graph.n, graph.edges),
    }
