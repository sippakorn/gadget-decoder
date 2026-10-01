"""Linear algebra over F2 on 0/1 numpy arrays: rank, solving, kernels, span tests."""

from __future__ import annotations

import numpy as np


def as_f2(M) -> np.ndarray:
    """Copy of M reduced mod 2, as uint8."""
    return (np.asarray(M, dtype=np.int64) & 1).astype(np.uint8)


def row_reduce(M) -> tuple[np.ndarray, list[int]]:
    """Reduced row echelon form of M over F2.

    Returns (R, pivots): R holds the rank-many nonzero rows and pivots[r] is the
    pivot column of row r.
    """
    A = as_f2(M)
    if A.ndim != 2:
        raise ValueError("row_reduce expects a 2-D array")
    rows, cols = A.shape
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        below = np.flatnonzero(A[r:, c])
        if below.size == 0:
            continue
        p = r + int(below[0])
        if p != r:
            A[[r, p]] = A[[p, r]]
        hits = np.flatnonzero(A[:, c])
        hits = hits[hits != r]
        if hits.size:
            A[hits] ^= A[r]
        pivots.append(c)
        r += 1
    return A[:r], pivots


def rank(M) -> int:
    return len(row_reduce(M)[1])


def solve(A, b) -> np.ndarray | None:
    """One solution x of A x = b (mod 2), with every free variable 0; None if there is none."""
    A = as_f2(A)
    b = as_f2(b).reshape(-1)
    m, n = A.shape
    if b.shape[0] != m:
        raise ValueError("solve: b has the wrong length")
    R, pivots = row_reduce(np.hstack([A, b[:, None]]))
    if pivots and pivots[-1] == n:
        return None
    x = np.zeros(n, dtype=np.uint8)
    for row, c in enumerate(pivots):
        x[c] = R[row, n]
    return x


def kernel_basis(A) -> np.ndarray:
    """Basis of {x : A x = 0 (mod 2)} as the rows of a (dim × n) array."""
    A = as_f2(A)
    n = A.shape[1]
    R, pivots = row_reduce(A)
    pivot_set = set(pivots)
    free = [c for c in range(n) if c not in pivot_set]
    K = np.zeros((len(free), n), dtype=np.uint8)
    for i, f in enumerate(free):
        K[i, f] = 1
        for row, c in enumerate(pivots):
            K[i, c] = R[row, f]
    return K


def in_colspace(M, v) -> bool:
    return solve(M, v) is not None


def in_rowspace(M, v) -> bool:
    return solve(as_f2(M).T, v) is not None


class ColumnSpaceTest:
    """Fast repeated test of v in colspace(M), through a basis of the left kernel of M."""

    def __init__(self, M):
        self.Y = kernel_basis(as_f2(M).T)

    @property
    def full_row_rank(self) -> bool:
        return self.Y.shape[0] == 0

    def contains(self, v) -> bool:
        if self.full_row_rank:
            return True
        v = as_f2(v).reshape(-1)
        return not ((self.Y.astype(np.int64) @ v) & 1).any()
