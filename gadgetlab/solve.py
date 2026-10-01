"""Constrained decode (§3.8), the single primitive: backends, memoisation and verify().

Works only with a check matrix H and integer indices of qubits and checks.
"""

from __future__ import annotations

import importlib.util
import os
import time
from dataclasses import dataclass

import numpy as np
import pulp

from gadgetlab import f2
from gadgetlab.errors import InputError

BACKENDS = ("ilp_cbc", "brute", "ilp_highs")
DEFAULT_BACKEND = "ilp_cbc"
BRUTE_MAX_KERNEL_DIM = 20


class SolverError(RuntimeError):
    """A solver returned a status other than Optimal or Infeasible, or an inconsistent answer."""


class BackendRefused(InputError):
    """The brute backend refuses instances whose solution space has dimension above 20."""


class BackendUnavailable(InputError):
    pass


class VerificationError(AssertionError):
    """verify() found a witness that does not satisfy its problem."""


def highs_available() -> bool:
    return importlib.util.find_spec("highspy") is not None


def available_backends() -> list[str]:
    return [b for b in BACKENDS if b != "ilp_highs" or highs_available()]


@dataclass(frozen=True)
class Problem:
    """Minimise sum(x_j, j in Q) s.t. x = 0 outside Q, x_j = phi(j) on T, and parity sigma_i on enforced checks."""

    H: np.ndarray
    sigma: np.ndarray
    Q: tuple[int, ...] | None  # None means global scope (Q = all qubits)
    phi: tuple[tuple[int, int], ...]  # (qubit, value), ascending by qubit


@dataclass(frozen=True)
class Result:
    status: str  # "optimal" or "infeasible"
    cost: int | None
    witness: np.ndarray | None  # 0/1 vector over all N qubits
    solver: str
    ms: float

    def ones(self) -> list[int]:
        return [] if self.witness is None else [int(j) for j in np.flatnonzero(self.witness)]


def verify(problem: Problem, result: Result) -> None:
    """Recompute every enforced parity check and the weight from the witness; raise VerificationError on any mismatch.

    Independent of the solver models: it works from H, sigma, Q and phi alone.
    """
    H = np.asarray(problem.H, dtype=np.int64)
    m, N = H.shape
    sigma = np.asarray(problem.sigma, dtype=np.int64) & 1
    if result.status == "infeasible":
        if result.cost is not None or result.witness is not None:
            raise VerificationError("an infeasible result must have cost null and no witness")
        return
    if result.status != "optimal":
        raise VerificationError(f"unknown status {result.status!r}")
    if result.witness is None or result.cost is None:
        raise VerificationError("an optimal result needs a cost and a witness")
    x = np.asarray(result.witness, dtype=np.int64)
    if x.shape != (N,) or not np.isin(x, (0, 1)).all():
        raise VerificationError("the witness must be a 0/1 vector over all qubits")
    if problem.Q is None:
        inside = np.ones(N, dtype=bool)
        enforced = np.arange(m)
    else:
        inside = np.zeros(N, dtype=bool)
        inside[list(problem.Q)] = True
        enforced = np.flatnonzero(H[:, inside].any(axis=1))
    outside_errors = np.flatnonzero(x * ~inside)
    if outside_errors.size:
        raise VerificationError(f"the witness has an error on qubit {int(outside_errors[0])}, outside the patch")
    for j, value in problem.phi:
        if x[j] != value:
            raise VerificationError(f"terminal qubit {j} should be {value} but the witness has {int(x[j])}")
    parities = (H[enforced] @ x) % 2
    bad = enforced[parities != sigma[enforced]]
    if bad.size:
        raise VerificationError(f"the witness violates enforced check {int(bad[0])}")
    weight = int(x[inside].sum())
    if weight != result.cost:
        raise VerificationError(f"the witness has weight {weight} but the reported cost is {result.cost}")


class Decoder:
    """Constrained decodes for one check matrix and one scope, memoised by
    (scope, patch, sigma restricted to the enforced checks, phi)."""

    def __init__(self, H, Q=None, backend: str = DEFAULT_BACKEND, check: bool = True):
        self.H = f2.as_f2(H)
        self.m, self.N = self.H.shape
        if backend not in BACKENDS:
            raise InputError(f"unknown solver {backend!r}; expected one of {', '.join(BACKENDS)}")
        if backend == "ilp_highs" and not highs_available():
            raise BackendUnavailable("the ilp_highs solver needs highspy (pip install highspy)")
        self.backend = backend
        self.check = check
        self.supports = [np.flatnonzero(self.H[i]) for i in range(self.m)]
        if Q is None:
            self.Q = None
            self.in_Q = np.ones(self.N, dtype=bool)
            self.enforced = np.arange(self.m)
        else:
            Q = sorted({int(j) for j in Q})
            if Q and not (0 <= Q[0] and Q[-1] < self.N):
                raise InputError("the patch contains a qubit index out of range")
            self.Q = tuple(Q)
            self.in_Q = np.zeros(self.N, dtype=bool)
            self.in_Q[list(Q)] = True
            self.enforced = np.flatnonzero(self.H[:, self.in_Q].any(axis=1))
        self.enforced_set = frozenset(int(i) for i in self.enforced)
        self.cache: dict = {}
        self.solver_calls = 0
        self.cache_hits = 0

    @property
    def scope(self) -> str:
        return "global" if self.Q is None else "local_closed"

    def normalise_phi(self, phi) -> tuple[tuple[int, int], ...]:
        items = phi.items() if isinstance(phi, dict) else phi
        out = {}
        for j, value in items:
            j, value = int(j), int(value)
            if not 0 <= j < self.N:
                raise InputError(f"terminal qubit {j} is out of range 0..{self.N - 1}")
            if value not in (0, 1):
                raise InputError(f"terminal qubit {j} must be set to 0 or 1, got {value}")
            if j in out:
                raise InputError(f"terminal qubit {j} is given twice")
            if not self.in_Q[j]:
                raise InputError(f"terminal qubit {j} lies outside the patch")
            out[j] = value
        return tuple(sorted(out.items()))

    def solve(self, sigma, phi=()) -> Result:
        sigma = f2.as_f2(sigma).reshape(-1)
        if sigma.shape[0] != self.m:
            raise InputError(f"the syndrome has length {sigma.shape[0]}, expected {self.m}")
        phi = self.normalise_phi(phi)
        key = (self.scope, self.Q, sigma[self.enforced].tobytes(), phi)
        hit = self.cache.get(key)
        if hit is not None:
            self.cache_hits += 1
            return hit
        self.solver_calls += 1
        start = time.perf_counter()
        if self.backend == "brute":
            status, x = self._solve_brute(sigma, phi)
        else:
            status, x = self._solve_ilp(sigma, phi)
        ms = (time.perf_counter() - start) * 1000.0
        if status == "optimal":
            result = Result("optimal", int(x[self.in_Q].sum()), x, self.backend, ms)
        else:
            result = Result("infeasible", None, None, self.backend, ms)
        if self.check:
            verify(Problem(self.H, sigma, self.Q, phi), result)
        self.cache[key] = result
        return result

    # ------------------------------------------------------------ shared reduction

    def _reduce(self, sigma, phi):
        """Free variables (Q minus terminals) and, per enforced check, the free positions in its
        support, the sum of terminal values in its support, and its support size."""
        fixed = dict(phi)
        free_mask = self.in_Q.copy()
        if fixed:
            free_mask[list(fixed)] = False
        free = np.flatnonzero(free_mask)
        position = np.full(self.N, -1, dtype=np.int64)
        position[free] = np.arange(free.size)
        rows = []
        for i in self.enforced:
            supp = self.supports[i]
            const = sum(fixed.get(int(j), 0) for j in supp)
            pos = position[supp]
            rows.append((int(i), pos[pos >= 0], const, supp.size))
        return free, rows

    def _witness(self, free, values, phi) -> np.ndarray:
        x = np.zeros(self.N, dtype=np.uint8)
        x[free] = values
        for j, value in phi:
            x[j] = value
        return x

    # ------------------------------------------------------------ ILP backends (§4.4)

    def _solve_ilp(self, sigma, phi):
        free, rows = self._reduce(sigma, phi)
        prob = pulp.LpProblem("constrained_decode", pulp.LpMinimize)
        x = [_variable(prob, f"x{int(j)}", 0, 1, "Binary") for j in free]
        prob += pulp.lpSum(x)
        for i, pos, const, size in rows:
            # (sum of x_j over supp(i)) - 2 y_i = sigma_i, with the terminal values moved to the right
            y = _variable(prob, f"y{i}", 0, size // 2, "Integer")
            prob += pulp.lpSum(x[p] for p in pos) - 2 * y == int(sigma[i]) - const, f"check{i}"
        prob.solve(_cbc() if self.backend == "ilp_cbc" else pulp.HiGHS(msg=False, threads=1))
        status = pulp.LpStatus.get(prob.status, str(prob.status))
        if status == "Infeasible":
            return "infeasible", None
        if status != "Optimal":
            raise SolverError(f"{self.backend} returned status {status!r}")
        values = np.zeros(free.size, dtype=np.uint8)
        for k, var in enumerate(x):
            v = var.varValue
            if v is None or abs(v - round(v)) > 1e-6 or round(v) not in (0, 1):
                raise SolverError(f"{self.backend} returned a non-binary value {v!r} for x{int(free[k])}")
            values[k] = int(round(v))
        objective = pulp.value(prob.objective) if x else 0
        if objective is None or int(round(objective)) != int(values.sum()):
            raise SolverError(f"{self.backend} reported objective {objective!r} for a witness of free weight {int(values.sum())}")
        return "optimal", self._witness(free, values, phi)

    # ------------------------------------------------------------ brute backend (§4.4)

    def _solve_brute(self, sigma, phi):
        free, rows = self._reduce(sigma, phi)
        A = np.zeros((len(rows), free.size), dtype=np.uint8)
        b = np.zeros(len(rows), dtype=np.uint8)
        for r, (i, pos, const, _size) in enumerate(rows):
            A[r, pos] = 1
            b[r] = (int(sigma[i]) + const) & 1
        x0 = f2.solve(A, b)
        if x0 is None:
            return "infeasible", None
        K = f2.kernel_basis(A)
        if K.shape[0] > BRUTE_MAX_KERNEL_DIM:
            raise BackendRefused(
                f"the brute solver refuses: the solution space has dimension {K.shape[0]}, above {BRUTE_MAX_KERNEL_DIM}")
        return "optimal", self._witness(free, min_weight_in_coset(x0, K), phi)


def _variable(prob, name: str, low: int, up: int, cat: str):
    add = getattr(prob, "add_variable", None)  # PuLP >= 3.3
    return add(name, low, up, cat) if add is not None else pulp.LpVariable(name, low, up, cat)


def _cbc():
    """CBC as bundled with PuLP (what PULP_CBC_CMD runs), else a cbc found on the PATH."""
    bundled = getattr(getattr(pulp, "PULP_CBC_CMD", None), "pulp_cbc_path", None)
    if bundled and os.path.isfile(bundled):
        return pulp.COIN_CMD(path=bundled, msg=False)
    return pulp.COIN_CMD(msg=False)


def min_weight_in_coset(x0: np.ndarray, K: np.ndarray) -> np.ndarray:
    """A minimum-weight vector of x0 + span(K), enumerating all 2^dim combinations.

    Combination index h * 2^lo + g selects basis vectors by the bits of g (the first lo
    vectors) and h (the rest); the first minimum in index order is returned.
    """
    d, n = K.shape
    if d == 0:
        return x0.copy()
    lo = min(d, 12)
    low = _span_table(K[:lo], n)
    high = _span_table(K[lo:], n)
    best_w, best = None, None
    for h in range(high.shape[0]):
        V = low ^ (x0 ^ high[h])
        w = V.sum(axis=1, dtype=np.int64)
        g = int(np.argmin(w))
        if best_w is None or w[g] < best_w:
            best_w, best = int(w[g]), V[g].copy()
    return best


def _span_table(B: np.ndarray, n: int) -> np.ndarray:
    """All 2^k combinations of the k rows of B; row g is the sum of the rows selected by the bits of g."""
    k = B.shape[0]
    T = np.zeros((1 << k, n), dtype=np.uint8)
    for b in range(k):
        T[1 << b: 1 << (b + 1)] = T[: 1 << b] ^ B[b]
    return T
