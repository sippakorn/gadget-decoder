"""dX and dZ by exhaustive search up to weight 4, only for N <= 60 (E3)."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from gadgetlab import f2

WMAX = 4
NMAX = 60


def min_logical_weight(Hcheck, Hstab, wmax: int = WMAX) -> int | None:
    """Smallest weight of a vector in ker(Hcheck) outside rowspace(Hstab), searched up to wmax.

    Returns None when there is none of weight <= wmax. A vector v is outside rowspace(Hstab)
    exactly when K v != 0 for a basis K of ker(Hstab).
    """
    Hc = f2.as_f2(Hcheck)
    m, N = Hc.shape
    K = f2.kernel_basis(Hstab)
    # signature of qubit j: its syndrome bits, then its bits against the basis K
    sig = []
    for j in range(N):
        s = 0
        for i in np.flatnonzero(Hc[:, j]):
            s |= 1 << int(i)
        for r in np.flatnonzero(K[:, j]):
            s |= 1 << (m + int(r))
        sig.append(s)
    low = (1 << m) - 1
    by_syndrome: dict[int, list[int]] = defaultdict(list)
    for j, s in enumerate(sig):
        by_syndrome[s & low].append(j)

    def search(start: int, left: int, acc: int) -> bool:
        if left == 1:
            # the last qubit must cancel the syndrome of the others and leave a nonzero K-part
            for j in by_syndrome.get(acc & low, ()):
                if j >= start and ((acc ^ sig[j]) >> m) != 0:
                    return True
            return False
        for j in range(start, N - left + 1):
            if search(j + 1, left - 1, acc ^ sig[j]):
                return True
        return False

    for w in range(1, wmax + 1):
        if search(0, w, 0):
            return w
    return None


def distances(HX, HZ, wmax: int = WMAX, nmax: int = NMAX) -> dict:
    """{"dX", "dZ"}: an integer, ">4" when nothing is found up to weight 4, or "not computed" when N > 60."""
    N = np.asarray(HZ).shape[1]
    if N > nmax:
        return {"dX": "not computed", "dZ": "not computed"}

    def fmt(w):
        return f">{wmax}" if w is None else int(w)

    return {"dX": fmt(min_logical_weight(HZ, HX, wmax)), "dZ": fmt(min_logical_weight(HX, HZ, wmax))}
