# Usage: pip install pulp networkx numpy ; python3 warmup_check.py
# Exact MW decoding via ILP (CBC). Distances by exhaustive search up to weight 4 (small codes only).
"""
Plausibility check for the Day-1 warm-up:
  HGP(H1, R_L), H1 = vertex-edge incidence of a cubic graph G (i.e. the cycle code of G),
  R_L = repetition-code check matrix. Syndrome = all-ones on layer 1 of the Z-checks.
Claim: minimum-weight X-error decoding has OPT = 3n/2 - 2*alpha(G), for every L >= 2.
Also reports [[N, k, d]] so we can see what code family the hardness lives in.
"""
import itertools, numpy as np, networkx as nx, pulp

def gf2_rank(M):
    M = M.copy() % 2; r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return r

def in_rowspace(M, v, rankM):
    return gf2_rank(np.vstack([M, v])) == rankM

def hgp(H1, H2):
    r1, n1 = H1.shape; r2, n2 = H2.shape
    HX = np.hstack([np.kron(H1, np.eye(n2, dtype=int)), np.kron(np.eye(r1, dtype=int), H2.T)]) % 2
    HZ = np.hstack([np.kron(np.eye(n1, dtype=int), H2), np.kron(H1.T, np.eye(r2, dtype=int))]) % 2
    return HX, HZ

def rep(L):
    R = np.zeros((L-1, L), dtype=int)
    for i in range(L-1): R[i, i] = R[i, i+1] = 1
    return R

def min_logical_weight(Hcheck, Hstab, wmax):
    """smallest-weight vector in ker(Hcheck) minus rowspace(Hstab), searched up to weight wmax"""
    N = Hcheck.shape[1]; rS = gf2_rank(Hstab)
    for w in range(1, wmax+1):
        for supp in itertools.combinations(range(N), w):
            v = np.zeros(N, dtype=int); v[list(supp)] = 1
            if not (Hcheck @ v % 2).any() and not in_rowspace(Hstab, v, rS):
                return w
    return f">{wmax}"

def mw_decode(H, s):
    m, N = H.shape
    prob = pulp.LpProblem("mw", pulp.LpMinimize)
    x = [pulp.LpVariable(f"x{j}", cat="Binary") for j in range(N)]
    y = [pulp.LpVariable(f"y{i}", lowBound=0, cat="Integer") for i in range(m)]
    prob += pulp.lpSum(x)
    for i in range(m):
        prob += pulp.lpSum(x[j] for j in range(N) if H[i, j]) - 2*y[i] == int(s[i])
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    return int(round(pulp.value(prob.objective)))

def alpha(G):
    return max(len(c) for c in nx.find_cliques(nx.complement(G)))

rng = np.random.default_rng(0)
print(f"{'n':>3} {'L':>2} {'girth':>5} | {'N':>4} {'k':>3} {'dX':>4} {'dZ':>4} | {'OPT':>4} {'3n/2-2a':>7} match")
for n in [6, 8, 10, 12]:
    for trial in range(2):
        G = nx.random_regular_graph(3, n, seed=int(rng.integers(1e6)))
        if not nx.is_connected(G): continue
        E = list(G.edges()); V = list(G.nodes())
        H1 = np.zeros((len(V), len(E)), dtype=int)       # rows = vertices, cols = edges
        for j, (u, v) in enumerate(E): H1[u, j] = H1[v, j] = 1
        g = nx.girth(G); a = alpha(G)
        for L in [2, 3, 4]:
            HX, HZ = hgp(H1, rep(L))
            assert not (HX @ HZ.T % 2).any()
            N = HX.shape[1]; k = N - gf2_rank(HX) - gf2_rank(HZ)
            s = np.zeros(HZ.shape[0], dtype=int)
            s[[e*(L-1) + 0 for e in range(len(E))]] = 1     # all-ones on layer 1
            opt = mw_decode(HZ, s)
            dX = min_logical_weight(HZ, HX, 4) if N <= 60 else "-"
            dZ = min_logical_weight(HX, HZ, 4) if N <= 60 else "-"
            pred = 3*n//2 - 2*a
            print(f"{n:>3} {L:>2} {g:>5} | {N:>4} {k:>3} {str(dX):>4} {str(dZ):>4} | {opt:>4} {pred:>7} {opt==pred}")
