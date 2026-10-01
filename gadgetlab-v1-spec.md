# GadgetLab v1 — specification

A syndrome-gadget explorer for the quantum code HGP(cycle code of G, R_L).

- Status: scope approved by the project owner on 2026-10-01.
- Audience: the agent implementing v1.
- Hand-off bundle: this file plus `warmup_check.py` (reference code; §10 notes one issue in it).

## 0. Read this first

- Build v1 as specified. Requirements have IDs (E = engine, V = viewer, T = tests, A = acceptance); report progress against them.
- §3 is normative mathematics. Do not change definitions, index formulas or conventions. If something there looks wrong or ambiguous, stop and ask.
- Correctness beats features. The tool produces evidence for research lemmas, so a wrong number is worse than a missing button. If the cross-check (T3) or the regression tests (T2, T5) fail, stop and report; never edit an expected value to make a test pass.
- Ask before adding dependencies beyond §4.3, changing the JSON formats in §5, or dropping a requirement.
- Keep research content out of code, comments and docs beyond what this spec states. In particular, do not add explanations of why any cost formula holds.

## 1. Purpose and vocabulary

A research project studies how hard it is to find the lightest error that explains a given syndrome in a family of quantum codes. One planned method is a computer search for gadgets: small syndrome patterns whose decoding cost behaves like a logic component, such as a wire that carries a bit or a branch that copies it.

v1 is a sandbox for building intuition on small instances of one code family. It must be structured so it can later grow into the real search tool (§11). There is one user, working on a laptop. It is not a product: no accounts, no deployment, no telemetry.

| Term | Meaning in this spec |
|---|---|
| Qubit | One column of the check matrix H_Z; one position of an error vector |
| Check | One row of H_Z (a Z-check); it fires when an odd number of its qubits are in error |
| Error | A 0/1 vector over the qubits; its weight is its number of 1s |
| Syndrome | The 0/1 vector of fired checks, σ = H_Z·x mod 2 |
| Decode | Find a minimum-weight error with a given syndrome |
| Terminal | A qubit whose value the user fixes when building a cost table |
| Cost table | The decode cost for every setting of the terminals |
| Spec | The required shape of a cost table (§3.11) |
| Region | The checks whose syndrome bits a search varies |
| Patch | A set of qubits outside which the error is forced to be zero (local scope) |
| Margin | How far the penalised rows of a table exceed the minimum the spec demands |

## 2. Scope

### 2.1 In scope (v1)

1. Build HGP(cycle code of G, R_L) for a small cubic graph G and L layers, and show its parameters.
2. Click checks to set a syndrome, decode it, and draw the optimal error.
3. Pick terminal qubits and compute the cost table.
4. Write a target spec and search syndrome patterns on a set of checks for ones that meet it.
5. Two evaluation scopes: global (the whole code) and local_closed (a patch, with no error outside it).
6. Sessions saved and loaded as JSON, and a CLI that can compute everything the viewer can.

### 2.2 Out of scope (v1)

- Any code other than HGP(H1, R_L) with H1 the incidence matrix of a connected cubic graph.
- Z-error decoding, depolarizing or noisy-syndrome models, maximum-likelihood decoding.
- Robustness against arbitrary errors outside a patch (v2, §11).
- Symmetry reduction, SAT or MaxSAT backends, warm-started solvers.
- Phones, multiple users, remote hosting.

## 3. Mathematics (normative)

### 3.1 Conventions

- F2 = {0, 1} with arithmetic mod 2. Vectors are 0/1 arrays. The weight |x| is the number of 1s.
- n is the number of vertices of G, not the code length. nE = 3n/2 is the number of edges. The code length is N (§3.4).

### 3.2 The graph G

- A connected cubic graph (every vertex has degree 3) with vertices 0, …, n−1; n is even and at least 4.
- Edges are stored canonically: each edge as (u, w) with u < w, and the list sorted lexicographically. The edge index of an edge is its position in this list.
- girth(G) is the length of a shortest cycle. α(G) is the size of a largest independent set; it is used only in T5 and behind the flag in V12.

### 3.3 Classical matrices

- H1 is n × nE, with H1[v, e] = 1 exactly when v is an endpoint of edge e.
- R_L is (L−1) × L, with R_L[i, i] = R_L[i, i+1] = 1 and 0 elsewhere. The engine accepts 2 ≤ L ≤ 8; the viewer offers 2 ≤ L ≤ 5.

### 3.4 The quantum code

With ⊗ the Kronecker product, exactly as `hgp(H1, rep(L))` in warmup_check.py:

- H_X = [ H1 ⊗ I_L  |  I_n ⊗ R_Lᵀ ]
- H_Z = [ I_nE ⊗ R_L  |  H1ᵀ ⊗ I_(L−1) ]

Sizes: N = nE·L + n·(L−1) qubits (columns), m = nE·(L−1) Z-checks (rows of H_Z), and n·L X-checks (rows of H_X).

### 3.5 Index map

| Object | Typed address | Ranges | Index |
|---|---|---|---|
| Edge-qubit | (E, e, l) | 0 ≤ e < nE, 0 ≤ l ≤ L−1 | e·L + l |
| Vertex-qubit | (V, v, t) | 0 ≤ v < n, 0 ≤ t ≤ L−2 | nE·L + v·(L−1) + t |
| Z-check | (e, t) | 0 ≤ e < nE, 0 ≤ t ≤ L−2 | e·(L−1) + t |
| X-check | (v, l) | 0 ≤ v < n, 0 ≤ l ≤ L−1 | v·L + l |

Supports, for an edge e = (u, w):

- Z-check (e, t) contains exactly the four qubits (E, e, t), (E, e, t+1), (V, u, t) and (V, w, t).
- X-check (v, l) contains the edge-qubits (E, e, l) for the three edges e at v, and the vertex-qubits (V, v, t) for t in {l−1, l} with 0 ≤ t ≤ L−2.

The engine builds H_Z and H_X both ways (Kronecker formula and explicit supports) and asserts that they are equal (T1).

### 3.6 Geometric picture (what the viewer draws)

Take L copies of G, called layers 0 to L−1, and join copy t to copy t+1 by a vertical edge at every vertex. This is the product graph G □ P_L.

- Edge-qubit (E, e, l) is edge e in layer l.
- Vertex-qubit (V, v, t) is the vertical edge at v between layers t and t+1.
- Z-check (e, t) is the square face bounded by edge e in layers t and t+1 and by the vertical edges at its endpoints u and w.

So an error is a set of edges of G □ P_L, and its syndrome is the set of square faces with an odd number of error edges on their boundary.

### 3.7 Syndrome and decoding

- The syndrome of an error x is σ = H_Z·x (mod 2), a vector of length m.
- MW(σ) = min { |x| : H_Z·x = σ }. The problem is infeasible when σ is not in the column space of H_Z.

### 3.8 Constrained decode, the single primitive

Inputs:

- a syndrome σ;
- a patch Q, which is a set of qubits (in global scope, Q is all qubits);
- a terminal assignment φ, giving a value in {0, 1} to each qubit of a set T ⊆ Q (T may be empty).

Enforced checks: in global scope, every Z-check. In local_closed scope, C_Q = the Z-checks whose support contains at least one qubit of Q.

Problem: minimise the sum of x_j over j in Q, subject to

- x_j in {0, 1} for every qubit j;
- x_j = 0 for every j not in Q;
- x_j = φ(j) for every j in T;
- for every enforced check i, the sum of x_j over j in supp(i) ≡ σ_i (mod 2).

Output: status (optimal or infeasible), cost (an integer, or null if infeasible), and witness x (null if infeasible).

Scopes: "global" (Q = all qubits) and "local_closed" (Q = a patch from §3.9, meaning "assume no error outside the patch").

### 3.9 Patch (v1)

Given a vertex set U and a layer range [l0, l1] with 0 ≤ l0 ≤ l1 ≤ L−1, the patch Q(U, l0, l1) consists of

- every edge-qubit (E, e, l) where e has at least one endpoint in U and l0 ≤ l ≤ l1, and
- every vertex-qubit (V, v, t) where v is in U and l0 ≤ t ≤ l1−1.

Helper: ball(s, r) is the set of vertices at graph distance at most r from s in G.

### 3.10 Cost table

- Terminals T = (t_1, …, t_k) is an ordered list of distinct qubits, with 1 ≤ k ≤ 6 (the viewer warns above 4). In local_closed scope every terminal must lie in Q.
- Rows: every a in {0,1}^k, written as the bit string a_1 a_2 … a_k, sorted lexicographically (all zeros first).
- Cost(a) is the constrained-decode cost with φ(t_i) = a_i, or infinity (JSON null) if infeasible.
- Runtime invariants (a violation is an error, not a warning): (i) the minimum of Cost(a) over all rows equals the cost with no terminals fixed; (ii) the row matching the terminal values of the unconstrained witness has that same cost.

### 3.11 Spec and matching

A spec gives each row a exactly one entry:

- ["eq", k]: Cost(a) = c + k;
- ["ge", k]: Cost(a) ≥ c + k (infinity satisfies this);
- ["any"]: no condition.

Here k is a non-negative integer, and at least one row must be "eq". A table matches the spec if a single integer c makes every entry true. c is read off the eq rows as Cost(a) − k; every eq row must have finite cost and give the same c.

Margin = the minimum, over ge rows with finite cost, of Cost(a) − c − k. If there is no such row, the margin is infinity (JSON null).

Presets (illustrative only; whether any pattern matches them is exactly what the tool explores):

- wire (k = 2): 00 eq 0, 11 eq 0, 01 ge 1, 10 ge 1;
- copy3 (k = 3): 000 eq 0, 111 eq 0, every other row ge 1;
- custom: every row "any" until edited.

### 3.12 Search

Inputs: a base syndrome σ0; a region R, which is an ordered list of Z-checks; terminals T; a spec S; a scope (and patch); a budget B (default 4096); a seed (default 0); top_k (default 50); and an optional max_weight.

Candidates are patterns p in {0,1}^|R|. The candidate syndrome is σ0 with its entries on R replaced by p.

- If 2^|R| ≤ B, take every pattern, in lexicographic order. Otherwise take B distinct patterns sampled uniformly with the seed.
- Skip patterns with more than max_weight ones (counted as filtered).
- In global scope, skip candidate syndromes outside the column space of H_Z (counted as infeasible_syndrome). For this code family H_Z has full row rank, so the count will be 0; the filter exists for the codes of v2.
- In local_closed scope, drop region checks outside C_Q, with a warning.

Evaluation is lazy: compute the eq rows first and fix c, rejecting at the first inconsistency; then compute the ge rows, rejecting at the first violation. For candidates that match, also compute the "any" rows so that the stored table is complete.

Ranking of matches: margin descending (infinity first), then c ascending, then |p| ascending, then p lexicographically. Keep the top_k, each with its full table and one witness per row.

Summary counts: candidates, filtered, infeasible_syndrome, rejected, matched, solver_calls, cache_hits, elapsed_s, and exhaustive (true or false).

## 4. Architecture

### 4.1 Shape

Hybrid: a Python engine (package, CLI and local HTTP server) plus a browser viewer (one static HTML file).

- Live mode: `gadgetlab serve` serves the viewer at http://127.0.0.1:8765 together with a JSON API, and the viewer computes through the API.
- Static mode: opening `viewer/index.html` from disk lets the user load a saved session and inspect it. Compute controls are disabled, with a one-line explanation.

### 4.2 Repository layout

```text
gadgetlab/
  pyproject.toml
  README.md               install steps and a 10-minute walkthrough
  gadgetlab/
    f2.py                 GF(2) rank, column-space test, solving, kernel basis
    graphs.py             presets, canonical edges, layouts, ball()
    geometry.py           typed addresses <-> indices, supports, panel data (adapter)
    code.py               H1, R_L, hgp(), code info (N, m, k, girth)
    distance.py           dX and dZ by exhaustive search up to weight 4
    solve.py              constrained decode, backends, verify()
    table.py              cost tables and their invariants
    spec.py               spec model, presets, matching, margin
    search.py             candidates, lazy evaluation, cache, process pool, jobs
    session.py            JSON schema, load, save, validation
    server.py             HTTP API and static viewer
    cli.py                command-line entry points
  viewer/index.html       HTML, CSS and JS inline; no build step; no CDN
  tests/
  examples/k4_L2.session.json
  examples/petersen_L3.session.json
```

Separation rule, needed for v2: `solve`, `table`, `spec` and `search` work only with matrices and integer indices. They must not import `graphs` or `geometry`, and must not know what an edge or a vertex is. `geometry` translates typed addresses at the boundary (session, server, CLI).

### 4.3 Dependencies

- Python 3.10 or later with numpy, networkx and pulp (which bundles CBC). Optional: highspy. Development: pytest.
- Server: standard library only (`http.server.ThreadingHTTPServer`, `json`, `multiprocessing`).
- Viewer: vanilla JavaScript and SVG. No frameworks, no external fonts or scripts, and no network requests other than to the local API.

### 4.4 Solver backends

All backends return the same result type.

- `ilp_cbc` (default). Binary x_j; for each enforced check i an integer slack y_i with 0 ≤ y_i ≤ floor(|supp(i)|/2) and (sum of x_j over supp(i)) − 2·y_i = σ_i, as in warmup's `mw_decode`. Read the objective only when the solver status is Optimal; map Infeasible to status "infeasible"; treat any other status as an error.
- `brute` (exact, small instances; the reference for T3). Substitute the fixed variables, drop the variables outside Q, then find a particular solution and a kernel basis by Gaussian elimination. Enumerate the kernel if its dimension is at most 20; otherwise refuse.
- `ilp_highs` (optional, when highspy is installed).

`verify(result)` recomputes every enforced parity check and the weight from the witness, using code independent of the solver's model. Tests call it on every result; the app calls it on every result while debug is on (the default in v1).

### 4.5 Jobs, cache and parallelism

- A search runs in a process pool (default workers: max(1, cpu_count − 1)). Progress is reported as (done, total), and a cancel flag is checked between candidates.
- Within one job, constrained decodes are memoised by (scope, patch, σ restricted to the enforced checks, φ).
- The CLI prints a progress line; the API exposes progress by polling.

## 5. Data formats (JSON, normative)

Typed addresses:

- edge-qubit: {"q": "E", "e": 3, "l": 0}
- vertex-qubit: {"q": "V", "v": 2, "t": 0}
- Z-check: {"e": 3, "t": 0}

JSON has no infinity. An infeasible cost is "cost": null together with "status": "infeasible". An infinite margin is "margin": null.

Session file (schema "gadgetlab.session/1"; lists shortened):

```json
{
  "schema": "gadgetlab.session/1",
  "graph": {"preset": "petersen", "seed": null, "n": 10,
            "edges": [[0, 1], [0, 4], [0, 5]],
            "layout": {"0": [0.5, 0.0], "1": [0.98, 0.35]}},
  "L": 3,
  "syndrome": [{"e": 0, "t": 0}],
  "terminals": [{"q": "E", "e": 3, "l": 1}, {"q": "V", "v": 2, "t": 0}],
  "scope": "global",
  "patch": null,
  "region": [{"e": 5, "t": 1}],
  "spec": {"preset": "wire",
           "rows": {"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1], "11": ["eq", 0]}},
  "results": {"decode": null, "table": null, "search": null},
  "notes": ""
}
```

When set, "patch" is {"vertices": [0, 1, 4], "layers": [0, 2]}.

Results:

- decode: {"status", "cost", "error": [addresses], "verified", "solver", "ms", "inputs_hash"}
- table: {"terminals", "rows": [{"a": "01", "status", "cost", "error"}], "unconstrained_cost", "match": {"matches", "c", "margin"}, "inputs_hash"}
- search: {"params", "summary", "matches": [{"rank", "pattern", "c", "margin", "table"}], "inputs_hash"}

"pattern" is aligned with the order of "region". inputs_hash is a hash of everything a result depends on; the viewer uses it to detect stale results (V11).

## 6. Engine requirements

- E1 Build. Presets: k4, k33, prism, cube, petersen, heawood and random_cubic(n, seed), via networkx complete_graph(4), complete_bipartite_graph(3, 3), circular_ladder_graph(3), hypercube_graph(3) relabelled to integers, petersen_graph(), heawood_graph() and random_regular_graph(3, n, seed). Also accept a user edge list. Reject graphs that are not connected or not cubic, with a clear message. For random_cubic, if the graph is disconnected retry with seed + 1 (up to 100 times) and record the seed actually used. Layouts are deterministic and normalised to [0, 1]: standard drawings for named presets (e.g. Petersen as an outer pentagon and inner pentagram), otherwise a spring layout with seed 0.
- E2 Index map and supports exactly as §3.5, with both constructions and the equality assertion.
- E3 Code info: n, nE, L, N, m, k = N − rank H_X − rank H_Z, and girth. On request, dX (minimum weight in ker H_Z outside the row space of H_X) and dZ (the same with X and Z swapped), by exhaustive search up to weight 4 and only when N ≤ 60. Report ">4" when nothing is found up to weight 4, and "not computed" when N > 60.
- E4 Constrained decode (§3.8) in both scopes, with every available backend and verify().
- E5 Cost tables (§3.10) with the invariant checks.
- E6 Spec matching, margin and presets (§3.11).
- E7 Search (§3.12): budget, sampling, filters, lazy evaluation, ranking, cache, pool, progress and cancel.
- E8 Sessions: load, save and validate against §5. Reject out-of-range addresses with a message naming the field.
- E9 CLI: `gadgetlab build | decode | table | search | distances | serve`. `build` creates a session file from a preset or an edge list; the other commands read a session file and write their results back into it (or to --out). Output is deterministic for a fixed seed, and the CLI and the API produce identical result JSON for the same session.
- E10 Server: binds to 127.0.0.1 only; port 8765 by default (--port). Flags: --r0-preset (V12) and --debug. Endpoints:

| Method and path | Body | Returns |
|---|---|---|
| GET / | none | the viewer |
| GET /api/config | none | {"r0_preset": false, "debug": true} |
| POST /api/build | {graph, L} | a session skeleton with code info and layout |
| POST /api/decode | session | decode result |
| POST /api/table | session | table result |
| POST /api/distances | session | {"dX", "dZ"} |
| POST /api/random_syndrome | session plus {weight, seed} | {"syndrome"}: H_Z·x for a random x of that weight |
| POST /api/search | session plus {budget, seed, top_k, max_weight, workers} | {"job_id"} |
| GET /api/jobs/ID | none | {"state", "progress": {"done", "total"}, "summary", "result"} |
| POST /api/jobs/ID/cancel | none | {"state"} |

Job states: running, done, cancelled, error. Invalid input returns HTTP 400 with {"error": message}; tracebacks appear only with --debug.

## 7. Viewer requirements

### 7.1 Functional

- V1 Page structure: a top bar with code info (preset, n, L, N, m, k, girth, and dX/dZ with a "Compute distances" button that is disabled with an explanation when N > 60); a left column of controls (graph, L, scope and patch, modes); the canvas in the centre; a right column for results (decode, cost table, spec, search); and a legend that is always visible.
- V2 The canvas is a sequence of panels, Layer 0, Gap 0, Layer 1, Gap 1, …, Layer L−1, each drawing G with the same layout.
  - In layer panel l, edges are the edge-qubits (E, e, l); vertices are plain dots, not interactive in v1.
  - In gap panel t, edges are the Z-checks (e, t); vertices are the vertex-qubits (V, v, t).
  - Each panel title says what its marks are, e.g. "Gap 0: edges are checks, vertices are qubits". Panels sit in one horizontally scrolling row, or wrap into a grid; the user can switch.
  - Optional: drag vertices to tidy the layout; the layout is saved in the session and shared by all panels.
- V3 The visual encoding never relies on colour alone:

| Thing | Encoding |
|---|---|
| Check with σ = 1 | thick solid stroke in the syndrome colour |
| Check with σ = 0 | thin stroke in the rule colour |
| Error on an edge-qubit | thick stroke in the error colour, in its layer panel |
| Error on a vertex-qubit | enlarged filled dot in the error colour, in its gap panel |
| Terminal | numbered badge (1, 2, …) beside the qubit |
| Region check | dashed halo around the edge |
| Patch | light shading under the included qubits |

- V4 Modes (a toolbar with one mode active; keys 1 to 4):
  - Syndrome: clicking a gap-panel edge toggles σ for that Z-check.
  - Terminals: clicking any qubit adds it to T or removes it. The order is shown and can be changed; at most 6.
  - Region: clicking a gap-panel edge toggles it in R.
  - Patch: clicking a vertex in any panel toggles it in U; a control sets the layer range; a helper fills U with ball(seed, radius).
  - Every mark has a hover tooltip with its typed address and endpoints, e.g. "Z-check e = 3 (1–4), t = 0".
  - Syndrome tools: "Clear syndrome", and "Random syndrome" (from a random error of chosen weight, seeded, via the API).
- V5 Decode: the "Decode" button shows status, cost, solver and time, and draws the witness. The viewer then recomputes the witness's syndrome itself, from the §3.5 supports, and shows either "Syndrome verified" or a warning listing the mismatched checks. Keep this code separate from the engine: it is an independent second implementation.
- V6 Cost table: lists the terminals; one row per assignment showing the bit string, the cost ("—" when infeasible), the spec entry and a pass or fail mark; a "Show" button draws that row's witness. Shows c and the margin when the spec is valid.
- V7 Spec editor: a preset menu (wire, copy3, custom) and, per row, an eq / ge / any selector with an integer offset. Inline validation, e.g. "Add at least one eq row". Changing the number of terminals resets the spec to custom, with a notice.
- V8 Search: shows |R|, 2^|R| and whether the run will be exhaustive or sampled; inputs for budget, seed, top_k, max_weight and workers; "Run search" and "Cancel search"; a progress bar with live counts; a results list (rank, c, margin, |p|). Clicking a result loads its syndrome onto the canvas, labelled e.g. "From search result 3", and loads its cost table.
- V9 Sessions: "New session", "Open session" (file picker) and "Save session" (download JSON). Undo and redo for edits are recommended (last 50 edits).
- V10 Static mode (file://): load a session JSON and display everything; compute buttons are disabled with a one-line explanation.
- V11 Staleness: any edit to the syndrome, terminals, scope, patch, region or spec marks the shown results as stale (dimmed and labelled), using inputs_hash. A stale result is never shown as current.
- V12 Spoiler guard: a flag r0_preset, default false (server flag --r0-preset; in static mode the URL parameter ?r0_preset=1). While it is false, no control fills a whole gap with σ = 1 and no α(G) comparison is shown. When it is true, a "Fill gap 0" preset appears together with a panel comparing MW(σ) with 3n/2 − 2α(G). Manual clicking is never restricted.
- V13 Self-test: opening index.html?selftest=1 runs the T4 vectors through the viewer's own support code and prints pass or fail (no Node needed).
- V14 Performance and platforms: interactions respond within 50 ms for n ≤ 20 and L ≤ 5. Current desktop Chrome, Firefox and Safari. Designed for 1280 px and wider, usable down to 1024 px; phones are not a target.

### 7.2 Visual direction

The subject is a lab instrument for decoding experiments on a layered graph. The canvas is the instrument; everything around it stays quiet.

Palette for the light theme. Derive the dark theme from the same hues, follow prefers-color-scheme, and keep WCAG AA contrast in both.

| Token | Hex | Use |
|---|---|---|
| ground | #F3F5F7 | page background |
| ink | #1D2733 | text and graph vertices |
| rule | #98A3AF | idle edges and dividers |
| syndrome | #B5176F | fired checks |
| error | #0E7A86 | error qubits |
| patch | #DCE7F3 | patch shading |

Type: the platform's system UI font (the viewer must work offline), with tabular numerals for every index and cost. Scale 13 / 15 / 18 / 24 px; sentence case throughout.

Layout:

```text
+--------------------------------------------------------------------------+
| Petersen, L = 3     n 10   N 65   checks 30   k 6   girth 5   [Compute]  |
+------------+-------------------------------------------------+-----------+
| Graph      |  Layer 0 | Gap 0 | Layer 1 | Gap 1 | Layer 2    | Decode    |
| Layers     |                                                 | Table     |
| Scope      |   panels scroll horizontally or wrap            | Spec      |
| Modes 1-4  |                                                 | Search    |
+------------+-------------------------------------------------+-----------+
| Legend                                                                   |
+--------------------------------------------------------------------------+
```

Principles:

- Every mark carries meaning; no decorative cards, shadows or gradients.
- Controls use the researcher's words: checks, qubits, syndrome, terminals. Buttons say what they do ("Decode", "Compute table", "Run search").
- Motion only as feedback to an action, e.g. one brief highlight when a search result is loaded; respect reduced-motion settings.
- Avoid template defaults: all-caps labels, a label above every block, middle-dot meta strings, arrows appended to button text, monospace for small labels.
- Visible keyboard focus on every control.

## 8. Tests

- T1 Construction: the Kronecker and explicit-support matrices are equal; H_X·H_Zᵀ = 0 (mod 2); the sizes N, m and n·L are correct. All presets, L in {2, 3, 4}.
- T2 Code parameters (all presets, plus random_cubic for n in {6, 8, 10, 12} with two seeds each): k = n/2 + 1, and girth matches networkx. Where N ≤ 60, dX and dZ match warmup's `min_logical_weight`, and min(dX, dZ) = min(L, girth) whenever that value is at most 4. Record dX and dZ separately in the test output. A failure here contradicts a recorded research claim: report it and do not change the expectation. Mark cases with N > 40 as slow.
- T3 Solver cross-check: for K4 with L = 2 and L = 3, and K3,3 with L = 2, generate 200 syndromes each (150 as H_Z·x for a random x of weight 1 to 6, and 50 uniformly random). For each, try random terminals (k in {0, 1, 2}) with random values, in both scopes (random patches for local_closed). ilp_cbc and brute must agree exactly on status and cost, and verify() must pass on every witness.
- T4 Single-qubit syndromes, for every qubit of K4 with L = 3: an error on (E, e, l) fires exactly the checks (e, t) with t in {l−1, l} and 0 ≤ t ≤ L−2; an error on (V, v, t) fires exactly the checks (e, t) for the three edges e at v.
- T5 Regression of a recorded claim: with σ = 1 on every Z-check (e, 0) and 0 elsewhere, MW(σ) = 3n/2 − 2α(G). Instances: every preset with L in {2, 3, 4}, plus random_cubic for n in {6, 8, 10, 12} with two seeds each; global scope, ilp_cbc. A failure is reported, never "fixed".
- T6 Table invariants (§3.10) on random cases.
- T7 Spec matching on hand-made tables: null costs; inconsistent eq rows; no ge rows (infinite margin); an infeasible ge row.
- T8 Search: in exhaustive mode the candidate count is 2^|R|; sampled patterns are distinct, and identical for equal seeds; ranking follows §3.12; cancel stops a job within 2 seconds; lazy evaluation finds the same matches as full evaluation (small case, |R| ≤ 8).
- T9 Sessions: saving then loading gives an equal session; invalid addresses are rejected.
- T10 API smoke test: start the server on a free port and call every endpoint once.
- T11 Infeasibility: build cases that must be infeasible (e.g. a check with σ = 1 whose four qubits are all terminals set to 0, and local_closed cases with inconsistent parities) and check status infeasible, null cost and no exception. Note that in global scope with no terminals every syndrome is feasible for this family, because H_Z has full row rank.
- T12 Determinism: two CLI runs with the same seed and workers = 1 give identical output; with more workers, the matches and their ranking are identical.

## 9. Acceptance (definition of done)

- A1 `pip install -e .` followed by `pytest` passes; slow tests may be skipped by default but must pass when run.
- A2 Run `gadgetlab serve`, open the browser, build Petersen with L = 3, click three checks and decode: the result arrives in under a second, the witness is drawn, and "Syndrome verified" appears.
- A3 Choose two terminals: the cost table has four rows, and its minimum equals the decode cost.
- A4 With the wire spec and a region of 10 checks, the exhaustive search (1024 candidates) finishes within 5 minutes on a 4-core laptop, with live progress and a working "Cancel search"; clicking a result loads it.
- A5 Save the session, open the viewer from disk (static mode) and load the session: everything displays.
- A6 The CLI reproduces the decode, table and search JSON of that saved session.
- A7 The README walkthrough reproduces A2 to A6 in about 10 minutes.

## 10. Pitfalls

- `mw_decode` in warmup_check.py reads the objective without checking the solver status. That is harmless there (its syndromes are feasible), but here terminals and patches often make a problem infeasible, so always check the status (§4.4).
- The index formulas follow Kronecker-product order (§3.5); T1 guards them.
- n counts vertices, not the code length.
- warmup_check.py uses networkx's edge order, while this tool uses canonical order (§3.2). Costs do not depend on edge order; addresses do.
- In local_closed scope, checks whose whole support lies outside Q are not enforced, and a terminal outside Q is an input error.
- Use null for infinity, as in §5.
- Build PuLP problems inside the pool workers and never share solver objects between processes. CBC runs as a subprocess, so per-call overhead (roughly 10 to 50 ms) dominates small problems; that is expected in v1.
- Never display a result computed for different inputs (V11).

## 11. Roadmap for v2 (do not build, but do not block)

- Other codes: HGP(H1, H2) from arbitrary matrices, including quantum expander codes HGP(C, C). When H1 has a column of weight above 2 there is no graph G, so the viewer will need a Tanner-graph view; this is why geometry is an adapter (§4.2).
- Robust local mode: the spec must hold for every setting of the qubits just outside the patch. This is exponential in the boundary size; plan a propose-and-refute loop in which the solver finds a breaking boundary setting, which is added as a constraint before the next proposal.
- Symmetry reduction using automorphisms of G that fix the terminals.
- Faster backends: warm-started HiGHS, MaxSAT or SAT.
- Noise variants: joint X/Z (depolarizing) and noisy-syndrome (phenomenological) decoding.
- Larger expander instances through the CLI only, including explicit high-girth cubic graphs.
- Export of a single result as a standalone, shareable HTML page.

## 12. Report back

When done, deliver: (1) a checklist of every requirement ID with its status; (2) the test output; (3) anything not met, and why; (4) every place where §3 seemed wrong or ambiguous, as an exact question. Do not resolve research-level questions yourself.
