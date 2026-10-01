# GadgetLab v1

A syndrome-gadget explorer for the quantum code HGP(cycle code of G, R_L), where G is a small
connected cubic graph and L is the number of layers. It builds the code, decodes syndromes,
computes cost tables for chosen terminal qubits, and searches syndrome patterns on a region of
checks for ones whose cost table meets a spec. The specification is
[gadgetlab-v1-spec.md](gadgetlab-v1-spec.md).

The tool has two parts: a Python engine (package, CLI and a local HTTP server) and a browser
viewer, `viewer/index.html`, a single file with no build step and no network access beyond the
local engine.

For a guided tour of the mathematics and of gadget searches, with worked examples, see
[tutorial.md](tutorial.md).

## Install

Python 3.10 or later.

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -e ".[dev]"          # numpy, networkx, pulp (bundles the CBC solver), pytest
pip install -e ".[dev,highs]"    # optional: also the ilp_highs solver backend
pytest                           # about 30 s; slow cases are skipped
pytest --runslow                 # everything
```

## Ten-minute walkthrough

This reproduces acceptance checks A2 to A6 of the spec.

**1. Start the engine and the viewer.**

```sh
gadgetlab serve
```

Open http://127.0.0.1:8765. The viewer builds Petersen with L = 3 on start. The top bar shows
the code: n 10, N 65, checks 30, k 6, girth 5. (To build another graph, pick it under Graph and
Layers L and press New session.)

The canvas shows the panels Layer 0, Gap 0, Layer 1, Gap 1, Layer 2. In a layer panel the edges
are edge-qubits. In a gap panel the edges are Z-checks and the vertices are vertex-qubits. Hover
over any mark for its address, for example "Z-check e = 3 (1–4), t = 0".

**2. Decode (A2).** In Syndrome mode (key 1), click three edges in the gap panels; each turns
into a thick magenta line (σ = 1). Press Decode. The result shows the status, cost, solver and
time; the witness is drawn in teal (thick edges in layer panels, large dots in gap panels); and
"Syndrome verified" confirms that the viewer recomputed the witness's syndrome itself.

**3. Cost table (A3).** In Terminals mode (key 2), click two qubits: an edge in a layer panel or
a vertex in a gap panel. Numbered badges mark them, and the Cost table section lists them (use
the arrows to reorder). Press Compute table: four rows appear, and the smallest cost equals the
decode cost.

**4. Search (A4).** In the Spec section choose the preset wire. In Region mode (key 3), click
10 checks; dashed halos mark them, and the Search section reads "2^10 = 1,024 patterns:
exhaustive". Press Run search: the progress bar and counts update live. Press Cancel search to
see a run stop, then Run search again and let it finish (about 10 s on a 4-core laptop). Click a
row of the results list: its syndrome is loaded onto the canvas, labelled "From search result 1",
together with its cost table. Use Show on a table row to draw that row's witness.

Results are tied to their inputs. Edit the syndrome, terminals, scope, patch, region or spec and
the affected results are dimmed and labelled stale; Undo (Ctrl+Z) brings them back.

**5. Save and inspect offline (A5).** Press Save session; the browser downloads
`petersen_L3.session.json`. Open `viewer/index.html` straight from disk (for example
`file:///path/to/gadget-decoder/viewer/index.html`), press Open session and choose the file.
Everything is displayed; compute buttons are disabled in this static mode.

**6. Reproduce with the CLI (A6).** The CLI recomputes each stored result and compares it with
the saved one, ignoring timing fields:

```sh
gadgetlab decode petersen_L3.session.json --check
gadgetlab table  petersen_L3.session.json --check
gadgetlab search petersen_L3.session.json --check
```

Each prints "identical to the stored result" (exit code 0), or lists the differences (exit code
1). The same works for the sessions in `examples/`.

## Command line

```text
gadgetlab build (--preset NAME [--n N] [--seed S] | --edges EDGES) -L L [--out FILE]
gadgetlab decode    SESSION [--solver ilp_cbc|brute|ilp_highs] [--out FILE] [--check] [--json]
gadgetlab table     SESSION [--solver ...] [--out FILE] [--check] [--json]
gadgetlab search    SESSION [--budget B] [--seed S] [--top-k K] [--max-weight W] [--workers W]
                            [--solver ...] [--out FILE] [--check] [--json] [--quiet]
gadgetlab distances SESSION [--out FILE]
gadgetlab serve [--port 8765] [--r0-preset] [--debug | --no-debug]
```

- Presets: `k4`, `k33`, `prism`, `cube`, `petersen`, `heawood`, `random_cubic` (needs `--n`;
  if the graph is disconnected the seed is increased, and the seed used is recorded). `--edges`
  takes a string such as `"0-1 0-2 0-3 1-2 1-3 2-3"` or a file (JSON `[[u, w], ...]` or text).
  The engine accepts 2 ≤ L ≤ 8.
- `decode`, `table` and `search` write their result into the session file, or into `--out`.
  `--check` writes nothing. `search --check` reuses the parameters stored with the result.
- Defaults: budget 4096, seed 0, top-k 50, no max weight, workers = CPU count − 1, solver ilp_cbc.
- Output is deterministic for a fixed seed, apart from the timing fields `ms` and `elapsed_s`.
  The CLI and the HTTP API produce the same result JSON for the same session.
- Sessions can be written by hand; see §5 of the spec. Any invalid field is rejected with a
  message that names it, such as `terminals[1].v: 12 is out of range 0..9`.

## Viewer

| Key | Mode | Click |
|---|---|---|
| 1 | Syndrome | a check (gap-panel edge) toggles σ |
| 2 | Terminals | a qubit adds or removes a terminal (at most 6) |
| 3 | Region | a check adds or removes it from the region |
| 4 | Patch | a vertex in any panel adds or removes it from the patch vertex set U |

- Ctrl+Z and Ctrl+Shift+Z (or Ctrl+Y) undo and redo the last 50 edits. Shift-drag a vertex to
  move it; all panels share the layout, which is saved in the session.
- Scope global decodes the whole code. Scope local_closed allows errors only inside the patch,
  given by the vertex set U and a layer range; "Fill patch with ball" sets U to ball(v, r).
- Panels can sit in a grid or in one scrolling row.
- `index.html?selftest=1` runs the viewer's own support code on the T4 vectors and checks its
  hash code against the engine; it works offline.
- `gadgetlab serve --r0-preset`, or `?r0_preset=1` in static mode, adds a "Fill gap 0" button
  and a panel comparing MW(σ) with 3n/2 − 2α(G). Without the flag neither is shown.

## HTTP API

The server binds to 127.0.0.1 only. POST bodies are JSON (`Content-Type: application/json`).

| Method and path | Body | Returns |
|---|---|---|
| GET / | | the viewer |
| GET /api/config | | `{"r0_preset", "debug"}` |
| POST /api/build | `{"graph": {"preset": "petersen"}, "L": 3}` (or `"n"`, `"seed"`, `"edges"`) | `{"session", "info"}` |
| POST /api/decode | session | decode result |
| POST /api/table | session | table result |
| POST /api/distances | session | `{"dX", "dZ"}` |
| POST /api/random_syndrome | session plus `"weight"`, `"seed"` | `{"syndrome"}` |
| POST /api/search | session plus `"budget"`, `"seed"`, `"top_k"`, `"max_weight"`, `"workers"` | `{"job_id"}` |
| GET /api/jobs/ID | | `{"state", "progress", "summary", "result"}` |
| POST /api/jobs/ID/cancel | | `{"state"}` |

Decode, table and search also accept an optional `"solver"`. Invalid input answers HTTP 400 with
`{"error": message}`. With debug on (the default; `--no-debug` turns it off) every result is
checked by `verify()` and unexpected errors include a traceback.

## Code layout

```text
gadgetlab/
  f2.py        GF(2) rank, column-space test, solving, kernel basis
  graphs.py    presets, canonical edges, layouts, ball(), girth, alpha
  geometry.py  typed addresses <-> indices, supports, patches, panel data
  code.py      H1, R_L, hgp() and code info; checks the Kronecker matrices against the supports
  distance.py  dX and dZ by exhaustive search up to weight 4
  solve.py     constrained decode: ilp_cbc, brute, ilp_highs; memoisation; verify()
  table.py     cost tables and their invariants
  spec.py      specs, presets, matching and margin
  search.py    candidates, lazy evaluation, process pool, progress, cancel, jobs
  session.py   session schema, validation, load, save, inputs_hash
  compute.py   session-level computations shared by the CLI and the server
  server.py    HTTP API and viewer
  cli.py       command line
  errors.py    input errors (HTTP 400, CLI exit code 2)
viewer/index.html   the viewer
tests/              T1 to T12, see below
examples/           sessions with results: k4_L2 and petersen_L3, plus three used by tutorial.md
```

`solve`, `table`, `spec` and `search` work only with matrices and integer indices; `geometry`
translates typed addresses at the boundary.

## Tests

| Spec | File |
|---|---|
| T1 construction, T4 single-qubit syndromes | `tests/test_construction.py` |
| T2 code parameters (dX and dZ are listed in the test summary) | `tests/test_code_params.py` |
| T3 solver cross-check | `tests/test_solver_crosscheck.py` |
| T5 regression of the recorded claim | `tests/test_regression.py` |
| T6 table invariants, T11 infeasibility | `tests/test_tables.py` |
| T7 spec matching | `tests/test_spec.py` |
| T8 search | `tests/test_search.py` |
| T9 sessions | `tests/test_session.py` |
| T10 API smoke test | `tests/test_api.py` |
| T12 determinism, CLI | `tests/test_cli.py` |
| graph presets and layouts | `tests/test_graphs.py` |
| viewer self-test vectors and its hash code | `tests/test_viewer.py` |

T2 compares against `min_logical_weight` from `warmup_check.py`, which the tests load without
running its script part. After changing the engine's index map or hashing, refresh the vectors
embedded in the viewer with `python tests/selftest_vectors.py`.
