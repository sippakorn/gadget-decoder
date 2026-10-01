"""Vectors embedded in viewer/index.html for its self-test (V13), computed by the engine.

Run `python tests/selftest_vectors.py` to rewrite them; tests/test_viewer.py checks they are current.
"""

import json
import re
from pathlib import Path

import numpy as np

from gadgetlab.code import build_code
from gadgetlab.graphs import from_preset
from gadgetlab.session import Patch, inputs_hash, new_session

VIEWER = Path(__file__).resolve().parent.parent / "viewer" / "index.html"
MARKERS = re.compile(r"/\*SELFTEST-BEGIN\*/(.*?)/\*SELFTEST-END\*/", re.S)


def build_vectors() -> dict:
    g = from_preset("k4")
    L = 3
    c = build_code(g.n, g.edges, L)
    geo = c.geometry
    t4 = []
    for j in range(geo.N):
        x = np.zeros(geo.N, dtype=np.int64)
        x[j] = 1
        fired = np.flatnonzero((c.HZ.astype(np.int64) @ x) % 2)
        t4.append({"qubit": geo.qubit_json(j), "fires": [geo.check_json(int(i)) for i in fired]})
    s = new_session(g, L)
    s.syndrome = sorted([geo.z_check(0, 0), geo.z_check(3, 1), geo.z_check(5, 0)])
    s.terminals = [geo.edge_qubit(1, 1), geo.vertex_qubit(2, 0)]
    s.region = [geo.z_check(2, 1), geo.z_check(0, 0), geo.z_check(4, 0)]
    s.spec = {"preset": "custom", "rows": {"00": ["eq", 0], "01": ["ge", 2], "10": ["any"], "11": ["eq", 1]}}
    s.scope = "local_closed"
    s.patch = Patch([0, 2], (0, 2))
    params = {"budget": 100, "seed": 7, "top_k": 5, "max_weight": None}
    return {
        "n": g.n, "edges": [list(e) for e in g.edges], "L": L, "t4": t4,
        "session": s.to_json(), "search_params": params,
        "hashes": {kind: inputs_hash(kind, s, params) for kind in ("decode", "table", "search")},
    }


def embedded() -> dict | None:
    match = MARKERS.search(VIEWER.read_text(encoding="utf-8"))
    return json.loads(match.group(1)) if match else None


def write() -> None:
    text = VIEWER.read_text(encoding="utf-8")
    payload = json.dumps(build_vectors(), separators=(",", ":"))
    new = MARKERS.sub(lambda _: f"/*SELFTEST-BEGIN*/{payload}/*SELFTEST-END*/", text, count=1)
    VIEWER.write_text(new, encoding="utf-8")


if __name__ == "__main__":
    write()
    print(f"wrote self-test vectors into {VIEWER}")
