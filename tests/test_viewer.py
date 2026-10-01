"""The viewer's self-test vectors (V13) match the engine; its own support and hash code agree with
the engine when run under Node (skipped without Node); the page loads nothing from the network."""

import json
import re
import shutil
import subprocess

import pytest

from selftest_vectors import VIEWER, build_vectors, embedded

HTML = VIEWER.read_text(encoding="utf-8")


def test_embedded_selftest_vectors_are_current():
    assert embedded() == build_vectors(), "run: python tests/selftest_vectors.py"


def test_viewer_loads_nothing_from_the_network():
    assert not re.search(r"<(link|img|iframe)\b", HTML)
    assert not re.search(r"<script[^>]*\bsrc=", HTML)
    assert "@import" not in HTML and "url(" not in HTML
    fetches = re.findall(r"""fetch\(([^)]*)\)""", HTML)
    assert fetches == ["path, opts"]
    assert set(re.findall(r'api\("([^"]+)"', HTML)) <= {"/api/config", "/api/build", "/api/decode", "/api/table",
                                                        "/api/distances", "/api/random_syndrome", "/api/search"}


CHECK = r"""
const Core = require(process.argv[2]);
const data = JSON.parse(process.argv[3]);
const geo = Core.makeGeometry(data.n, data.edges, data.L);
const t4 = data.t4.map(v => [...geo.syndromeOf([geo.qubitIndex(v.qubit)])].sort((a, b) => a - b));
const s = Core.validateSession(data.session);
const hashes = Object.fromEntries(Object.keys(data.hashes).map(k => [k, Core.inputsHash(k, s, k === "search" ? data.search_params : undefined)]));
console.log(JSON.stringify({t4, hashes, info: Core.codeInfo(s.graph, s.L)}));
"""


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is not installed")
def test_viewer_core_agrees_with_the_engine_under_node(tmp_path):
    core = re.search(r'<script id="core">(.*?)</script>', HTML, re.S).group(1)
    (tmp_path / "core.js").write_text(core, encoding="utf-8")
    (tmp_path / "check.js").write_text(CHECK, encoding="utf-8")
    data = build_vectors()
    out = subprocess.run(["node", str(tmp_path / "check.js"), str(tmp_path / "core.js"), json.dumps(data)],
                         capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    got = json.loads(out.stdout)
    from gadgetlab.session import from_json
    s = from_json(data["session"])
    geo = s.geometry
    assert got["t4"] == [sorted(geo.z_check(a["e"], a["t"]) for a in v["fires"]) for v in data["t4"]]
    assert got["hashes"] == data["hashes"]
    from gadgetlab.code import code_info
    engine = code_info(s.graph, s.L)
    assert {k: got["info"][k] for k in ("n", "nE", "L", "N", "m", "x_checks", "k", "girth")} == \
           {k: engine[k] for k in ("n", "nE", "L", "N", "m", "x_checks", "k", "girth")}
