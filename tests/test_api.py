"""T10 API smoke test: start the server on a free port and call every endpoint; the API and the
CLI layer give identical result JSON for the same session (E9)."""

import json
import threading
import time
import urllib.error
import urllib.request

import pytest

from gadgetlab import compute
from gadgetlab.search import SearchParams
from gadgetlab.server import make_server
from gadgetlab.session import from_json


@pytest.fixture(scope="module")
def server():
    httpd = make_server(0, quiet=True)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.port}", httpd
    httpd.cancel_all()
    httpd.shutdown()
    httpd.server_close()


def call(base, method, path, body=None, headers=None, raw=None):
    data = raw if raw is not None else (None if body is None else json.dumps(body).encode())
    req = urllib.request.Request(base + path, data=data, method=method)
    if data is not None:
        req.add_header("Content-Type", "application/json")
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            ctype = resp.headers.get("Content-Type", "")
            payload = resp.read()
            return resp.status, (json.loads(payload) if "json" in ctype else payload.decode())
    except urllib.error.HTTPError as err:
        return err.code, json.loads(err.read())


def _wait(base, job_id, timeout=120):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        status, snap = call(base, "GET", f"/api/jobs/{job_id}")
        assert status == 200
        if snap["state"] != "running":
            return snap
        time.sleep(0.1)
    raise AssertionError("job did not finish")


def test_t10_every_endpoint_once(server):
    base, httpd = server
    status, page = call(base, "GET", "/")
    assert status == 200 and "<title>GadgetLab</title>" in page

    assert call(base, "GET", "/api/config") == (200, {"r0_preset": False, "debug": True})

    status, built = call(base, "POST", "/api/build", {"graph": {"preset": "k4"}, "L": 2})
    assert status == 200
    assert built["info"] == {"name": "K4", "preset": "k4", "n": 4, "nE": 6, "L": 2, "N": 16, "m": 6,
                             "x_checks": 8, "k": 3, "girth": 3}
    session = built["session"]
    session["syndrome"] = [{"e": 0, "t": 0}, {"e": 3, "t": 0}]
    session["terminals"] = [{"q": "E", "e": 1, "l": 0}, {"q": "V", "v": 2, "t": 0}]
    session["spec"] = {"preset": "wire", "rows": {"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1], "11": ["eq", 0]}}
    session["region"] = [{"e": e, "t": 0} for e in range(4)]
    local = from_json(session)

    status, decode = call(base, "POST", "/api/decode", session)
    assert status == 200 and decode["status"] == "optimal" and decode["verified"] is True
    assert compute.strip_timing(decode) == compute.strip_timing(compute.decode(local))

    status, table = call(base, "POST", "/api/table", session)
    assert status == 200 and [r["a"] for r in table["rows"]] == ["00", "01", "10", "11"]
    assert compute.strip_timing(table) == compute.strip_timing(compute.table(local))
    assert min(r["cost"] for r in table["rows"] if r["cost"] is not None) == decode["cost"]

    assert call(base, "POST", "/api/distances", session) == (200, {"dX": 2, "dZ": 3})

    status, rs = call(base, "POST", "/api/random_syndrome", {**session, "weight": 3, "seed": 7})
    assert status == 200 and rs == compute.random_syndrome(local, 3, 7)

    status, started = call(base, "POST", "/api/search",
                           {**session, "budget": 4096, "seed": 0, "top_k": 50, "max_weight": None, "workers": 1})
    assert status == 200
    snap = _wait(base, started["job_id"])
    assert snap["state"] == "done" and snap["progress"] == {"done": 16, "total": 16}
    cli = compute.search(local, SearchParams(workers=1))
    assert compute.strip_timing(snap["result"]) == compute.strip_timing(cli)
    assert snap["summary"] == snap["result"]["summary"]

    assert call(base, "POST", f"/api/jobs/{started['job_id']}/cancel") == (200, {"state": "done"})


def test_cancel_endpoint_stops_a_running_search(server):
    base, _ = server
    status, built = call(base, "POST", "/api/build", {"graph": {"preset": "petersen"}, "L": 3})
    session = built["session"]
    session["terminals"] = [{"q": "E", "e": 3, "l": 1}, {"q": "V", "v": 2, "t": 0}]
    session["spec"] = {"preset": "wire", "rows": {"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1], "11": ["eq", 0]}}
    session["region"] = [{"e": e, "t": t} for e in range(15) for t in (0, 1)]
    status, started = call(base, "POST", "/api/search", {**session, "budget": 4096, "workers": 2})
    assert status == 200
    time.sleep(1.0)
    status, state = call(base, "POST", f"/api/jobs/{started['job_id']}/cancel")
    assert (status, state) == (200, {"state": "cancelled"})
    assert call(base, "GET", f"/api/jobs/{started['job_id']}")[1]["state"] == "cancelled"


def test_errors(server):
    base, _ = server
    status, body = call(base, "POST", "/api/build", {"graph": {"edges": [[0, 1], [1, 2], [2, 0]]}, "L": 2})
    assert status == 400 and "even number of vertices" in body["error"]
    status, body = call(base, "POST", "/api/decode", {"schema": "gadgetlab.session/1"})
    assert status == 400 and "missing key 'graph'" in body["error"]
    _, built = call(base, "POST", "/api/build", {"graph": {"preset": "k4"}, "L": 2})
    bad = {**built["session"], "syndrome": [{"e": 9, "t": 0}]}
    status, body = call(base, "POST", "/api/decode", bad)
    assert status == 400 and body["error"].startswith("syndrome[0].e: 9 is out of range")
    status, body = call(base, "POST", "/api/table", built["session"])
    assert status == 400 and "terminal" in body["error"]
    status, body = call(base, "POST", "/api/search", {**built["session"], "budget": 0})
    assert status == 400 and "budget" in body["error"]
    assert call(base, "GET", "/api/jobs/nope")[0] == 404
    assert call(base, "GET", "/nothing")[0] == 404
    status, body = call(base, "POST", "/api/decode", raw=b"{", headers={})
    assert status == 400 and "not valid JSON" in body["error"]
    req = urllib.request.Request(base + "/api/decode", data=b"{}", method="POST", headers={"Content-Type": "text/plain"})
    with pytest.raises(urllib.error.HTTPError) as err:
        urllib.request.urlopen(req, timeout=10)
    assert err.value.code == 415
    status, body = call(base, "GET", "/api/config", headers={"Host": "evil.example:8765"})
    assert status == 403
    assert "traceback" not in body
