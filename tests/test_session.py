"""T9 sessions: save then load gives an equal session; invalid input is rejected naming the field (E8)."""

import copy
import json

import pytest

from gadgetlab import compute
from gadgetlab.graphs import from_preset
from gadgetlab.session import SCHEMA, Session, SessionError, from_json, inputs, inputs_hash, load, new_session, save


def _session() -> Session:
    s = new_session(from_preset("petersen"), 3)
    geo = s.geometry
    s.syndrome = sorted([geo.z_check(0, 0), geo.z_check(4, 1), geo.z_check(9, 0)])
    s.terminals = [geo.edge_qubit(3, 1), geo.vertex_qubit(2, 0)]
    s.region = [geo.z_check(5, 1), geo.z_check(1, 0)]
    s.spec = {"preset": "wire", "rows": {"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1], "11": ["eq", 0]}}
    s.notes = "round trip"
    return s


def test_save_then_load_gives_an_equal_session(tmp_path):
    s = _session()
    s.results["decode"] = compute.decode(s)
    s.results["table"] = compute.table(s)
    path = tmp_path / "s.json"
    save(s, path)
    loaded = load(path)
    assert loaded == s
    assert loaded.to_json() == s.to_json()
    assert json.loads(path.read_text())["schema"] == SCHEMA


def test_section_5_example_shape_loads():
    obj = _session().to_json()
    obj["scope"] = "local_closed"
    obj["patch"] = {"vertices": [0, 1, 4], "layers": [0, 2]}
    s = from_json(obj)
    assert s.patch.vertices == [0, 1, 4] and s.patch.layers == (0, 2)
    assert obj["graph"]["edges"][:3] == [[0, 1], [0, 4], [0, 5]]


def _bad(mutate):
    obj = _session().to_json()
    mutate(obj)
    return obj


@pytest.mark.parametrize("mutate,field", [
    (lambda o: o["syndrome"].append({"e": 15, "t": 0}), r"syndrome\[3\]\.e: 15 is out of range 0\.\.14"),
    (lambda o: o["syndrome"].append({"e": 1, "t": 2}), r"syndrome\[3\]\.t: 2 is out of range 0\.\.1"),
    (lambda o: o["syndrome"].append(dict(o["syndrome"][0])), r"syndrome\[3\]: check .* appears twice"),
    (lambda o: o["terminals"].append({"q": "V", "v": 10, "t": 0}), r"terminals\[2\]\.v: 10 is out of range"),
    (lambda o: o["terminals"].__setitem__(0, {"q": "E", "e": 3, "l": 3}), r"terminals\[0\]\.l: 3 is out of range"),
    (lambda o: o["terminals"].__setitem__(1, {"q": "Z", "v": 1, "t": 0}), r"terminals\[1\]\.q"),
    (lambda o: o["region"].append({"e": -1, "t": 0}), r"region\[2\]\.e: -1 is out of range"),
    (lambda o: o.__setitem__("patch", {"vertices": [0, 12], "layers": [0, 1]}), r"patch\.vertices\[1\]: 12 is out of range"),
    (lambda o: o.__setitem__("patch", {"vertices": [0], "layers": [2, 1]}), r"patch\.layers: need l0 <= l1"),
    (lambda o: o.__setitem__("patch", {"vertices": [0], "layers": [0, 3]}), r"patch\.layers\[1\]: 3 is out of range"),
    (lambda o: o.__setitem__("scope", "local"), r"scope"),
    (lambda o: o.__setitem__("L", 9), r"L: 9 is out of range 2\.\.8"),
    (lambda o: o["spec"]["rows"].pop("11"), r"spec\.rows: missing row '11'"),
    (lambda o: o["spec"]["rows"].__setitem__("01", ["gt", 1]), r"spec\.rows\.01"),
    (lambda o: o["graph"]["edges"].reverse(), r"graph\.edges: edges must be in canonical order"),
    (lambda o: o["graph"]["layout"].__setitem__("3", [1.5, 0.2]), r"graph\.layout\.3"),
    (lambda o: o.__setitem__("schema", "gadgetlab.session/0"), r"schema"),
    (lambda o: o.__setitem__("extra", 1), r"unexpected key 'extra'"),
    (lambda o: o["results"].__setitem__("decode", {"status": "optimal", "cost": 1, "error": [{"q": "E", "e": 99, "l": 0}],
                                                   "verified": True, "solver": "x", "ms": 1, "inputs_hash": "h"}),
     r"results\.decode\.error\[0\]\.e: 99 is out of range"),
])
def test_invalid_input_is_rejected_naming_the_field(mutate, field):
    with pytest.raises(SessionError, match=field):
        from_json(_bad(mutate))


def test_terminal_count_and_spec_must_agree():
    obj = _session().to_json()
    obj["terminals"].append({"q": "E", "e": 0, "l": 0})
    with pytest.raises(SessionError, match="3 terminals"):
        from_json(obj)


def test_inputs_hash_tracks_what_each_result_depends_on():
    s = _session()
    base = {k: inputs_hash(k, s, {"budget": 10, "seed": 0, "top_k": 5, "max_weight": None}) for k in ("decode", "table", "search")}
    t = copy.deepcopy(s)
    t.region = list(reversed(t.region))  # region order matters for patterns
    assert inputs_hash("decode", t) == base["decode"] and inputs_hash("table", t) == base["table"]
    assert inputs_hash("search", t, {"budget": 10, "seed": 0, "top_k": 5}) != base["search"]
    u = copy.deepcopy(s)
    u.syndrome = sorted(set(u.syndrome) ^ {u.region[0]})  # an edit on a region check
    assert inputs_hash("decode", u) != base["decode"]
    assert inputs_hash("search", u, {"budget": 10, "seed": 0, "top_k": 5}) == base["search"]
    v = copy.deepcopy(s)
    v.spec = {"preset": "custom", "rows": {**s.spec["rows"], "01": ["ge", 2]}}
    assert inputs_hash("decode", v) == base["decode"] and inputs_hash("table", v) != base["table"]
    w = copy.deepcopy(s)
    w.patch = None
    w.notes = "changed"
    w.graph.layout = {k: (0.5, 0.5) for k in w.graph.layout}
    assert all(inputs_hash(k, w, {"budget": 10, "seed": 0, "top_k": 5}) == base[k] for k in base)
    assert inputs("decode", s)["patch"] is None  # the patch is irrelevant in global scope
