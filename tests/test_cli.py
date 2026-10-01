"""E9 CLI and T12 determinism: two runs with the same seed and one worker give identical output;
with more workers the matches and their ranking are identical."""

import json
import subprocess
import sys

import pytest

from gadgetlab import compute
from gadgetlab.cli import main


def run(*args, cwd=None):
    return subprocess.run([sys.executable, "-m", "gadgetlab", *args], cwd=cwd, capture_output=True, text=True, timeout=600)


def _prepare(path):
    obj = json.loads(path.read_text())
    obj["syndrome"] = [{"e": 0, "t": 0}, {"e": 5, "t": 0}]
    obj["terminals"] = [{"q": "E", "e": 1, "l": 1}, {"q": "V", "v": 3, "t": 0}]
    obj["spec"] = {"preset": "custom", "rows": {"00": ["eq", 0], "01": ["any"], "10": ["ge", 1], "11": ["ge", 0]}}
    obj["region"] = [{"e": e, "t": 0} for e in range(6)]
    path.write_text(json.dumps(obj, indent=2))


@pytest.fixture
def k4(tmp_path):
    path = tmp_path / "k4.json"
    result = run("build", "--preset", "k4", "-L", "2", "--out", str(path))
    assert result.returncode == 0, result.stderr
    assert "K4, L = 2: n 4, N 16, checks 6, k 3, girth 3" in result.stdout
    _prepare(path)
    return path


def test_t12_same_seed_one_worker_identical_output(k4, tmp_path):
    outs = []
    for name in ("a.json", "b.json"):
        out = tmp_path / name
        r = run("search", str(k4), "--budget", "40", "--seed", "3", "--workers", "1", "--quiet", "--out", str(out))
        assert r.returncode == 0, r.stderr
        outs.append(json.loads(out.read_text()))
    a, b = outs
    assert a["results"]["search"]["summary"]["exhaustive"] is False
    assert a["results"]["search"]["summary"]["candidates"] == 40
    assert compute.strip_timing(a) == compute.strip_timing(b)
    assert a["results"]["search"]["matches"], "the case should have matches"

    out = tmp_path / "c.json"
    r = run("search", str(k4), "--budget", "40", "--seed", "3", "--workers", "3", "--quiet", "--out", str(out))
    assert r.returncode == 0, r.stderr
    c = json.loads(out.read_text())["results"]["search"]
    assert compute.strip_timing(c["matches"]) == compute.strip_timing(a["results"]["search"]["matches"])
    assert c["params"]["workers"] == 3

    d = tmp_path / "d.json"
    r = run("search", str(k4), "--budget", "40", "--seed", "4", "--workers", "1", "--quiet", "--out", str(d))
    assert r.returncode == 0
    assert json.loads(d.read_text())["results"]["search"]["inputs_hash"] != a["results"]["search"]["inputs_hash"]


def test_t12_decode_and_table_runs_are_identical(k4, tmp_path):
    for command in ("decode", "table"):
        first, second = tmp_path / f"{command}1.json", tmp_path / f"{command}2.json"
        assert run(command, str(k4), "--out", str(first)).returncode == 0
        assert run(command, str(k4), "--out", str(second)).returncode == 0
        assert compute.strip_timing(json.loads(first.read_text())) == compute.strip_timing(json.loads(second.read_text()))


def test_results_are_written_back_and_check_reproduces_them(k4, capsys):
    assert main(["decode", str(k4)]) == 0
    assert main(["table", str(k4)]) == 0
    assert main(["search", str(k4), "--workers", "1", "--quiet"]) == 0
    stored = json.loads(k4.read_text())["results"]
    assert all(stored[k] is not None for k in ("decode", "table", "search"))
    capsys.readouterr()
    assert main(["decode", str(k4), "--check"]) == 0
    assert main(["table", str(k4), "--check"]) == 0
    assert main(["search", str(k4), "--check", "--quiet"]) == 0
    out = capsys.readouterr().out
    assert out.count("identical to the stored result") == 3

    obj = json.loads(k4.read_text())
    obj["results"]["decode"]["cost"] += 1
    k4.write_text(json.dumps(obj))
    assert main(["decode", str(k4), "--check"]) == 1
    assert "decode.cost" in capsys.readouterr().out


def test_build_from_edge_list_and_errors(tmp_path, capsys):
    out = tmp_path / "custom.json"
    assert main(["build", "--edges", "0-1 0-2 0-3 1-2 1-3 2-3", "-L", "3", "--out", str(out)]) == 0
    obj = json.loads(out.read_text())
    assert obj["graph"]["preset"] == "custom" and obj["graph"]["edges"][0] == [0, 1]
    assert main(["build", "--edges", "0-1 1-2 2-0", "-L", "3", "--out", str(out)]) == 2
    assert "even number of vertices" in capsys.readouterr().err
    assert main(["build", "--preset", "random_cubic", "--n", "10", "--seed", "2", "-L", "2", "--out", str(out)]) == 0
    assert json.loads(out.read_text())["graph"]["seed"] == 2
    assert main(["table", str(out)]) == 2
    assert "terminal" in capsys.readouterr().err


def test_distances_command(tmp_path, capsys):
    out = tmp_path / "p.json"
    assert main(["build", "--preset", "petersen", "-L", "2", "--out", str(out)]) == 0
    capsys.readouterr()
    assert main(["distances", str(out)]) == 0
    assert json.loads(capsys.readouterr().out) == {"dX": 2, "dZ": ">4"}
    assert main(["build", "--preset", "petersen", "-L", "3", "--out", str(out)]) == 0
    capsys.readouterr()
    assert main(["distances", str(out)]) == 0
    assert json.loads(capsys.readouterr().out) == {"dX": "not computed", "dZ": "not computed"}
