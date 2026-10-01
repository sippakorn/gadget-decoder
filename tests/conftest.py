import ast
from collections import defaultdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
RECORDS: dict[str, list[str]] = defaultdict(list)


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="also run tests marked slow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="slow; run pytest with --runslow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)


@pytest.fixture
def record():
    """record(section, line): a line printed in the terminal summary, e.g. dX and dZ per case (T2)."""

    def add(section: str, line: str) -> None:
        RECORDS[section].append(line)

    return add


def pytest_terminal_summary(terminalreporter):
    for section in sorted(RECORDS):
        terminalreporter.write_sep("-", section)
        for line in sorted(RECORDS[section]):
            terminalreporter.write_line(line)


@pytest.fixture(scope="session")
def warmup():
    """The functions of warmup_check.py, without running its script part."""
    path = ROOT / "warmup_check.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    keep = [node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
    namespace: dict = {}
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"), namespace)
    return namespace
