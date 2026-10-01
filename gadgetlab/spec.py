"""Specs for cost tables: model, presets, matching and margin (§3.11)."""

from __future__ import annotations

from dataclasses import dataclass

from gadgetlab.errors import InputError

PRESET_NAMES = ("wire", "copy3", "custom")
PRESET_SIZES = {"wire": 2, "copy3": 3}


class SpecError(InputError):
    pass


def bit_rows(k: int) -> list[str]:
    """All bit strings of length k, lexicographic (all zeros first)."""
    return [format(i, f"0{k}b") for i in range(1 << k)] if k > 0 else []


def preset_rows(name: str, k: int) -> dict[str, list]:
    if name == "custom":
        return {a: ["any"] for a in bit_rows(k)}
    if name not in PRESET_SIZES:
        raise SpecError(f"unknown spec preset {name!r}; expected one of {', '.join(PRESET_NAMES)}")
    if k != PRESET_SIZES[name]:
        raise SpecError(f"the {name} preset needs {PRESET_SIZES[name]} terminals, not {k}")
    if name == "wire":
        return {"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1], "11": ["eq", 0]}
    return {a: (["eq", 0] if a in ("000", "111") else ["ge", 1]) for a in bit_rows(3)}


def parse_entry(entry, field: str) -> list:
    if not isinstance(entry, (list, tuple)) or not entry:
        raise SpecError(f'{field}: expected ["eq", k], ["ge", k] or ["any"], got {entry!r}')
    op = entry[0]
    if op == "any":
        if len(entry) != 1:
            raise SpecError(f'{field}: "any" takes no offset')
        return ["any"]
    if op not in ("eq", "ge"):
        raise SpecError(f'{field}: unknown entry {op!r}; expected "eq", "ge" or "any"')
    if len(entry) != 2:
        raise SpecError(f'{field}: expected ["{op}", k]')
    k = entry[1]
    if not isinstance(k, int) or isinstance(k, bool) or k < 0:
        raise SpecError(f"{field}: the offset must be a non-negative integer, got {k!r}")
    return [op, k]


def parse_rows(rows, k: int, field: str = "spec.rows") -> dict[str, list]:
    """Validate a rows object for k terminals; returns it with entries normalised, in bit order."""
    if not isinstance(rows, dict):
        raise SpecError(f"{field}: expected an object mapping bit strings to entries")
    expected = bit_rows(k)
    for a in rows:
        if a not in expected:
            raise SpecError(f"{field}: unexpected row {a!r} for {k} terminal{'s' if k != 1 else ''}")
    for a in expected:
        if a not in rows:
            raise SpecError(f"{field}: missing row {a!r}")
    return {a: parse_entry(rows[a], f"{field}.{a}") for a in expected}


def problems(rows: dict[str, list]) -> list[str]:
    """Reasons a well-formed spec cannot be matched; empty when it can."""
    if not rows:
        return ["Choose at least one terminal"]
    if not any(entry[0] == "eq" for entry in rows.values()):
        return ["Add at least one eq row"]
    return []


@dataclass(frozen=True)
class Match:
    matches: bool
    c: int | None
    margin: int | None  # None means infinity

    def to_json(self) -> dict:
        return {"matches": self.matches, "c": self.c, "margin": self.margin}


def match(rows: dict[str, list], costs: dict[str, int | None]) -> Match:
    """Match a table (costs, None = infinity) against a spec (§3.11).

    c is read off the eq rows; when they are infeasible or disagree, c and margin are None and
    the table does not match. Otherwise margin is the minimum of Cost(a) - c - k over ge rows
    with finite cost (None if there is none), and the table matches when no ge row falls short.
    """
    if problems(rows):
        raise SpecError(problems(rows)[0])
    c = None
    for a, entry in rows.items():
        if entry[0] != "eq":
            continue
        cost = costs[a]
        if cost is None:
            return Match(False, None, None)
        if c is None:
            c = cost - entry[1]
        elif cost - entry[1] != c:
            return Match(False, None, None)
    ok = True
    margin = None
    for a, entry in rows.items():
        if entry[0] != "ge" or costs[a] is None:
            continue
        slack = costs[a] - c - entry[1]
        ok = ok and slack >= 0
        margin = slack if margin is None else min(margin, slack)
    return Match(ok, c, margin)
