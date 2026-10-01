"""T7 spec matching on hand-made tables, presets and parsing."""

import pytest

from gadgetlab.spec import SpecError, bit_rows, match, parse_rows, preset_rows, problems

WIRE = preset_rows("wire", 2)


def test_presets():
    assert WIRE == {"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1], "11": ["eq", 0]}
    copy3 = preset_rows("copy3", 3)
    assert copy3["000"] == copy3["111"] == ["eq", 0]
    assert all(copy3[a] == ["ge", 1] for a in bit_rows(3) if a not in ("000", "111"))
    assert preset_rows("custom", 2) == {a: ["any"] for a in ("00", "01", "10", "11")}
    with pytest.raises(SpecError):
        preset_rows("wire", 3)


def test_wire_matches_with_margin():
    m = match(WIRE, {"00": 4, "01": 6, "10": 5, "11": 4})
    assert (m.matches, m.c, m.margin) == (True, 4, 0)


def test_null_costs_in_ge_rows_satisfy_and_do_not_count_toward_margin():
    m = match(WIRE, {"00": 4, "01": None, "10": 7, "11": 4})
    assert (m.matches, m.c, m.margin) == (True, 4, 2)


def test_no_finite_ge_row_gives_infinite_margin():
    m = match(WIRE, {"00": 2, "01": None, "10": None, "11": 2})
    assert (m.matches, m.c, m.margin) == (True, 2, None)
    only_eq = {"0": ["eq", 0], "1": ["eq", 1]}
    m = match(only_eq, {"0": 3, "1": 4})
    assert (m.matches, m.c, m.margin) == (True, 3, None)


def test_inconsistent_eq_rows_do_not_match():
    m = match(WIRE, {"00": 3, "01": 9, "10": 9, "11": 7})
    assert (m.matches, m.c, m.margin) == (False, None, None)


def test_infeasible_eq_row_does_not_match():
    m = match(WIRE, {"00": None, "01": 5, "10": 5, "11": 3})
    assert (m.matches, m.c, m.margin) == (False, None, None)


def test_violated_ge_row_gives_negative_margin():
    m = match(WIRE, {"00": 4, "01": 4, "10": 6, "11": 4})
    assert (m.matches, m.c, m.margin) == (False, 4, -1)


def test_offsets_shift_c():
    rows = {"00": ["eq", 2], "01": ["ge", 3], "10": ["any"], "11": ["eq", 0]}
    m = match(rows, {"00": 7, "01": 9, "10": 0, "11": 5})
    assert (m.matches, m.c, m.margin) == (True, 5, 1)


def test_spec_without_eq_rows_is_rejected():
    rows = preset_rows("custom", 2)
    assert problems(rows) == ["Add at least one eq row"]
    with pytest.raises(SpecError, match="eq row"):
        match(rows, {a: 1 for a in rows})


@pytest.mark.parametrize("rows,message", [
    ({"00": ["eq", 0], "01": ["ge", 1], "10": ["ge", 1]}, "missing row '11'"),
    ({**WIRE, "111": ["any"]}, "unexpected row"),
    ({**WIRE, "01": ["ge", -1]}, "non-negative"),
    ({**WIRE, "01": ["gt", 1]}, "unknown entry"),
    ({**WIRE, "01": ["any", 1]}, "takes no offset"),
    ({**WIRE, "01": ["eq"]}, r'\["eq", k\]'),
])
def test_parse_rows_rejects_bad_specs(rows, message):
    with pytest.raises(SpecError, match=message):
        parse_rows(rows, 2)
