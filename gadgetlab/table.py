"""Cost tables (§3.10) and their runtime invariants.

Works only with a check matrix (through a Decoder) and integer qubit indices.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from gadgetlab.errors import InputError
from gadgetlab.solve import Decoder, Result
from gadgetlab.spec import bit_rows

MAX_TERMINALS = 6


class InvariantError(RuntimeError):
    """A cost table violates an invariant of §3.10. This is an error, never a warning."""


@dataclass(frozen=True)
class CostTable:
    terminals: tuple[int, ...]
    rows: dict[str, Result]  # bit string -> result, in lexicographic order
    unconstrained: Result

    def costs(self) -> dict[str, int | None]:
        return {a: r.cost for a, r in self.rows.items()}


def validate_terminals(decoder: Decoder, terminals) -> tuple[int, ...]:
    terminals = tuple(int(t) for t in terminals)
    if not 1 <= len(terminals) <= MAX_TERMINALS:
        raise InputError(f"a cost table needs between 1 and {MAX_TERMINALS} terminals, got {len(terminals)}")
    if len(set(terminals)) != len(terminals):
        raise InputError("the terminals must be distinct qubits")
    for k, t in enumerate(terminals, start=1):
        if not 0 <= t < decoder.N:
            raise InputError(f"terminal {k} (qubit {t}) is out of range 0..{decoder.N - 1}")
        if not decoder.in_Q[t]:
            raise InputError(f"terminal {k} (qubit {t}) lies outside the patch")
    return terminals


def assignment(terminals, a: str) -> tuple[tuple[int, int], ...]:
    """phi with phi(t_i) = a_i."""
    return tuple(sorted((t, int(bit)) for t, bit in zip(terminals, a)))


def compute_table(decoder: Decoder, sigma, terminals) -> CostTable:
    terminals = validate_terminals(decoder, terminals)
    rows = {a: decoder.solve(sigma, assignment(terminals, a)) for a in bit_rows(len(terminals))}
    table = CostTable(terminals, rows, decoder.solve(sigma, ()))
    check_invariants(table)
    return table


def check_invariants(table: CostTable) -> None:
    """(i) min over rows of Cost(a) equals the cost with no terminals fixed; (ii) the row matching
    the terminal values of the unconstrained witness has that same cost."""

    def cost(r: Result) -> float:
        return math.inf if r.status != "optimal" else r.cost

    row_min = min(cost(r) for r in table.rows.values())
    free = cost(table.unconstrained)
    if row_min != free:
        raise InvariantError(f"invariant (i) fails: the minimum over rows is {row_min}, the unconstrained cost is {free}")
    if table.unconstrained.status == "optimal":
        a = "".join(str(int(table.unconstrained.witness[t])) for t in table.terminals)
        if cost(table.rows[a]) != free:
            raise InvariantError(
                f"invariant (ii) fails: row {a} of the unconstrained witness has cost {cost(table.rows[a])}, not {free}")
