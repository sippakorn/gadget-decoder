"""Command-line entry points: gadgetlab build | decode | table | search | distances | serve (E9)."""

from __future__ import annotations

import argparse
import json
import re
import signal
import sys
import threading
import time
from pathlib import Path

from gadgetlab import compute
from gadgetlab.errors import InputError
from gadgetlab.graphs import PRESETS, from_edges, from_preset
from gadgetlab.search import SearchCancelled
from gadgetlab.session import dumps, load, new_session, save
from gadgetlab.solve import BACKENDS, DEFAULT_BACKEND


def main(argv=None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="gadgetlab", description="Syndrome-gadget explorer for HGP(cycle code of G, R_L).")
    sub = p.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build", help="create a session file from a preset or an edge list")
    src = b.add_mutually_exclusive_group(required=True)
    src.add_argument("--preset", choices=PRESETS, help="graph preset")
    src.add_argument("--edges", help='edge list: a file (JSON [[u, w], ...] or text) or a string such as "0-1 0-2 0-3 1-2 1-3 2-3"')
    b.add_argument("--n", type=int, help="number of vertices (random_cubic)")
    b.add_argument("--seed", type=int, default=0, help="seed for random_cubic (default 0)")
    b.add_argument("-L", "--L", dest="L", type=int, required=True, help="number of layers, 2 to 8")
    b.add_argument("--out", help="session file to write (default: <preset>_L<L>.session.json)")
    b.set_defaults(func=cmd_build)

    for name, func, text in (("decode", cmd_decode, "minimum-weight decode of the session syndrome"),
                             ("table", cmd_table, "cost table for the session terminals")):
        c = sub.add_parser(name, help=text)
        _session_args(c)
        c.add_argument("--solver", choices=BACKENDS, default=DEFAULT_BACKEND, help="solver backend (default ilp_cbc)")
        c.set_defaults(func=func)

    s = sub.add_parser("search", help="search syndrome patterns on the region for ones that meet the spec")
    _session_args(s)
    s.add_argument("--budget", type=int, help="candidate budget (default 4096)")
    s.add_argument("--seed", type=int, help="sampling seed (default 0)")
    s.add_argument("--top-k", dest="top_k", type=int, help="matches to keep (default 50)")
    s.add_argument("--max-weight", dest="max_weight", type=int, help="skip patterns with more ones than this")
    s.add_argument("--workers", type=int, help="worker processes (default: CPU count - 1)")
    s.add_argument("--solver", choices=BACKENDS, help="solver backend (default ilp_cbc)")
    s.add_argument("--quiet", action="store_true", help="no progress line")
    s.set_defaults(func=cmd_search)

    d = sub.add_parser("distances", help="dX and dZ by exhaustive search up to weight 4 (N <= 60)")
    d.add_argument("session", help="session file")
    d.add_argument("--out", help="write {dX, dZ} as JSON to this file")
    d.set_defaults(func=cmd_distances)

    v = sub.add_parser("serve", help="serve the viewer and the JSON API on http://127.0.0.1:PORT")
    v.add_argument("--port", type=int, default=8765, help="port (default 8765)")
    v.add_argument("--r0-preset", dest="r0_preset", action="store_true", help="enable the r0 preset (V12)")
    v.add_argument("--debug", action=argparse.BooleanOptionalAction, default=True,
                   help="verify every result and include tracebacks in error responses (default on)")
    v.set_defaults(func=cmd_serve)
    return p


def _session_args(c) -> None:
    c.add_argument("session", help="session file; results are written back into it unless --out is given")
    c.add_argument("--out", help="write the session with its new result here instead")
    c.add_argument("--check", action="store_true",
                   help="recompute, compare with the stored result (timing fields ignored) and write nothing")
    c.add_argument("--json", action="store_true", help="also print the result JSON")


# ---------------------------------------------------------------- build


def parse_edge_text(text: str) -> list[list[int]]:
    numbers = [int(t) for t in re.findall(r"\d+", text)]
    if len(numbers) % 2:
        raise InputError("the edge list has an odd number of vertex labels")
    return [numbers[k:k + 2] for k in range(0, len(numbers), 2)]


def cmd_build(args) -> int:
    if args.preset:
        graph = from_preset(args.preset, n=args.n, seed=args.seed)
    else:
        path = Path(args.edges)
        if path.is_file():
            text = path.read_text(encoding="utf-8")
            edges = json.loads(text) if path.suffix == ".json" else parse_edge_text(text)
        else:
            edges = parse_edge_text(args.edges)
        graph = from_edges(edges, args.n)
    session = new_session(graph, args.L)
    out = args.out or f"{graph.preset}_L{args.L}.session.json"
    save(session, out)
    info = compute.info(session)
    print(f"{info['name']}, L = {info['L']}: n {info['n']}, N {info['N']}, checks {info['m']}, "
          f"k {info['k']}, girth {info['girth']}")
    print(f"wrote {out}")
    return 0


# ---------------------------------------------------------------- decode, table


def _finish(args, session, kind: str, result: dict) -> int:
    if args.check:
        stored = session.results.get(kind)
        if stored is None:
            print(f"error: the session has no stored {kind} result to check against", file=sys.stderr)
            return 2
        diffs = compute.differences(stored, result, kind)
        if diffs:
            print(f"{kind}: the recomputed result differs from the stored one:")
            for line in diffs:
                print(f"  {line}")
            return 1
        print(f"{kind}: identical to the stored result (timing fields ignored)")
        return 0
    if args.json:
        print(dumps(result), end="")
    session.results[kind] = result
    out = args.out or args.session
    save(session, out)
    print(f"wrote {kind} result to {out}", file=sys.stderr if args.json else sys.stdout)
    return 0


def cmd_decode(args) -> int:
    session = load(args.session)
    result = compute.decode(session, args.solver)
    if not args.json:
        if result["status"] == "optimal":
            print(f"decode: optimal, cost {result['cost']} ({result['solver']}, {result['ms']} ms, verified)")
        else:
            print(f"decode: infeasible ({result['solver']}, {result['ms']} ms)")
    return _finish(args, session, "decode", result)


def cmd_table(args) -> int:
    session = load(args.session)
    result = compute.table(session, args.solver)
    if not args.json:
        for row in result["rows"]:
            print(f"  {row['a']}  {'—' if row['cost'] is None else row['cost']}")
        m = result["match"]
        if m is not None and m["c"] is None:
            print("  spec does not match: the eq rows are infeasible or give different values of c")
        elif m is not None:
            margin = "infinite" if m["margin"] is None else m["margin"]
            print(f"  spec {'matches' if m['matches'] else 'does not match'}: c = {m['c']}, margin = {margin}")
    return _finish(args, session, "table", result)


# ---------------------------------------------------------------- search


def cmd_search(args) -> int:
    session = load(args.session)
    given = {"budget": args.budget, "seed": args.seed, "top_k": args.top_k, "max_weight": args.max_weight,
             "workers": args.workers, "solver": args.solver}
    if args.check:
        stored = session.results.get("search")
        if stored is None:
            print("error: the session has no stored search result to check against", file=sys.stderr)
            return 2
        base = dict(stored["params"])
        base.update({k: v for k, v in given.items() if v is not None})
        params = compute.search_params(base)
    else:
        params = compute.search_params(given)

    cancel = threading.Event()
    previous = signal.signal(signal.SIGINT, lambda *_: cancel.set())
    reporter = _Progress(enabled=not args.quiet)
    try:
        result = compute.search(session, params, progress=reporter.update, should_cancel=cancel.is_set)
    except SearchCancelled as exc:
        reporter.close()
        s = exc.summary
        print(f"search cancelled after {s['rejected'] + s['matched'] + s['filtered'] + s['infeasible_syndrome']} "
              f"of {s['candidates']} candidates; nothing written", file=sys.stderr)
        return 130
    finally:
        signal.signal(signal.SIGINT, previous)
    reporter.close()
    s = result["summary"]
    if result["params"]["dropped_region"]:
        dropped = ", ".join(f"(e = {c['e']}, t = {c['t']})" for c in result["params"]["dropped_region"])
        print(f"warning: region checks outside the patch were dropped: {dropped}", file=sys.stderr)
    if not args.json:
        mode = "exhaustive" if s["exhaustive"] else "sampled"
        print(f"search ({mode}): {s['candidates']} candidates, {s['matched']} matched, {s['rejected']} rejected, "
              f"{s['filtered']} filtered, {s['infeasible_syndrome']} infeasible syndromes; "
              f"{s['solver_calls']} solver calls, {s['elapsed_s']} s")
        for m in result["matches"][:10]:
            margin = "infinite" if m["margin"] is None else m["margin"]
            print(f"  {m['rank']:>3}  pattern {m['pattern']}  c = {m['c']}  margin = {margin}")
    return _finish(args, session, "search", result)


class _Progress:
    def __init__(self, enabled: bool):
        self.enabled = enabled
        self.tty = sys.stderr.isatty()
        self.last = 0.0
        self.shown = False

    def update(self, done: int, total: int, summary: dict) -> None:
        if not self.enabled:
            return
        now = time.monotonic()
        if done < total and now - self.last < (0.2 if self.tty else 5.0):
            return
        self.last = now
        line = (f"search {done}/{total}  matched {summary['matched']}  rejected {summary['rejected']}  "
                f"solver calls {summary['solver_calls']}  {summary['elapsed_s']:.1f} s")
        if self.tty:
            print("\r" + line, end="", file=sys.stderr, flush=True)
            self.shown = True
        else:
            print(line, file=sys.stderr, flush=True)

    def close(self) -> None:
        if self.shown:
            print(file=sys.stderr)
            self.shown = False


# ---------------------------------------------------------------- distances, serve


def cmd_distances(args) -> int:
    session = load(args.session)
    result = compute.distances(session)
    text = json.dumps(result)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    return 0


def cmd_serve(args) -> int:
    from gadgetlab.server import serve

    serve(port=args.port, r0_preset=args.r0_preset, debug=args.debug)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
