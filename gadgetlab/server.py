"""Local HTTP server: the viewer at / and the JSON API under /api (E10). Standard library only."""

from __future__ import annotations

import json
import sys
import threading
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

from gadgetlab import compute
from gadgetlab.errors import InputError
from gadgetlab.search import Evaluator, Job
from gadgetlab.session import from_json

HOST = "127.0.0.1"
DEFAULT_PORT = 8765
MAX_BODY = 64 * 1024 * 1024
MAX_JOBS = 32
VIEWER = Path(__file__).resolve().parent.parent / "viewer" / "index.html"


class GadgetLabServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, port: int = DEFAULT_PORT, r0_preset: bool = False, debug: bool = True,
                 viewer: Path | None = None, quiet: bool = False):
        super().__init__((HOST, port), Handler)
        self.r0_preset = r0_preset
        self.debug = debug
        self.viewer = Path(viewer) if viewer else VIEWER
        self.quiet = quiet
        self.jobs: dict[str, Job] = {}
        self.jobs_lock = threading.Lock()

    @property
    def port(self) -> int:
        return self.server_address[1]

    def add_job(self, job: Job) -> None:
        with self.jobs_lock:
            self.jobs[job.id] = job
            if len(self.jobs) > MAX_JOBS:
                for job_id, old in sorted(self.jobs.items(), key=lambda kv: kv[1].created):
                    if len(self.jobs) <= MAX_JOBS:
                        break
                    if old.state != "running":
                        del self.jobs[job_id]

    def cancel_all(self) -> None:
        with self.jobs_lock:
            for job in self.jobs.values():
                job.cancel()


class HTTPError(Exception):
    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status


def _split(body: dict, extras: tuple) -> tuple[object, dict]:
    """The session part of a request body and the extra top-level keys of this endpoint."""
    if not isinstance(body, dict):
        raise InputError("the request body must be a JSON object")
    body = dict(body)
    taken = {k: body.pop(k) for k in extras if k in body}
    return from_json(body), taken


class Handler(BaseHTTPRequestHandler):
    server: GadgetLabServer
    server_version = "GadgetLab/1.0"
    protocol_version = "HTTP/1.1"
    disable_nagle_algorithm = True  # headers and body go out in two writes; do not let the body wait for an ACK

    # ------------------------------------------------------------ plumbing

    def log_message(self, fmt, *args) -> None:
        if self.server.quiet or self.path.startswith("/api/jobs/"):
            return
        sys.stderr.write(f"{self.address_string()} {fmt % args}\n")

    def _send(self, status: int, payload: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(payload)

    def _json(self, status: int, obj) -> None:
        self._send(status, json.dumps(obj).encode("utf-8"), "application/json")

    def _check_host(self) -> None:
        host = (self.headers.get("Host") or "").rsplit(":", 1)[0].strip("[]")
        if host not in ("127.0.0.1", "localhost"):
            raise HTTPError(403, "requests must address 127.0.0.1 or localhost")

    def _read_raw(self) -> bytes | None:
        """Read the whole request body first, so every reply leaves the connection clean."""
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            length = -1
        if length < 0 or length > MAX_BODY:
            self.close_connection = True
            return None
        return self.rfile.read(length) if length else b""

    def _body(self):
        if self._raw is None:
            raise HTTPError(413, "the request body is too large or has an invalid Content-Length")
        ctype = (self.headers.get("Content-Type") or "").split(";")[0].strip().lower()
        if ctype != "application/json":
            raise HTTPError(415, "send the request body as application/json")
        try:
            return json.loads(self._raw.decode("utf-8") or "null")
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise HTTPError(400, f"the request body is not valid JSON: {exc}") from None

    def _handle(self, method: str) -> None:
        try:
            self._check_host()
            path = urlsplit(self.path).path
            if method == "GET":
                self._get(path)
            else:
                self._post(path)
        except HTTPError as exc:
            self._json(exc.status, {"error": str(exc)})
        except InputError as exc:
            self._json(400, {"error": str(exc)})
        except Exception as exc:  # engine failures are reported, never hidden
            payload = {"error": f"{type(exc).__name__}: {exc}"}
            if self.server.debug:
                payload["traceback"] = traceback.format_exc()
            self._json(500, payload)

    def do_GET(self) -> None:
        self._handle("GET")

    def do_POST(self) -> None:
        self._raw = self._read_raw()
        self._handle("POST")

    # ------------------------------------------------------------ routes

    def _get(self, path: str) -> None:
        if path in ("/", "/index.html"):
            try:
                page = self.server.viewer.read_bytes()
            except OSError:
                raise HTTPError(404, f"viewer not found at {self.server.viewer}") from None
            self._send(200, page, "text/html; charset=utf-8")
        elif path == "/api/config":
            self._json(200, {"r0_preset": self.server.r0_preset, "debug": self.server.debug})
        elif path.startswith("/api/jobs/"):
            self._json(200, self._job(path[len("/api/jobs/"):]).snapshot())
        else:
            raise HTTPError(404, f"no such path: {path}")

    def _job(self, job_id: str) -> Job:
        with self.server.jobs_lock:
            job = self.server.jobs.get(job_id)
        if job is None:
            raise HTTPError(404, f"unknown job {job_id!r}")
        return job

    def _post(self, path: str) -> None:
        if path.startswith("/api/jobs/") and path.endswith("/cancel"):
            job = self._job(path[len("/api/jobs/"):-len("/cancel")])
            self._json(200, {"state": job.cancel()})
            return
        body = self._body()
        debug = self.server.debug
        if path == "/api/build":
            if not isinstance(body, dict):
                raise InputError("the request body must be a JSON object")
            session, info = compute.build(body.get("graph"), body.get("L"))
            self._json(200, {"session": session.to_json(), "info": info})
        elif path == "/api/decode":
            session, extra = _split(body, ("solver",))
            self._json(200, compute.decode(session, extra.get("solver") or "ilp_cbc", debug))
        elif path == "/api/table":
            session, extra = _split(body, ("solver",))
            self._json(200, compute.table(session, extra.get("solver") or "ilp_cbc", debug))
        elif path == "/api/distances":
            session, _ = _split(body, ())
            self._json(200, compute.distances(session))
        elif path == "/api/random_syndrome":
            session, extra = _split(body, ("weight", "seed"))
            self._json(200, compute.random_syndrome(session, extra.get("weight"), extra.get("seed", 0)))
        elif path == "/api/search":
            session, extra = _split(body, ("budget", "seed", "top_k", "max_weight", "workers", "solver"))
            params = compute.search_params(extra).resolved()
            problem = compute.search_problem(session)
            Evaluator(problem, params.backend, debug)  # validate now, so bad input answers 400

            def target(progress, should_cancel):
                return compute.search(session, params, debug, progress, should_cancel)

            job = Job(target, debug=debug).start()
            self.server.add_job(job)
            self._json(200, {"job_id": job.id})
        else:
            raise HTTPError(404, f"no such path: {path}")


def make_server(port: int = DEFAULT_PORT, r0_preset: bool = False, debug: bool = True, quiet: bool = False,
                viewer: Path | None = None) -> GadgetLabServer:
    return GadgetLabServer(port, r0_preset=r0_preset, debug=debug, quiet=quiet, viewer=viewer)


def serve(port: int = DEFAULT_PORT, r0_preset: bool = False, debug: bool = True) -> None:
    httpd = make_server(port, r0_preset, debug)
    if not httpd.viewer.is_file():
        print(f"warning: viewer not found at {httpd.viewer}", file=sys.stderr)
    print(f"GadgetLab is serving http://{HOST}:{httpd.port}  (debug {'on' if debug else 'off'}, "
          f"r0 preset {'on' if r0_preset else 'off'}; Ctrl+C stops)", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping", file=sys.stderr)
    finally:
        httpd.cancel_all()
        httpd.server_close()
