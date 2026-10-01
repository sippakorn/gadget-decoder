"""Error types shared by the engine, the CLI and the server."""


class InputError(ValueError):
    """Invalid user input. The server answers it with HTTP 400, the CLI with exit code 2."""
