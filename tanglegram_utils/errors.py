"""Custom exceptions: let the CLI print a clear message without a traceback."""


class ErrorDatos(Exception):
    """Invalid or inconsistent input (names that do not match, missing column, etc.)."""
