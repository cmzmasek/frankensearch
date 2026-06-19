"""User-facing error types.

Anything that is the *user's* problem (bad input, a missing dependency, a missing
database, ...) should be raised as a `FrankensearchError` so the CLI can print a
clean, actionable message instead of a Python traceback. Tracebacks are reserved
for genuine internal bugs and are only shown when the user passes ``--debug``.
"""

from __future__ import annotations


class FrankensearchError(Exception):
    """Base class for expected, user-facing errors.

    ``hint`` is an optional follow-up line suggesting how to fix the problem.
    """

    def __init__(self, message: str, hint: str | None = None) -> None:
        super().__init__(message)
        self.hint = hint


class UserError(FrankensearchError):
    """The user supplied something invalid (arguments, input files, ...)."""


class DependencyError(FrankensearchError):
    """A required external tool (e.g. BLAST+) is missing or not working."""
