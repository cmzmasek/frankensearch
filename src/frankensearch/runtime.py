"""Process-wide runtime state shared across the CLI."""

from dataclasses import dataclass


@dataclass
class RuntimeState:
    """Global flags set once from top-level CLI options."""

    debug: bool = False


STATE = RuntimeState()
