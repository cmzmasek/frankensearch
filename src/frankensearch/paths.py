"""Filesystem locations for FRANKENSEARCH data (caches, databases).

Everything lives under a single app home directory so it is easy to find and to
point ``--db-dir`` / diagnostics at. Override the root with the
``FRANKENSEARCH_HOME`` environment variable (handy for tests and CI).
"""

from __future__ import annotations

import os
from pathlib import Path

ENV_HOME = "FRANKENSEARCH_HOME"


def app_home() -> Path:
    override = os.environ.get(ENV_HOME)
    return Path(override).expanduser() if override else Path.home() / ".frankensearch"


def taxonomy_dir() -> Path:
    return app_home() / "taxonomy"


def database_dir() -> Path:
    return app_home() / "blastdb"
