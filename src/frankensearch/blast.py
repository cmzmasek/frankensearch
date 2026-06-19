"""Thin helpers for locating and querying the BLAST+ command-line tools."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from . import scoring

REQUIRED_TOOLS = ("blastp", "makeblastdb", "blastdbcmd")

# A short sequence used for the matrix self-test (doctor).
_SELFTEST_SEQUENCE = "MQIFVKTLTGKTITLEVEPSDT"


def find_tool(name: str) -> str | None:
    """Return the full path to a BLAST+ tool, or None if it is not on PATH."""
    return shutil.which(name)


def tool_version(name: str) -> str | None:
    """Return the first line of ``<tool> -version``, or None if unavailable."""
    if find_tool(name) is None:
        return None
    try:
        result = subprocess.run(
            [name, "-version"], capture_output=True, text=True, timeout=15
        )
    except (OSError, subprocess.SubprocessError):
        return None
    lines = result.stdout.strip().splitlines()
    return lines[0].strip() if lines else None


def selftest_matrix(matrix: str = "identity") -> tuple[bool, str]:
    """Confirm blastp can load a matrix by running a tiny self-alignment.

    Returns ``(ok, detail)`` for display by ``doctor``.
    """
    if find_tool("blastp") is None or find_tool("makeblastdb") is None:
        return False, "BLAST+ tools not available"

    with tempfile.TemporaryDirectory() as raw_dir:
        work = Path(raw_dir)
        (work / "s.fa").write_text(f">s\n{_SELFTEST_SEQUENCE}\n")
        (work / "q.fa").write_text(f">q\n{_SELFTEST_SEQUENCE}\n")
        made = subprocess.run(
            ["makeblastdb", "-in", str(work / "s.fa"), "-dbtype", "prot", "-out", str(work / "db")],
            capture_output=True, text=True,
        )
        if made.returncode != 0:
            return False, "makeblastdb failed during self-test"

        cmd = [
            "blastp", "-query", str(work / "q.fa"), "-db", str(work / "db"),
            "-comp_based_stats", "0", "-evalue", "200000", "-word_size", "2",
            "-outfmt", "6 pident", *scoring.blast_args(matrix, ungapped=False),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            err = result.stderr.strip().splitlines()
            return False, err[0][:160] if err else "blastp failed"
        pident = result.stdout.strip().splitlines()
        detail = f"{scoring.matrix_blast_name(matrix)} loads"
        if pident:
            detail += f" (self-match {pident[0]}%)"
        return True, detail
