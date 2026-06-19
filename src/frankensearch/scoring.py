"""Map FRANKENSEARCH scoring choices to valid ``blastp`` parameters.

BLAST+ only accepts its built-in matrices and a fixed set of gap costs per matrix;
a custom (e.g. +5/-4) matrix cannot be loaded from the command line in this build.
We therefore use BLAST's built-in **IDENTITY** matrix by default — pure identity
scoring, which is exactly what this non-homology search wants — and let the other
built-in matrices be selectable.

For IDENTITY the only gapped costs BLAST+ allows are 15/2 (or ungapped), so we set
those explicitly. For the other matrices we omit gap costs and let blastp use its
own (valid) defaults.
"""

from __future__ import annotations

MATRICES = ("identity", "pam30", "blosum45", "blosum62")

# Explicit gap costs where we need a specific supported setting; otherwise blastp's
# per-matrix default is used.
_GAP_COSTS: dict[str, tuple[int, int]] = {
    "identity": (15, 2),
}


def matrix_blast_name(matrix: str) -> str:
    return matrix.upper()


def gap_costs(matrix: str) -> tuple[int, int] | None:
    return _GAP_COSTS.get(matrix)


def gap_description(matrix: str, *, ungapped: bool) -> str:
    if ungapped:
        return "ungapped"
    costs = gap_costs(matrix)
    return f"open {costs[0]}, extend {costs[1]}" if costs else "matrix default"


def blast_args(matrix: str, *, ungapped: bool) -> list[str]:
    """The ``-matrix`` / gap-related arguments for a local blastp run."""
    args = ["-matrix", matrix_blast_name(matrix)]
    if ungapped:
        args.append("-ungapped")
    else:
        costs = gap_costs(matrix)
        if costs:
            args += ["-gapopen", str(costs[0]), "-gapextend", str(costs[1])]
    return args


def remote_blast_args(matrix: str, *, ungapped: bool) -> tuple[list[str], str | None]:
    """``-matrix`` / gap arguments for a remote NCBI run, plus an optional warning.

    NCBI's remote service doesn't offer the IDENTITY matrix, so we fall back to
    PAM30 (the closest identity-favouring built-in) and report why. Gap costs are
    left to NCBI's per-matrix defaults.
    """
    warning: str | None = None
    name = matrix_blast_name(matrix)
    if matrix == "identity":
        name = "PAM30"
        warning = (
            "Remote NCBI BLAST does not support the IDENTITY matrix; "
            "falling back to PAM30 for this run."
        )
    args = ["-matrix", name]
    if ungapped:
        args.append("-ungapped")
    return args, warning
