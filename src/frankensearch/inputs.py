"""Read and validate query sequences.

Accepts three formats, auto-detected:

* **FASTA**  -- ``>id description`` lines followed by sequence lines.
* **TSV**    -- two columns: name, amino-acid sequence (tab separated).
* **CSV**    -- two columns: name, amino-acid sequence (comma separated).

The goal is to be forgiving and to explain problems in plain language: bad
records are skipped with a warning rather than crashing, and only a genuinely
unusable file (missing, empty, or with no valid sequences at all) raises an
error.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from .errors import UserError

# Standard 20 amino acids plus the ambiguity / non-standard codes BLAST accepts.
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")
NUCLEOTIDE_ALPHABET = set("ACGTUN")

# Labels that signal a header row in a TSV/CSV table. A row is only treated as a
# header when its FIRST cell looks like a name column AND its SECOND cell looks
# like a sequence column — so a real first record (e.g. a protein named "query",
# or a 2-residue sequence "AA") is not mistaken for a header.
_HEADER_NAME_TOKENS = {"name", "id", "identifier", "accession", "protein", "query", "label", "gene"}
_HEADER_SEQ_TOKENS = {
    "seq", "sequence", "aa", "aaseq", "aa_seq", "aa_sequence",
    "peptide", "residues", "protein_sequence",
}

# Sequences longer than this are flagged as suspicious (likely not a short franken peptide).
_LONG_SEQUENCE_WARN = 10_000


class InputFormat(str, Enum):
    fasta = "fasta"
    tsv = "tsv"
    csv = "csv"


@dataclass(frozen=True)
class Query:
    """A single validated query sequence."""

    id: str
    sequence: str
    description: str = ""


@dataclass
class ParseResult:
    fmt: InputFormat
    records: list[Query]
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def parse(path: Path) -> ParseResult:
    """Parse and validate an input file, returning the valid records and messages."""
    text = _read_text(path)
    fmt, delimiter = _detect_format(text, path)

    if fmt is InputFormat.fasta:
        raw_records = _parse_fasta(text)
    else:
        raw_records = _parse_table(text, delimiter)

    result = ParseResult(fmt=fmt, records=[])
    seen_ids: set[str] = set()
    nucleotide_like = 0

    for rec_id, description, raw_seq in raw_records:
        rec_id = rec_id.strip()
        if not rec_id:
            result.warnings.append("Skipped a record with an empty name/ID.")
            continue

        sequence, problem = _clean_and_check(raw_seq)
        if problem:
            result.warnings.append(f"Skipped '{rec_id}': {problem}.")
            continue

        if rec_id in seen_ids:
            unique_id = _make_unique_id(rec_id, seen_ids)
            result.warnings.append(
                f"Duplicate query ID '{rec_id}' renamed to '{unique_id}' to keep results distinct."
            )
            rec_id = unique_id
        seen_ids.add(rec_id)

        if len(sequence) > _LONG_SEQUENCE_WARN:
            result.warnings.append(
                f"'{rec_id}' is {len(sequence):,} residues long — "
                "unusually long for a franken peptide."
            )
        if _looks_nucleotide(sequence):
            nucleotide_like += 1

        result.records.append(Query(id=rec_id, sequence=sequence, description=description.strip()))

    if nucleotide_like:
        result.warnings.append(
            f"{nucleotide_like} sequence(s) look like nucleotides (only A/C/G/T/U/N) — "
            "FRANKENSEARCH expects amino-acid sequences."
        )

    if not result.records:
        raise UserError(
            f"No valid sequences were found in {path}.",
            hint=(
                "Check that the file is FASTA, or a TSV/CSV with name in "
                "column 1 and sequence in column 2."
            ),
        )

    return result


# --------------------------------------------------------------------------- #
# Reading & format detection
# --------------------------------------------------------------------------- #
def _read_text(path: Path) -> str:
    if not path.exists():
        raise UserError(
            f"Input file not found: {path}",
            hint="Check the path (and spelling) and try again.",
        )
    if path.is_dir():
        raise UserError(
            f"Input path is a directory, not a file: {path}",
            hint="Point --input at a FASTA, TSV, or CSV file.",
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise UserError(f"Input file is empty: {path}")
    return text


def _detect_format(text: str, path: Path) -> tuple[InputFormat, str | None]:
    first_line = next((ln for ln in text.splitlines() if ln.strip()), "")

    if first_line.lstrip().startswith(">"):
        return InputFormat.fasta, None

    sample = "\n".join(text.splitlines()[:20])
    tabs, commas = sample.count("\t"), sample.count(",")
    if tabs and tabs >= commas:
        return InputFormat.tsv, "\t"
    if commas:
        return InputFormat.csv, ","

    raise UserError(
        f"Could not recognize the format of {path}.",
        hint=(
            "Expected FASTA (lines starting with '>') or a two-column TSV/CSV "
            "(name then amino-acid sequence)."
        ),
    )


# --------------------------------------------------------------------------- #
# Per-format parsing
# --------------------------------------------------------------------------- #
def _parse_fasta(text: str) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        if header is not None:
            ident, _, description = header.partition(" ")
            records.append((ident, description, "".join(seq_parts)))

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            seq_parts = []
        elif header is not None:
            seq_parts.append(line.strip())
    flush()
    return records


def _parse_table(text: str, delimiter: str | None) -> list[tuple[str, str, str]]:
    rows = [
        row
        for row in csv.reader(text.splitlines(), delimiter=delimiter or ",")
        if row and any(cell.strip() for cell in row)
    ]
    if rows and _looks_like_header(rows[0]):
        rows = rows[1:]

    records: list[tuple[str, str, str]] = []
    for row in rows:
        if len(row) < 2 or not row[1].strip():
            ident = row[0].strip() if row else "?"
            # A name with no sequence column; report against the name we have.
            records.append((ident, "", ""))
            continue
        records.append((row[0], "", row[1]))
    return records


def _looks_like_header(row: list[str]) -> bool:
    if len(row) < 2:
        return False
    name_cell, seq_cell = row[0].strip().lower(), row[1].strip().lower()
    return name_cell in _HEADER_NAME_TOKENS and seq_cell in _HEADER_SEQ_TOKENS


# --------------------------------------------------------------------------- #
# Validation helpers
# --------------------------------------------------------------------------- #
def _clean_and_check(raw_seq: str) -> tuple[str, str | None]:
    """Return (cleaned uppercase sequence, problem-description-or-None)."""
    sequence = "".join(raw_seq.split()).upper()
    if not sequence:
        return "", "no sequence in the second column"
    invalid = sorted({c for c in sequence if c not in AA_ALPHABET})
    if invalid:
        shown = " ".join(repr(c) for c in invalid[:8])
        return sequence, f"contains invalid amino-acid character(s): {shown}"
    return sequence, None


def _looks_nucleotide(sequence: str) -> bool:
    return len(sequence) >= 30 and set(sequence) <= NUCLEOTIDE_ALPHABET


def _make_unique_id(base: str, taken: set[str]) -> str:
    """Return base with a numeric suffix that is not already in ``taken``."""
    n = 2
    while f"{base}__{n}" in taken:
        n += 1
    return f"{base}__{n}"
