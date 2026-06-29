#!/usr/bin/env python3

# Copyright (C) 2026 Gene S. Tan
# Copyright (C) 2026 Christian M. Zmasek
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import argparse
import csv
import os
import re
import sys

# Standard 20 amino acids plus the ambiguity / non-standard codes BLAST accepts.
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")

# Column labels that mark a header row in a TSV/CSV (col 1 = name, col 2 = seq).
_HEADER_NAME_TOKENS = {"name", "id", "identifier", "accession", "protein", "query", "label", "gene"}
_HEADER_SEQ_TOKENS = {
    "seq", "sequence", "aa", "aaseq", "aa_seq", "aa_sequence",
    "peptide", "residues", "protein_sequence",
}


def read_fasta(filepath):
    """
    Reads a FASTA file and returns a dict:
    {header: sequence}
    """
    sequences = {}
    header = None
    seq_chunks = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq_chunks).upper()
                header = line[1:].split()[0]  # take first token
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header:
            sequences[header] = "".join(seq_chunks).upper()

    return sequences


def _looks_like_header(row):
    """
    A row is treated as a header only if column 1 looks like a name label AND
    column 2 looks like a sequence label -- so a real first record (or a header
    word such as "sequence", which is itself all valid amino-acid letters) is
    never misclassified.
    """
    if len(row) < 2:
        return False
    return (row[0].strip().lower() in _HEADER_NAME_TOKENS
            and row[1].strip().lower() in _HEADER_SEQ_TOKENS)


def read_table(filepath, delimiter):
    """
    Reads a CSV/TSV file and returns a dict: {name: sequence}.
    Column 1 = name, column 2 = sequence; a header row and any further
    columns are ignored.
    """
    sequences = {}
    with open(filepath, newline="") as f:
        rows = [row for row in csv.reader(f, delimiter=delimiter)
                if row and any(cell.strip() for cell in row)]

    if rows and _looks_like_header(rows[0]):
        rows = rows[1:]

    for row in rows:
        if len(row) < 2:
            continue
        name_tokens = row[0].split()
        name = name_tokens[0] if name_tokens else ""  # first token (matches read_fasta)
        seq = "".join(row[1].split()).upper()  # drop whitespace, normalise case
        if not (name and seq):
            continue
        invalid = sorted({c for c in seq if c not in AA_ALPHABET})
        if invalid:
            shown = " ".join(repr(c) for c in invalid[:8])
            print(
                f"Warning: skipped '{name}': invalid amino-acid character(s): {shown}",
                file=sys.stderr,
            )
            continue
        sequences[name] = seq

    return sequences


def read_sequences(filepath):
    """
    Reads sequences from a FASTA, TSV, or CSV file (chosen by file extension)
    and returns a dict: {name: sequence}.
    """
    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".csv":
        return read_table(filepath, delimiter=",")
    if ext in (".tsv", ".tab"):
        return read_table(filepath, delimiter="\t")
    return read_fasta(filepath)


def extract_polypeptide_fragments(sequence, motif="GGSGGGGSGG", flanking_length=14):
    """
    Extract fragments centered on motif with flanking amino acids.
    Returns list of tuples: (fragment, start, end)
    """
    fragments = []

    for match in re.finditer(re.escape(motif), sequence):  # match the motif literally
        start = max(match.start() - flanking_length, 0)
        end = min(match.end() + flanking_length, len(sequence))
        fragments.append((sequence[start:end], start, end))

    return fragments


def process_file(input_path, output_path, motif, flanking_length):
    sequences = read_sequences(input_path)
    basename = os.path.basename(input_path)
    # Keep the input format in the output name so files that share a stem
    # (e.g. sars.fasta and sars.csv) don't overwrite each other.
    stem, ext = os.path.splitext(basename)
    suffix = ext.lstrip(".").lower()
    out_name = f"{stem}_{suffix}_extracted.fasta" if suffix else f"{stem}_extracted.fasta"
    out_file = os.path.join(output_path, out_name)

    with open(out_file, "w") as out:
        for seq_id, seq in sequences.items():
            fragments = extract_polypeptide_fragments(
                seq,
                motif=motif,
                flanking_length=flanking_length
            )

            for i, (fragment, start, end) in enumerate(fragments, 1):
                header = f">{basename}|{seq_id}|motif_{i}|{start}_{end}"
                out.write(header + "\n")
                out.write(fragment + "\n")

    return out_file


def main():
    parser = argparse.ArgumentParser(
        description="Extract motif-centered peptides from FASTA, TSV, or CSV files"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input folder of FASTA, TSV, or CSV files (TSV/CSV: col 1 = name, "
             "col 2 = sequence; header and extra columns ignored)",
    )
    parser.add_argument("-o", "--output", required=True, help="Output folder")
    parser.add_argument("-m", "--motif", default="GGSGGGGSGG", help="Motif sequence")
    parser.add_argument("-f", "--flank", type=int, default=14, help="Flanking length")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    input_files = [
        f for f in os.listdir(args.input)
        if f.lower().endswith((".fasta", ".fa", ".faa", ".tsv", ".tab", ".csv"))
    ]

    if not input_files:
        print("No FASTA/TSV/CSV files found in input directory.")
        return

    for name in input_files:
        input_path = os.path.join(args.input, name)
        output_file = process_file(
            input_path,
            args.output,
            args.motif,
            args.flank
        )
        print(f"Processed {name} -> {output_file}")


if __name__ == "__main__":
    main()
