"""Tests for input parsing, format detection, and validation."""

from pathlib import Path

import pytest

from frankensearch import inputs
from frankensearch.errors import UserError


def _write(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content)
    return path


def test_parse_fasta_multiline(tmp_path):
    path = _write(tmp_path, "q.fasta", ">a description here\nMKTAY\nIAKQR\n>b\nMQIFV\n")
    result = inputs.parse(path)
    assert result.fmt is inputs.InputFormat.fasta
    assert [q.id for q in result.records] == ["a", "b"]
    assert result.records[0].sequence == "MKTAYIAKQR"
    assert result.records[0].description == "description here"


def test_parse_tsv(tmp_path):
    path = _write(tmp_path, "q.tsv", "a\tMKTAY\nb\tMQIFV\n")
    result = inputs.parse(path)
    assert result.fmt is inputs.InputFormat.tsv
    assert [q.id for q in result.records] == ["a", "b"]


def test_parse_csv_with_header(tmp_path):
    path = _write(tmp_path, "q.csv", "name,sequence\na,MKTAY\nb,MQIFV\n")
    result = inputs.parse(path)
    assert result.fmt is inputs.InputFormat.csv
    assert [q.id for q in result.records] == ["a", "b"]


def test_first_row_named_like_a_header_word_is_kept(tmp_path):
    # 'query' as a NAME (with a real sequence) must not be dropped as a header.
    path = _write(tmp_path, "q.tsv", "query\tMKTAYIAKQR\nq2\tMQIFVKTLTG\n")
    assert [q.id for q in inputs.parse(path).records] == ["query", "q2"]


def test_short_sequence_first_row_is_kept(tmp_path):
    # A 2-residue sequence 'AA' in row 1 must not be mistaken for a header cell.
    path = _write(tmp_path, "q.csv", "x,AA\ny,MKTAY\n")
    assert [q.id for q in inputs.parse(path).records] == ["x", "y"]


def test_lowercase_is_normalised(tmp_path):
    path = _write(tmp_path, "q.fasta", ">a\nmktay\n")
    result = inputs.parse(path)
    assert result.records[0].sequence == "MKTAY"


def test_invalid_characters_skipped_with_warning(tmp_path):
    path = _write(tmp_path, "q.fasta", ">good\nMKTAY\n>bad\nMK1@Y\n")
    result = inputs.parse(path)
    ids = [q.id for q in result.records]
    assert ids == ["good"]
    assert any("bad" in w for w in result.warnings)


def test_duplicate_ids_are_renamed_to_stay_distinct(tmp_path):
    path = _write(tmp_path, "q.fasta", ">a\nMKTAY\n>a\nMQIFV\n>a\nMKTAY\n")
    result = inputs.parse(path)
    assert [q.id for q in result.records] == ["a", "a__2", "a__3"]
    assert any("renamed" in w.lower() for w in result.warnings)


def test_nucleotide_sequence_warning(tmp_path):
    path = _write(tmp_path, "q.fasta", f">a\n{'ACGT' * 20}\n")
    result = inputs.parse(path)
    assert any("nucleotide" in w.lower() for w in result.warnings)


def test_empty_file_raises(tmp_path):
    with pytest.raises(UserError):
        inputs.parse(_write(tmp_path, "q.fasta", "   \n"))


def test_missing_file_raises(tmp_path):
    with pytest.raises(UserError):
        inputs.parse(tmp_path / "does_not_exist.fasta")


def test_no_valid_sequences_raises(tmp_path):
    with pytest.raises(UserError):
        inputs.parse(_write(tmp_path, "q.fasta", ">a\n123456\n"))


def test_unrecognised_format_raises(tmp_path):
    with pytest.raises(UserError):
        inputs.parse(_write(tmp_path, "q.txt", "just one column of text\nmore text\n"))
