"""Tests for building/tracking BLAST databases.

These run the real ``makeblastdb`` (available in the dev env) on a tiny FASTA, but
stub the UniProt download so no network is used.
"""

from pathlib import Path

import pytest

from frankensearch import database, uniprot
from frankensearch.errors import UserError
from frankensearch.taxonomy import Taxon

TINY_FASTA = (
    ">sp|P1|A_TEST test protein one OS=Test\nMKTAYIAKQRQISFVKSHFSRQLEERLGL\n"
    ">sp|P2|B_TEST test protein two OS=Test\nMQIFVKTLTGKTITLEVEPSDTIENVKAK\n"
)

TAXON = Taxon(9606, "Homo sapiens", "species")


@pytest.fixture
def stub_download(monkeypatch):
    monkeypatch.setattr(
        uniprot, "build_query", lambda taxid, ps, *, client: f"(organism_id:{taxid})"
    )

    def fake_download(query, dest, *, client):
        Path(dest).write_text(TINY_FASTA)
        return uniprot.DownloadResult(2, release="2026_02", release_date="10-June-2026")

    monkeypatch.setattr(uniprot, "download_fasta", fake_download)


def test_build_creates_db_and_metadata(tmp_path, stub_download):
    meta = database.build(TAXON, "all", db_dir=tmp_path, client=None)
    assert meta.sequence_count == 2
    assert meta.scientific_name == "Homo sapiens"
    assert meta.uniprot_release == "2026_02"
    assert database.is_built(9606, tmp_path)
    assert database.load_metadata(9606, tmp_path) == meta
    assert [m.taxid for m in database.list_built(tmp_path)] == [9606]


def test_build_skips_when_already_present(tmp_path, monkeypatch):
    calls: list[int] = []
    monkeypatch.setattr(uniprot, "build_query", lambda taxid, ps, *, client: "q")

    def fake_download(query, dest, *, client):
        calls.append(1)
        Path(dest).write_text(TINY_FASTA)
        return uniprot.DownloadResult(2)

    monkeypatch.setattr(uniprot, "download_fasta", fake_download)

    database.build(TAXON, "all", db_dir=tmp_path, client=None)
    database.build(TAXON, "all", db_dir=tmp_path, client=None)  # should reuse, not re-download
    assert calls == [1]


def test_build_rebuilds_when_proteome_set_changes(tmp_path, monkeypatch):
    calls: list[str] = []
    monkeypatch.setattr(uniprot, "build_query", lambda taxid, ps, *, client: "q")

    def fake_download(query, dest, *, client):
        calls.append("dl")
        Path(dest).write_text(TINY_FASTA)
        return uniprot.DownloadResult(2)

    monkeypatch.setattr(uniprot, "download_fasta", fake_download)

    database.build(TAXON, "swissprot", db_dir=tmp_path, client=None)
    database.build(TAXON, "reference", db_dir=tmp_path, client=None)  # different set → rebuild
    assert calls == ["dl", "dl"]
    assert database.load_metadata(9606, tmp_path).proteome_set == "reference"


def test_is_current_matches_only_same_set(tmp_path, stub_download):
    database.build(TAXON, "all", db_dir=tmp_path, client=None)
    assert database.is_current(9606, "all", tmp_path)
    assert not database.is_current(9606, "swissprot", tmp_path)
    assert not database.is_current(10090, "all", tmp_path)


def test_build_zero_sequences_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(uniprot, "build_query", lambda taxid, ps, *, client: "q")

    def empty_download(query, dest, *, client):
        Path(dest).write_text("")
        return uniprot.DownloadResult(0)

    monkeypatch.setattr(uniprot, "download_fasta", empty_download)
    with pytest.raises(UserError):
        database.build(TAXON, "all", db_dir=tmp_path, client=None)
