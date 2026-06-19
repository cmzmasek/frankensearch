"""Smoke tests for the CLI scaffolding."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from frankensearch import __version__, database, uniprot
from frankensearch.cli import _check_no_overwrite, app
from frankensearch.errors import UserError
from frankensearch.taxonomy import Taxon

runner = CliRunner()


def test_help_lists_commands():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    for command in ("search", "setup", "doctor"):
        assert command in result.stdout


def test_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.stdout


def test_search_dry_run(tmp_path, monkeypatch):
    # Avoid any network: stub taxonomy resolution.
    monkeypatch.setattr("frankensearch.cli._resolve_species", lambda taxid_list, **kwargs: [])
    query_file = tmp_path / "queries.fasta"
    query_file.write_text(">a\nMKTAYIAKQR\n>b\nMQIFVKTLTG\n")
    result = runner.invoke(app, ["search", str(query_file), "--taxids", "9606,10029", "--dry-run"])
    assert result.exit_code == 0
    assert "Search plan" in result.stdout
    assert "9606" in result.stdout


def test_search_rejects_bad_taxid():
    # Bad taxids are rejected before the input file is even read.
    result = runner.invoke(app, ["search", "queries.fasta", "--taxids", "not_a_number"])
    assert result.exit_code != 0


def test_doctor_runs(tmp_path, monkeypatch):
    # Point the app home at a temp dir so we never touch the real cache.
    monkeypatch.setenv("FRANKENSEARCH_HOME", str(tmp_path))
    result = runner.invoke(app, ["doctor"])
    # Exit code depends on whether BLAST+ is on PATH; either way it should render.
    assert "frankensearch doctor" in result.stdout


def test_check_no_overwrite(tmp_path):
    existing = tmp_path / "a.tsv"
    existing.write_text("x")
    with pytest.raises(UserError):
        _check_no_overwrite((existing,), force=False)
    _check_no_overwrite((existing,), force=True)  # --force allows it
    _check_no_overwrite((tmp_path / "missing.tsv",), force=False)  # absent is fine


def test_resolve_species_prefers_db_metadata(tmp_path, monkeypatch):
    from frankensearch import cli, taxonomy

    monkeypatch.setattr(uniprot, "build_query", lambda taxid, ps, *, client: "q")

    def fake_download(query, dest, *, client):
        Path(dest).write_text(">sp|P1|A_TEST a OS=Testus\nMKTAYIAKQR\n")
        return uniprot.DownloadResult(1)

    monkeypatch.setattr(uniprot, "download_fasta", fake_download)
    taxon = Taxon(9606, "Homo sapiens", "species")
    database.build(taxon, "swissprot", db_dir=tmp_path, client=None)

    def boom(self, taxid):
        raise AssertionError("UniProt resolve must not be called when DB metadata exists")

    monkeypatch.setattr(taxonomy.TaxonomyResolver, "resolve", boom)
    taxa = cli._resolve_species([9606], db_dir=tmp_path, prefer_metadata=True)
    assert [t.name for t in taxa] == ["Homo sapiens"]


def test_databases_empty(tmp_path, monkeypatch):
    monkeypatch.setenv("FRANKENSEARCH_HOME", str(tmp_path))
    result = runner.invoke(app, ["databases"])
    assert result.exit_code == 0
    assert "No databases" in result.stdout


def test_databases_lists_built(tmp_path, monkeypatch):
    monkeypatch.setattr(uniprot, "build_query", lambda taxid, ps, *, client: "q")

    def fake_download(query, dest, *, client):
        Path(dest).write_text(">sp|P1|A_TEST alpha OS=Testus\nMKTAYIAKQRQISFV\n")
        return uniprot.DownloadResult(1, release="2026_02", release_date="10-June-2026")

    monkeypatch.setattr(uniprot, "download_fasta", fake_download)
    taxon = Taxon(9606, "Homo sapiens", "species")
    database.build(taxon, "swissprot", db_dir=tmp_path, client=None)

    result = runner.invoke(app, ["databases", "--db-dir", str(tmp_path)])
    assert result.exit_code == 0
    assert "Homo sapiens" in result.stdout
    assert "swissprot" in result.stdout
