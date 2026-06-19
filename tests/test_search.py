"""Tests for the search engine: parsing, alignment rendering, and a real run."""

from pathlib import Path

import pytest

from frankensearch import database, uniprot
from frankensearch.errors import UserError
from frankensearch.inputs import Query
from frankensearch.search import (
    SearchParams,
    blastp_command,
    parse_subject,
    render_alignment,
    run_search,
)
from frankensearch.taxonomy import Taxon

DB_FASTA = (
    ">sp|P1|A_TEST alpha test protein OS=Testus testus OX=1 GN=A\n"
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGK\n"
    ">sp|P2|B_TEST beta test protein OS=Testus testus OX=1 GN=B\n"
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDN\n"
)

TAXON = Taxon(9606, "Homo sapiens", "species")


# --- unit tests (no BLAST) ------------------------------------------------- #
def test_parse_subject_uniprot():
    acc, name = parse_subject(
        "sp|P62987|RL40_HUMAN",
        "sp|P62987|RL40_HUMAN Ubiquitin-ribosomal protein eL40 OS=Homo sapiens OX=9606 GN=UBA52",
    )
    assert acc == "P62987"
    assert name == "Ubiquitin-ribosomal protein eL40"


def test_parse_subject_non_uniprot_fallback():
    acc, name = parse_subject("weird_id", "weird_id some description")
    assert acc == "weird_id"
    assert name == "some description"


def test_blastp_command_local_defaults():
    cmd = blastp_command(Path("q.fa"), TAXON, SearchParams(matrix="identity"), None)
    assert "-remote" not in cmd
    assert "IDENTITY" in cmd
    assert "-seg" not in cmd  # off by default
    assert "-comp_based_stats" in cmd


def test_blastp_command_seg_on():
    cmd = blastp_command(Path("q.fa"), TAXON, SearchParams(seg=True), None)
    assert cmd[cmd.index("-seg") + 1] == "yes"


def test_blastp_command_remote_falls_back_to_pam30():
    cmd = blastp_command(Path("q.fa"), TAXON, SearchParams(remote=True, matrix="identity"), None)
    assert "-remote" in cmd and "nr" in cmd
    assert "txid9606[ORGN]" in cmd
    assert "PAM30" in cmd


def test_render_alignment_marks_identities_and_coords():
    text = render_alignment("ABCD", "ABXD", qstart=1, sstart=1)
    lines = text.splitlines()
    assert lines[0].startswith("Query") and lines[0].rstrip().endswith("4")
    assert "AB D" in lines[1]  # midline: identical letters, space at the mismatch
    assert lines[2].startswith("Sbjct")


# --- integration test (real makeblastdb + blastp, no network) -------------- #
@pytest.fixture
def built_db(tmp_path, monkeypatch):
    monkeypatch.setattr(uniprot, "build_query", lambda taxid, ps, *, client: "q")

    def fake_download(query, dest, *, client):
        Path(dest).write_text(DB_FASTA)
        return uniprot.DownloadResult(2)

    monkeypatch.setattr(uniprot, "download_fasta", fake_download)
    database.build(TAXON, "all", db_dir=tmp_path, client=None)
    return tmp_path


def test_run_search_finds_and_ranks(built_db):
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")  # perfect substring of P1
    params = SearchParams(num_hits=5)
    hits = run_search([query], [TAXON], params, db_dir=built_db)
    assert hits
    top = hits[0]
    assert top.query_id == "frank1"
    assert top.accession == "P1"
    assert top.identity_over_alignment == pytest.approx(1.0)
    # 22 identical residues over a 22-residue query
    assert top.identity_over_query == pytest.approx(22 / 22)


def test_run_search_respects_top_n(built_db):
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    hits = run_search([query], [TAXON], SearchParams(num_hits=1), db_dir=built_db)
    assert len(hits) <= 1


def test_run_search_warns_when_max_target_seqs_hit(built_db):
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    warnings: list[str] = []
    run_search(
        [query],
        [TAXON],
        SearchParams(num_hits=5, max_target_seqs=1),
        db_dir=built_db,
        on_warning=warnings.append,
    )
    assert any("max-target-seqs" in w for w in warnings)


def test_run_search_missing_db_raises(tmp_path):
    query = Query("frank1", "MQIFVKTLTGKT")
    with pytest.raises(UserError):
        run_search([query], [TAXON], SearchParams(), db_dir=tmp_path)


def test_remote_skips_local_db_check(tmp_path, monkeypatch):
    # Remote mode must not require a local database; stub the blastp call (no network).
    from frankensearch import search as search_module

    monkeypatch.setattr(search_module, "_run_blastp", lambda *a, **k: [])
    query = Query("frank1", "MQIFVKTLTGKT")
    # tmp_path has no DB built; remote=True should not raise about a missing DB.
    hits = run_search([query], [TAXON], SearchParams(remote=True), db_dir=tmp_path)
    assert hits == []
