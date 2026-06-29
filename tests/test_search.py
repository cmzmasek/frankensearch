"""Tests for the search engine: parsing, alignment rendering, and a real run."""

from pathlib import Path

import pytest

from frankensearch import database, uniprot
from frankensearch.errors import UserError
from frankensearch.inputs import Query
from frankensearch.search import (
    Hit,
    SearchParams,
    blastp_command,
    parse_subject,
    render_alignment,
    run_search,
    top1_of_group,
)
from frankensearch.taxonomy import Taxon

DB_FASTA = (
    ">sp|P1|A_TEST alpha test protein OS=Testus testus OX=1 GN=A\n"
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGK\n"
    ">sp|P2|B_TEST beta test protein OS=Testus testus OX=1 GN=B\n"
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDN\n"
)

TAXON = Taxon(9606, "Homo sapiens", "species")


def _mk_hit(acc, nident, align_len, query_len=20, bitscore=50.0):
    return Hit(
        query_id="q", query_len=query_len, taxid=9606, species="Homo sapiens",
        subject_id=f"sp|{acc}|X", accession=acc, target_name=acc,
        nident=nident, align_len=align_len, bitscore=bitscore, evalue=1e-5,
        qstart=1, qend=align_len, sstart=1, send=align_len,
        qseq="A" * align_len, sseq="A" * align_len,
    )


def _ranked(hits, rank_by):
    return sorted(hits, key=lambda h: h.sort_key(rank_by), reverse=True)


# --- unit tests (no BLAST) ------------------------------------------------- #
def test_top1_is_single_best_when_no_dead_heat():
    # Shorter 100% matches do NOT tie the longest 100% match: the ranking prefers
    # the longer alignment, so only the rank-1 hit is returned.
    rank_by = "identity-alignment"
    hits = _ranked(
        [_mk_hit("P1", 20, 20), _mk_hit("P2", 8, 8), _mk_hit("P3", 5, 5), _mk_hit("P4", 18, 20)],
        rank_by,
    )
    assert [h.accession for h in top1_of_group(hits, rank_by)] == ["P1"]


def test_top1_reports_genuine_dead_heat():
    # Two hits identical on identity, length, and bit score are a true tie; a
    # lower-identity hit is excluded.
    rank_by = "identity-alignment"
    hits = _ranked(
        [_mk_hit("A", 10, 10, bitscore=50.0), _mk_hit("B", 10, 10, bitscore=50.0),
         _mk_hit("C", 9, 10, bitscore=45.0)],
        rank_by,
    )
    assert {h.accession for h in top1_of_group(hits, rank_by)} == {"A", "B"}


def test_top1_alignment_length_mode():
    # Best by alignment length; the equal-length, equal-score hit ties, the
    # shorter one does not.
    rank_by = "alignment-length"
    hits = _ranked(
        [_mk_hit("A", 10, 30, bitscore=50.0), _mk_hit("B", 30, 30, bitscore=50.0),
         _mk_hit("C", 25, 25, bitscore=50.0)],
        rank_by,
    )
    assert {h.accession for h in top1_of_group(hits, rank_by)} == {"A", "B"}


def test_top1_empty_group():
    assert top1_of_group([], "identity-alignment") == []


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


def _hit(*, nident, align_len, query_len=100):
    return Hit(
        query_id="q", query_len=query_len, taxid=9606, species="Homo sapiens",
        subject_id="s", accession="A", target_name="t",
        nident=nident, align_len=align_len, bitscore=float(nident), evalue=1.0,
        qstart=1, qend=align_len, sstart=1, send=align_len,
        qseq="A" * align_len, sseq="A" * align_len,
    )


def test_sort_key_ranking_modes():
    short_exact = _hit(nident=10, align_len=10)   # 100% over a short alignment
    long_partial = _hit(nident=40, align_len=80)  # 50% over a long alignment
    hits = [short_exact, long_partial]

    by_aln_id = sorted(hits, key=lambda h: h.sort_key("identity-alignment"), reverse=True)
    assert by_aln_id[0] is short_exact  # higher identity wins

    by_query_id = sorted(hits, key=lambda h: h.sort_key("identity-query"), reverse=True)
    assert by_query_id[0] is long_partial  # 40/100 > 10/100 of the query

    by_len = sorted(hits, key=lambda h: h.sort_key("alignment-length"), reverse=True)
    assert by_len[0] is long_partial  # longest alignment wins


def test_render_alignment_marks_identities_and_coords():
    text = render_alignment("ABCD", "ABXD", qstart=1, sstart=1)
    lines = text.splitlines()
    assert lines[0].startswith("Query") and lines[0].rstrip().endswith("4")
    assert lines[1].startswith("Match")  # the middle line is labelled "Match"
    assert "AB.D" in lines[1]  # match line: identical letters, dot at the mismatch
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
    hits = run_search([query], [TAXON], params, db_dir=built_db).hits
    assert hits
    top = hits[0]
    assert top.query_id == "frank1"
    assert top.accession == "P1"
    assert top.identity_over_alignment == pytest.approx(1.0)
    # 22 identical residues over a 22-residue query
    assert top.identity_over_query == pytest.approx(22 / 22)


def test_run_search_respects_top_n(built_db):
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    hits = run_search([query], [TAXON], SearchParams(num_hits=1), db_dir=built_db).hits
    assert len(hits) <= 1


def test_run_search_top1_is_best_hit(built_db):
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    results = run_search([query], [TAXON], SearchParams(num_hits=5), db_dir=built_db)
    assert results.top1  # at least the best hit per (query, species)
    # every top1 hit ties the best identity ratio for its group
    best = max(h.identity_over_alignment for h in results.top1)
    assert all(h.identity_over_alignment == pytest.approx(best) for h in results.top1)


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
    results = run_search([query], [TAXON], SearchParams(remote=True), db_dir=tmp_path)
    assert results.hits == []
    assert results.top1 == []
