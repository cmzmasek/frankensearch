"""Tests for the search engine: parsing, alignment rendering, and a real run."""

from pathlib import Path

import pytest

from frankensearch import database, uniprot
from frankensearch.errors import UserError
from frankensearch.inputs import Query
from frankensearch.search import (
    Hit,
    SearchParams,
    _sort_key_from_parts,
    blastp_command,
    parse_subject,
    rank_metric,
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


def test_rank_metric_per_mode():
    hit = _mk_hit("A", nident=9, align_len=10, query_len=20)  # qry ratio 9/20=0.45
    assert rank_metric(hit, "identity-alignment") == pytest.approx(0.9)
    assert rank_metric(hit, "identity-query") == pytest.approx(0.45)
    assert rank_metric(hit, "alignment-length") == 10.0


def test_default_rank_by_is_identity_query():
    assert SearchParams().rank_by == "identity-query"


def test_run_search_warns_on_num_hits_truncation(monkeypatch):
    from frankensearch import search as search_module

    # two hits in one (query, species) group; -n 1 -> the group is truncated
    row_a = ["q0", "sp|A|X", "10", "10", "20", "1", "10", "1", "10", "50.0", "1e-5",
             "AAAAAAAAAA", "AAAAAAAAAA", "sp|A|X protein A"]
    row_b = ["q0", "sp|B|Y", "8", "10", "20", "1", "10", "1", "10", "40.0", "1e-3",
             "AAAAAAAACC", "AAAAAAAAAA", "sp|B|Y protein B"]
    monkeypatch.setattr(search_module, "_run_blastp", lambda *a, **k: [row_a, row_b])
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    warnings: list[str] = []
    results = run_search(
        [query], [TAXON], SearchParams(remote=True, num_hits=1),
        db_dir=Path("unused"), on_warning=warnings.append,
    )
    assert results.truncated_groups == 1
    assert len(results.hits) == 1
    assert any("num-hits" in w for w in warnings)


def _blast_row(sub, nident, qlen=20, align_len=10, bitscore=50.0, evalue="1e-5"):
    # One tabular row in _OUTFMT_FIELDS order (qseqid first, stitle last).
    return ["q0", f"sp|{sub}|X", str(nident), str(align_len), str(qlen), "1", str(align_len),
            "1", str(align_len), f"{bitscore}", evalue, "A" * align_len, "A" * align_len,
            f"sp|{sub}|X protein {sub}"]


def test_sort_key_from_parts_matches_hit_sort_key():
    # The cheap pre-build key must be bit-identical to Hit.sort_key in every mode,
    # or lazy construction could reorder results.
    row = _blast_row("Z", 7, qlen=20, align_len=12, bitscore=33.5, evalue="2.5e-4")
    hit = Hit(
        query_id="q", query_len=20, taxid=9606, species="Homo sapiens",
        subject_id="sp|Z|X", accession="Z", target_name="Z",
        nident=7, align_len=12, bitscore=33.5, evalue=2.5e-4,
        qstart=1, qend=12, sstart=1, send=12, qseq="A" * 12, sseq="A" * 12,
    )
    for mode in ("identity-query", "identity-alignment", "alignment-length"):
        assert _sort_key_from_parts(row, mode) == hit.sort_key(mode)


def test_run_search_ranks_best_first_through_heap(monkeypatch):
    # Fed out of rank order; the bounded top-N heap must still emit best-first.
    from frankensearch import search as search_module

    rows = [_blast_row("M", 7), _blast_row("H", 10), _blast_row("L", 4)]
    monkeypatch.setattr(search_module, "_run_blastp", lambda *a, **k: rows)
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    results = run_search(
        [query], [TAXON], SearchParams(remote=True, num_hits=5), db_dir=Path("unused")
    )
    assert [h.accession for h in results.hits] == ["H", "M", "L"]


def test_boundary_tie_keeps_earliest_like_stable_sort(monkeypatch):
    # Two hits (A, B) share the full sort key; a strictly-better hit (C) follows.
    # With -n 2, a stable descending sort + [:2] keeps the EARLIEST tied hit (A),
    # so the bounded heap must too -- [C, A], never [C, B].
    from frankensearch import search as search_module

    rows = [_blast_row("A", 5), _blast_row("B", 5), _blast_row("C", 6)]
    monkeypatch.setattr(search_module, "_run_blastp", lambda *a, **k: rows)
    query = Query("frank1", "M" * 20)
    results = run_search(
        [query], [TAXON], SearchParams(remote=True, num_hits=2), db_dir=Path("unused")
    )
    assert [h.accession for h in results.hits] == ["C", "A"]


def test_top1_keeps_full_dead_heat_beyond_num_hits(monkeypatch):
    # Three hits with an identical full sort key are a genuine rank-1 dead heat:
    # top1 keeps all three even though -n 1 caps the main results at one.
    from frankensearch import search as search_module

    rows = [_blast_row("A", 10), _blast_row("B", 10), _blast_row("C", 10)]
    monkeypatch.setattr(search_module, "_run_blastp", lambda *a, **k: rows)
    query = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
    results = run_search(
        [query], [TAXON], SearchParams(remote=True, num_hits=1), db_dir=Path("unused")
    )
    assert len(results.hits) == 1
    assert {h.accession for h in results.top1} == {"A", "B", "C"}


def test_lazy_build_matches_brute_force_reference(monkeypatch):
    # Differential check: the bounded/lazy run_search must yield exactly what a naive
    # "build every hit, sort, take top-N" reference would -- same hits in the same
    # order (ties included) and the same top1 -- across several -n values. Small value
    # ranges make full-key ties common, deliberately stressing boundary eviction.
    import random

    from frankensearch import search as search_module
    from frankensearch.search import _build_hit

    rng = random.Random(20260630)
    rows = []
    for i in range(300):
        qid = f"q{i % 3}"  # three queries share one taxon
        nident = rng.randint(1, 6)          # tiny ranges -> frequent full-key ties
        align_len = rng.randint(nident, 8)
        bitscore = float(rng.randint(1, 4))
        evalue = rng.choice(["1e-5", "1e-3"])
        rows.append([
            qid, f"sp|S{i:03d}|X", str(nident), str(align_len), "20",
            "1", str(align_len), "1", str(align_len), f"{bitscore}", evalue,
            "A" * align_len, "A" * align_len, f"sp|S{i:03d}|X p",
        ])
    monkeypatch.setattr(search_module, "_run_blastp", lambda *a, **k: rows)
    queries = [Query("q0", "M" * 20), Query("q1", "M" * 20), Query("q2", "M" * 20)]
    by_id = {q.id: q for q in queries}
    rank_by = "identity-query"

    for num_hits in (1, 2, 3, 5, 10):
        results = run_search(
            queries, [TAXON], SearchParams(remote=True, num_hits=num_hits), db_dir=Path("unused")
        )
        grouped: dict[str, list[Hit]] = {}
        for row in rows:
            grouped.setdefault(row[0], []).append(_build_hit(row, by_id[row[0]], TAXON))
        ref_hits, ref_top1 = [], []
        for q in queries:
            hits = grouped[q.id]
            hits.sort(key=lambda h: h.sort_key(rank_by), reverse=True)
            ref_hits.extend(hits[:num_hits])
            ref_top1.extend(top1_of_group(hits, rank_by))

        assert [h.subject_id for h in results.hits] == [h.subject_id for h in ref_hits]
        assert {h.subject_id for h in results.top1} == {h.subject_id for h in ref_top1}


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
