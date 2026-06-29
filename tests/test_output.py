"""Tests for the .tsv and .txt writers."""

import csv

from frankensearch import output
from frankensearch.database import DbMetadata
from frankensearch.inputs import Query
from frankensearch.search import Hit, SearchParams, render_alignment
from frankensearch.taxonomy import Taxon

BLAST_VERSIONS = {"blastp": "blastp: 2.17.0+", "makeblastdb": "makeblastdb: 2.17.0+"}

QUERY = Query("frank1", "MQIFVKTLTGKTITLEVEPSDT")
TAXON = Taxon(9606, "Homo sapiens", "species")


def _hit(qseq="MQIFVKT", sseq="MQIFVKT", nident=7, align_len=7, query_len=22):
    return Hit(
        query_id="frank1",
        query_len=query_len,
        taxid=9606,
        species="Homo sapiens",
        subject_id="sp|P62987|RL40_HUMAN",
        accession="P62987",
        target_name="Ubiquitin-ribosomal protein eL40",
        nident=nident,
        align_len=align_len,
        bitscore=178.0,
        evalue=2.0e-47,
        qstart=1,
        qend=len(qseq),
        sstart=1,
        send=len(sseq),
        qseq=qseq,
        sseq=sseq,
    )


def test_write_tsv_header_and_row(tmp_path):
    path = tmp_path / "out.tsv"
    output.write_tsv([_hit()], path)
    rows = list(csv.reader(path.open(), delimiter="\t"))
    assert rows[0] == output.TSV_COLUMNS
    record = dict(zip(rows[0], rows[1], strict=True))
    assert record["query_id"] == "frank1"
    assert record["target_accession"] == "P62987"
    assert record["identity_ratio_alignment"] == "1.0000"
    # the .tsv carries the compact match line, not the full alignment
    assert "alignment" not in record
    assert record["match"] == "MQIFVKT"  # all 7 residues identical


def test_write_alignments_txt_has_pairwise(tmp_path):
    hit = _hit()
    path = tmp_path / "out_alignments.txt"
    output.write_alignments_txt(
        [hit],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    assert "ALIGNMENTS" in text
    assert "Query" in text and "Sbjct" in text
    # the rendered alignment appears, indented as _hit_block writes it
    rendered = render_alignment(hit.qseq, hit.sseq, hit.qstart, hit.sstart)
    expected = "\n".join("      " + line for line in rendered.splitlines())
    assert expected in text


def test_write_txt_groups_and_includes_alignment(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [_hit()],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(num_hits=3),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    assert "Query: frank1" in text
    assert "Homo sapiens  (taxid 9606)" in text
    assert "P62987" in text
    assert "identity (alignment): 100.0%" in text
    assert "Query" in text and "Sbjct" in text


def test_write_txt_has_hit_table_with_coordinates(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [_hit(qseq="MQIFVKT", sseq="MQIFVKT")],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(num_hits=3, rank_by="identity-query"),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    # The HITS table carries the alignment-length and coordinate columns and marks
    # the ranked column.
    header = next(line for line in text.splitlines() if line.startswith("Query"))
    for col in ("Accession", "Aln len", "Qstart", "Qend", "Sstart", "Send", "Bits", "Match"):
        assert col in header
    assert header.rstrip().endswith("Match")  # match line is the last column
    assert "%id/qry <" in header  # ranked column flagged
    assert "%id/aln <" not in header and "Aln len <" not in header  # only one column flagged
    assert "ALIGNMENTS" in text


def test_hit_table_ranked_by_alignment_length_flags_aln_column(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [_hit()],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(num_hits=3, rank_by="alignment-length"),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    assert "Aln len <" in text  # alignment-length column flagged as the ranking column
    assert "%id/aln <" not in text and "%id/qry <" not in text  # only one column flagged
    assert "Ranked by:  alignment length" in text


def test_hit_table_match_column_shows_dots(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [_hit(qseq="KLVEV-SDT", sseq="KLAEVGSDT", nident=7, align_len=9)],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(num_hits=3),
        input_path=tmp_path / "in.fasta",
    )
    # The Match column carries the dotted match string in the (single) data row.
    hits_section = path.read_text().split("ALIGNMENTS")[0]
    assert any(line.rstrip().endswith("KL.EV.SDT") for line in hits_section.splitlines())


def test_write_txt_notes_taxon_with_no_hits(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(),
        input_path=tmp_path / "in.fasta",
    )
    assert "(no hits)" in path.read_text()


def test_write_txt_lists_queried_databases(tmp_path):
    path = tmp_path / "out.txt"
    meta = DbMetadata(
        taxid=9606,
        scientific_name="Homo sapiens",
        rank="species",
        proteome_set="reference",
        query="(proteome:UP000005640)",
        sequence_count=147519,
        built_at="2026-06-18T00:00:00+00:00",
    )
    output.write_txt(
        [_hit()],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(),
        input_path=tmp_path / "in.fasta",
        db_metadata={9606: meta},
    )
    text = path.read_text()
    assert "Databases:  Homo sapiens (taxid 9606, 147,519 sequences, reference)" in text


def test_write_txt_remote_omits_database_list(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(remote=True),
        input_path=tmp_path / "in.fasta",
        db_metadata={},
    )
    assert "Databases:" not in path.read_text()


def test_write_txt_remote_shows_backend_and_fallback(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(remote=True, matrix="identity"),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    assert "NCBI remote (nr)" in text
    assert "PAM30 (remote fallback" in text
    assert "non-redundant" in text


def test_write_outputs_writes_six_files(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    tsv_path, txt_path, aln_path, summary_path, top1_tsv, top1_txt = output.write_outputs(
        [_hit()],
        tmp_path / "results",
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(),
        input_path=input_file,
        command="frankensearch search in.fasta --taxids 9606",
        frankensearch_version="0.1.0",
        blast_versions=BLAST_VERSIONS,
        db_metadata={},
        top1_hits=[_hit()],
    )
    assert tsv_path == tmp_path / "results.tsv" and tsv_path.exists()
    assert txt_path == tmp_path / "results.txt" and txt_path.exists()
    assert aln_path == tmp_path / "results_alignments.txt" and aln_path.exists()
    assert summary_path == tmp_path / "results.summary.md" and summary_path.exists()
    assert top1_tsv == tmp_path / "results_top1.tsv" and top1_tsv.exists()
    assert top1_txt == tmp_path / "results_top1.txt" and top1_txt.exists()
    # The main .txt is table-only; the alignments live in the companion file.
    main_text = txt_path.read_text()
    assert "Sbjct" not in main_text
    assert "_alignments.txt" in main_text
    assert "Sbjct" in aln_path.read_text()
    # The _top1.txt identifies itself, keeps its alignments inline, and uses the
    # best-hit selection wording.
    top1_text = top1_txt.read_text()
    assert "best hit per query/species" in top1_text
    assert "Sbjct" in top1_text
    assert "Top hits per (query, species):" not in top1_text


FRANK2 = Query("frank2", "MKTAYIAKQRQISFVKSHF")


def test_write_filtered_tsv_status_and_no_hit(tmp_path):
    path = tmp_path / "f.tsv"
    # frank1 has a 100% (over-alignment) hit; frank2 has none in the results.
    passing, _ = output._passing_by_group(
        [_hit()], [QUERY, FRANK2], [TAXON], "identity-alignment", 0.5
    )
    output.write_filtered_tsv(path, queries=[QUERY, FRANK2], taxa=[TAXON], passing_by_group=passing)
    rows = list(csv.reader(path.open(), delimiter="\t"))
    assert rows[0] == ["status", *output.TSV_COLUMNS]
    by_query = {r[1]: r for r in rows[1:]}
    assert by_query["frank1"][0] == "hit"
    no_hit = by_query["frank2"]
    assert no_hit[0] == "no_hit_above_threshold"
    assert len(no_hit) == len(output.TSV_COLUMNS) + 1  # status + every column present
    # query/species are filled; the hit fields are blank.
    assert no_hit[1:5] == ["frank2", str(len(FRANK2.sequence)), "9606", "Homo sapiens"]
    assert all(cell == "" for cell in no_hit[5:])


def test_write_filtered_txt_lists_no_hit_queries(tmp_path):
    path = tmp_path / "f.txt"
    passing, no_hit = output._passing_by_group(
        [_hit()], [QUERY, FRANK2], [TAXON], "identity-alignment", 0.5
    )
    output.write_filtered_txt(
        path,
        queries=[QUERY, FRANK2],
        taxa=[TAXON],
        params=SearchParams(rank_by="identity-alignment", filter_by=0.5),
        input_path=tmp_path / "in.fasta",
        passing_by_group=passing,
        no_hit=no_hit,
    )
    text = path.read_text()
    assert "QUERIES WITH NO HIT ABOVE THRESHOLD" in text
    assert "frank2   Homo sapiens (taxid 9606)" in text
    assert "keep hits with identity over alignment length >= 0.5" in text


def test_write_filtered_txt_all_pass_message(tmp_path):
    path = tmp_path / "f.txt"
    passing, no_hit = output._passing_by_group(
        [_hit()], [QUERY], [TAXON], "identity-alignment", 0.5
    )
    output.write_filtered_txt(
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(rank_by="identity-alignment", filter_by=0.5),
        input_path=tmp_path / "in.fasta",
        passing_by_group=passing,
        no_hit=no_hit,
    )
    assert "every query has a hit above the threshold" in path.read_text()


def test_write_outputs_with_filter_writes_seven_files(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    paths = output.write_outputs(
        [_hit()],
        tmp_path / "results",
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(filter_by=0.5),
        input_path=input_file,
        command="frankensearch search in.fasta --taxids 9606 --filter-by 0.5",
        frankensearch_version="0.3.0",
        blast_versions=BLAST_VERSIONS,
        db_metadata={},
        top1_hits=[_hit()],
    )
    assert len(paths) == 8
    assert (tmp_path / "results_filtered_by_0.5.tsv").exists()
    assert (tmp_path / "results_filtered_by_0.5.txt").exists()


def test_write_outputs_no_alignments(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    paths = output.write_outputs(
        [_hit()],
        tmp_path / "results",
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(include_alignments=False),
        input_path=input_file,
        command="frankensearch search in.fasta --taxids 9606 --no-alignments",
        frankensearch_version="0.4.0",
        blast_versions=BLAST_VERSIONS,
        db_metadata={},
        top1_hits=[_hit()],
    )
    assert len(paths) == 5  # no _alignments.txt
    assert not (tmp_path / "results_alignments.txt").exists()
    # main .txt has no alignments and no pointer; _top1.txt is table-only too
    main_text = (tmp_path / "results.txt").read_text()
    assert "Sbjct" not in main_text and "_alignments.txt" not in main_text
    assert "Sbjct" not in (tmp_path / "results_top1.txt").read_text()
    # the compact match column survives in the .tsv
    assert "match" in (tmp_path / "results.tsv").read_text().splitlines()[0]


def test_write_summary_contains_methods_info(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    meta = DbMetadata(
        taxid=9606,
        scientific_name="Homo sapiens",
        rank="species",
        proteome_set="swissprot",
        query="(taxonomy_id:9606) AND (reviewed:true)",
        sequence_count=20442,
        built_at="2026-06-18T00:00:00+00:00",
        makeblastdb_version="makeblastdb: 2.17.0+",
        uniprot_release="2026_02",
        uniprot_release_date="10-June-2026",
    )
    path = tmp_path / "r.summary.md"
    output.write_summary(
        path,
        hits=[_hit()],
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(),
        input_path=input_file,
        command="frankensearch search in.fasta --taxids 9606",
        frankensearch_version="0.1.0",
        blast_versions=BLAST_VERSIONS,
        db_metadata={9606: meta},
    )
    text = path.read_text()
    assert "FRANKENSEARCH version:** 0.1.0" in text
    assert "frankensearch search in.fasta --taxids 9606" in text
    assert "Input SHA-256" in text
    assert "swissprot" in text
    assert "(taxonomy_id:9606) AND (reviewed:true)" in text
    assert "2026_02" in text  # UniProt release
    assert "Suggested methods text" in text
    assert "BLAST+" in text  # a reference is present
