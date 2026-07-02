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
        query_seq=QUERY.sequence,
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


def test_write_tsv_output_query_adds_sequence_column(tmp_path):
    path = tmp_path / "out.tsv"
    output.write_tsv([_hit()], path, output_query=True)
    rows = list(csv.reader(path.open(), delimiter="\t"))
    assert rows[0] == output.tsv_columns(True)
    record = dict(zip(rows[0], rows[1], strict=True))
    assert record["query_sequence"] == QUERY.sequence
    # the column sits right after query_len
    assert rows[0].index("query_sequence") == rows[0].index("query_len") + 1


def test_write_txt_output_query_adds_query_seq_column(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [_hit()],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(output_query=True),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    # it's a table column (after Query), not a separate section
    assert "QUERY SEQUENCES" not in text
    header = next(line for line in text.splitlines() if line.startswith("Query"))
    assert "Query-Seq" in header
    assert header.index("Query-Seq") < header.index("Species")
    assert QUERY.sequence in text


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


def test_write_outputs_writes_seven_files(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    paths = output.write_outputs(
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
    tsv_path, txt_path, aln_path, summary_path, top1_tsv, top1_txt, top1_aln = paths
    assert tsv_path == tmp_path / "results.tsv" and tsv_path.exists()
    assert txt_path == tmp_path / "results.txt" and txt_path.exists()
    assert aln_path == tmp_path / "results_alignments.txt" and aln_path.exists()
    assert summary_path == tmp_path / "results.summary.md" and summary_path.exists()
    assert top1_tsv == tmp_path / "results_top1.tsv" and top1_tsv.exists()
    assert top1_txt == tmp_path / "results_top1.txt" and top1_txt.exists()
    assert top1_aln == tmp_path / "results_top1_alignments.txt" and top1_aln.exists()
    # The main .txt is table-only and points to its alignments file, which has them.
    main_text = txt_path.read_text()
    assert "Sbjct" not in main_text
    assert "results_alignments.txt" in main_text
    assert "Sbjct" in aln_path.read_text()
    # The _top1.txt is also table-only now, pointing to _top1_alignments.txt.
    top1_text = top1_txt.read_text()
    assert "best hit per query/species" in top1_text
    assert "Sbjct" not in top1_text
    assert "results_top1_alignments.txt" in top1_text
    assert "Sbjct" in top1_aln.read_text()


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
    # frank2 has no passing hit in any species -> no alignment block for it;
    # frank1 (which passes) keeps its block.
    assert "Query: frank2" not in text
    assert "Query: frank1" in text


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


def test_write_outputs_with_filter_writes_nine_files(tmp_path):
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
    assert len(paths) == 9
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
    assert len(paths) == 5  # no _alignments.txt or _top1_alignments.txt
    assert not (tmp_path / "results_alignments.txt").exists()
    assert not (tmp_path / "results_top1_alignments.txt").exists()
    # main .txt has no alignments and no pointer; _top1.txt is table-only too
    main_text = (tmp_path / "results.txt").read_text()
    assert "Sbjct" not in main_text and "_alignments.txt" not in main_text
    top1_text = (tmp_path / "results_top1.txt").read_text()
    assert "Sbjct" not in top1_text and "_alignments.txt" not in top1_text
    # the compact match column survives in the .tsv
    assert "match" in (tmp_path / "results.tsv").read_text().splitlines()[0]


def _meta_with_residues(residues=10_000_000):
    return DbMetadata(
        taxid=9606,
        scientific_name="Homo sapiens",
        rank="species",
        proteome_set="reference",
        query="(proteome:UP000005640)",
        sequence_count=20000,
        built_at="2026-06-29T00:00:00+00:00",
        residue_count=residues,
    )


def test_chance_match_length():
    assert output.chance_match_length(10_000_000, 77) == 7
    assert output.chance_match_length(0, 77) == 0  # graceful when M unknown
    assert output.chance_match_length(10_000_000, 0) == 0


def test_database_line_includes_k_star(tmp_path):
    path = tmp_path / "out.txt"
    output.write_txt(
        [_hit()],
        path,
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(),
        input_path=tmp_path / "in.fasta",
        db_metadata={9606: _meta_with_residues()},
    )
    text = path.read_text()
    assert "10,000,000 residues" in text
    assert "chance-match length k*≈" in text


def test_summary_chance_reference_and_truncation_note(tmp_path):
    path = tmp_path / "r.summary.md"
    output.write_summary(
        path,
        hits=[_hit()],
        queries=[QUERY],
        taxa=[TAXON],
        params=SearchParams(num_hits=10),
        input_path=tmp_path / "in.fasta",
        command="frankensearch search in.fasta --taxids 9606",
        frankensearch_version="0.7.0",
        blast_versions=BLAST_VERSIONS,
        db_metadata={9606: _meta_with_residues()},
        truncated_groups=3,
    )
    text = path.read_text()
    assert "## Chance reference" in text
    assert "10,000,000" in text and "Robinson" in text
    assert "3 (query, species) group(s) had more hits than the -n cap" in text


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


# --- Executive summary -----------------------------------------------------

MOUSE = Taxon(10090, "Mus musculus", "species")
FULL_Q = "ABCDEFGHIJKLMNOPQRST"  # 20 aa construct


def _jhit(
    motif="motif_1",
    *,
    seq="seq1",
    source="s.tsv",
    target_name="Prot X",
    accession="PX",
    taxid=9606,
    species="Homo sapiens",
    nident=10,
    align_len=10,
    qstart=1,
    qend=10,
    sstart=5,
    span=(1, 20),
):
    query_id = f"{source}|{seq}|{motif}|{span[0]}_{span[1]}"
    qseq = FULL_Q[qstart - 1 : qend]
    sseq = qseq[:nident] + "X" * (align_len - nident)
    return Hit(
        query_id=query_id,
        query_len=len(FULL_Q),
        taxid=taxid,
        species=species,
        subject_id=f"sp|{accession}|X",
        accession=accession,
        target_name=target_name,
        nident=nident,
        align_len=align_len,
        bitscore=50.0,
        evalue=1e-9,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=sstart + align_len - 1,
        qseq=qseq,
        sseq=sseq,
        query_seq=FULL_Q,
    )


def test_parse_junction_label():
    label = output.parse_junction_label("sample.tsv|seq2117|motif_10|936_974")
    assert label is not None
    assert (label.source, label.seq, label.motif) == ("sample.tsv", "seq2117", "motif_10")
    assert (label.start, label.end) == (936, 974)
    # Non-junction ids do not parse.
    assert output.parse_junction_label("frank1") is None
    assert output.parse_junction_label("a|b|c") is None  # only 3 fields
    assert output.parse_junction_label("a|b|motif_1|not_a_span") is None


def test_labels_are_junctions():
    assert output.labels_are_junctions([Query("s.tsv|seq1|motif_1|1_20", FULL_Q)])
    assert not output.labels_are_junctions([QUERY])  # "frank1"
    assert not output.labels_are_junctions([])


def test_write_exec_summary_collapses_and_shows_full_query(tmp_path):
    path = tmp_path / "r_executive_summary.txt"
    queries = [Query("s.tsv|seq1|motif_1|1_20", FULL_Q)]
    hits = [
        _jhit("motif_1", target_name="Prot X", accession="PX", nident=10),  # best, id 0.50
        _jhit("motif_1", target_name="Prot X", accession="PX", nident=8),  # dup target
        _jhit("motif_2", target_name="Prot Y", accession="PY", nident=6, span=(20, 40)),
    ]
    output.write_exec_summary(
        path, hits, queries=queries, taxa=[TAXON], params=SearchParams(exec_summary=True),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    assert "EXECUTIVE SUMMARY" in text
    assert "CONSTRUCT  seq1" in text
    assert "vs Homo sapiens" in text
    assert "STRONGEST MATCH" in text
    # Headline is the best Prot X hit (50% of the 20-aa query).
    assert "50% of the query is identical" in text
    # Collapsed by target: exactly one Prot X row and one Prot Y row in the list.
    body = text.split("TOP DISTINCT PROTEINS", 1)[1]
    assert body.count("Prot X") == 1
    assert "Prot Y" in body
    # The list carries an E-value column with each hit's value.
    assert "E-value" in body
    assert "1.0e-09" in body
    # Full-query alignment shows the un-aligned flank residues (11-20), not just the HSP.
    assert "KLMNOPQRST" in text


def test_exec_summary_lines_stay_within_100(tmp_path):
    # The document standardises on a 100-column width; a very long protein name and a
    # very long species name must be truncated so no rendered line exceeds it (both
    # the header "Species:" line and the per-construct "vs" line). The one exemption
    # is the echoed input path (data, not layout — like the Databases lines elsewhere).
    path = tmp_path / "r_executive_summary.txt"
    queries = [Query("s.tsv|seq1|motif_1|1_20", FULL_Q)]
    hits = [_jhit("motif_1", target_name="Very long protein name " * 6, accession="A0A0B4J1F4")]
    long_name = "Long organism scientific name with strain suffix " * 3  # ~145 chars
    long_species = Taxon(2697049, long_name, "species")
    output.write_exec_summary(
        path, hits, queries=queries, taxa=[TAXON, long_species],
        params=SearchParams(exec_summary=True), input_path=tmp_path / "in.fasta",
    )
    for line in path.read_text().splitlines():
        if line.startswith("Input:"):
            continue
        assert len(line) <= 100, f"line exceeds 100 chars ({len(line)}): {line!r}"


def test_write_exec_summary_per_species_and_no_hit(tmp_path):
    path = tmp_path / "r_executive_summary.txt"
    queries = [Query("s.tsv|seq1|motif_1|1_20", FULL_Q)]
    # A hit only in human; mouse should be reported as no-hit for this construct.
    hits = [_jhit("motif_1", taxid=9606, species="Homo sapiens")]
    output.write_exec_summary(
        path, hits, queries=queries, taxa=[TAXON, MOUSE], params=SearchParams(exec_summary=True),
        input_path=tmp_path / "in.fasta",
    )
    text = path.read_text()
    assert "vs Homo sapiens" in text
    assert "vs Mus musculus" in text
    assert "no hits in this species" in text


def test_full_query_alignment_covers_whole_query_with_flanks():
    # A partial match (query 3-12 of a 20-aa construct): the full query must still be
    # drawn, with the un-aligned flanks (1-2 and 13-20) on the Query line.
    hit = _jhit("motif_1", nident=10, align_len=10, qstart=3, qend=12, sstart=5)
    lines = output._full_query_alignment(hit)
    q_lines = [line for line in lines if line.startswith("Query")]
    assert q_lines[0].split()[1] == "1"  # numbering starts at the first query residue
    assert q_lines[-1].split()[-1] == "20"  # ...and ends at the last (whole query shown)
    assert "AB" in q_lines[0]  # left flank
    assert "MNOPQRST" in "\n".join(lines)  # right flank
    s_lines = [line for line in lines if line.startswith("Sbjct")]
    assert s_lines[0].split()[1] == "5"  # subject numbering starts at sstart


def test_full_query_alignment_uses_wide_lines_within_100():
    # A 70-residue alignment against a large-coordinate subject must fit on ONE line
    # (no 60/10 spill) and still stay within 100 columns once indented.
    seq70 = "ACDEFGHIKLMNPQRSTVWY" * 3 + "ACDEFGHIKL"  # 70 aa
    hit = Hit(
        query_id="s.tsv|seq1|motif_1|1_70", query_len=70, taxid=9606, species="Homo sapiens",
        subject_id="sp|Q8WZ42|TITIN", accession="Q8WZ42", target_name="Titin",
        nident=70, align_len=70, bitscore=100.0, evalue=1e-30,
        qstart=1, qend=70, sstart=34001, send=34070, qseq=seq70, sseq=seq70, query_seq=seq70,
    )
    lines = output._full_query_alignment(hit)
    query_lines = [line for line in lines if line.startswith("Query")]
    assert len(query_lines) == 1  # all 70 residues on one line
    # Indented into the block, the widest line still respects the 100-col budget.
    assert all(len("      " + line) <= 100 for line in lines)


def test_write_outputs_with_exec_summary_writes_extra_file(tmp_path):
    out_prefix = tmp_path / "r"
    queries = [Query("s.tsv|seq1|motif_1|1_20", FULL_Q)]
    written = output.write_outputs(
        [_jhit("motif_1")],
        out_prefix,
        queries=queries,
        taxa=[TAXON],
        params=SearchParams(exec_summary=True),
        input_path=tmp_path / "in.fasta",
        command="cmd",
        frankensearch_version="0.1.0",
        blast_versions=BLAST_VERSIONS,
        db_metadata={},
        top1_hits=[_jhit("motif_1")],
    )
    exec_path = output.exec_summary_output_path(out_prefix)
    assert exec_path in written
    assert exec_path.exists()
