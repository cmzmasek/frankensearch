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
    # alignment is a single field with newlines escaped
    assert "\n" not in record["alignment"]
    assert "\\n" in record["alignment"]


def test_tsv_alignment_round_trips(tmp_path):
    hit = _hit()
    path = tmp_path / "out.tsv"
    output.write_tsv([hit], path)
    record = dict(zip(*list(csv.reader(path.open(), delimiter="\t")), strict=True))
    restored = record["alignment"].replace("\\n", "\n")
    assert restored == render_alignment(hit.qseq, hit.sseq, hit.qstart, hit.sstart)


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


def test_write_outputs_writes_three_files(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    tsv_path, txt_path, summary_path = output.write_outputs(
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
    )
    assert tsv_path == tmp_path / "results.tsv" and tsv_path.exists()
    assert txt_path == tmp_path / "results.txt" and txt_path.exists()
    assert summary_path == tmp_path / "results.summary.md" and summary_path.exists()


def test_write_summary_contains_methods_info(tmp_path):
    input_file = tmp_path / "in.fasta"
    input_file.write_text(">frank1\nMQIFVKTLTGKTITLEVEPSDT\n")
    meta = DbMetadata(
        taxid=9606,
        scientific_name="Homo sapiens",
        rank="species",
        proteome_set="swissprot",
        query="(organism_id:9606) AND (reviewed:true)",
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
    assert "(organism_id:9606) AND (reviewed:true)" in text
    assert "2026_02" in text  # UniProt release
    assert "Suggested methods text" in text
    assert "BLAST+" in text  # a reference is present
