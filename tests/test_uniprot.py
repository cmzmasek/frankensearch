"""Tests for UniProt query building and download (no real network)."""

import httpx
import pytest

from frankensearch import uniprot
from frankensearch.errors import UserError


def _client(handler):
    return httpx.Client(transport=httpx.MockTransport(handler))


def test_build_query_swissprot():
    query = uniprot.build_query(9606, "swissprot", client=None)
    assert query == "(organism_id:9606) AND (reviewed:true)"


def test_build_query_all():
    assert uniprot.build_query(10090, "all", client=None) == "(organism_id:10090)"


def test_build_query_unknown_set():
    with pytest.raises(UserError):
        uniprot.build_query(1, "bogus", client=None)


def test_find_reference_proteome_prefers_reference():
    def handler(request):
        return httpx.Response(
            200,
            json={
                "results": [
                    {"id": "UP000005640", "proteomeType": "Reference proteome"},
                    {"id": "UP001055169", "proteomeType": "Excluded"},
                ]
            },
        )

    assert uniprot.find_reference_proteome(9606, client=_client(handler)) == "UP000005640"


def test_find_reference_proteome_none_usable_raises():
    def handler(request):
        return httpx.Response(200, json={"results": [{"id": "UPx", "proteomeType": "Excluded"}]})

    with pytest.raises(UserError):
        uniprot.find_reference_proteome(9606, client=_client(handler))


def test_build_query_reference_uses_lookup():
    def handler(request):
        return httpx.Response(
            200, json={"results": [{"id": "UP000000625", "proteomeType": "Reference proteome"}]}
        )

    query = uniprot.build_query(83333, "reference", client=_client(handler))
    assert query == "(proteome:UP000000625)"


def test_download_fasta_counts_and_writes(tmp_path):
    fasta = ">a\nMKT\n>b\nMQI\n"

    def handler(request):
        return httpx.Response(200, content=fasta.encode())

    dest = tmp_path / "out.fasta"
    result = uniprot.download_fasta("(organism_id:1)", dest, client=_client(handler))
    assert result.sequence_count == 2
    assert dest.read_text() == fasta


def test_download_fasta_captures_uniprot_release(tmp_path):
    def handler(request):
        return httpx.Response(
            200,
            content=b">a\nMKT\n",
            headers={"x-uniprot-release": "2026_02", "x-uniprot-release-date": "10-June-2026"},
        )

    result = uniprot.download_fasta("q", tmp_path / "o.fasta", client=_client(handler))
    assert result.sequence_count == 1
    assert result.release == "2026_02"
    assert result.release_date == "10-June-2026"


def test_download_fasta_bad_status_raises(tmp_path):
    def handler(request):
        return httpx.Response(500, text="boom")

    with pytest.raises(UserError):
        uniprot.download_fasta("q", tmp_path / "o.fasta", client=_client(handler))


def test_download_fasta_non_fasta_raises(tmp_path):
    def handler(request):
        return httpx.Response(200, text="<html>not fasta</html>")

    with pytest.raises(UserError):
        uniprot.download_fasta("q", tmp_path / "o.fasta", client=_client(handler))
