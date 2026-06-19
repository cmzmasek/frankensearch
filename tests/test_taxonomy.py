"""Tests for taxonomy resolution and caching (no real network)."""

import json

import httpx
import pytest

from frankensearch import taxonomy
from frankensearch.taxonomy import (
    Taxon,
    TaxonNotFound,
    TaxonomyResolver,
    TaxonomyUnavailable,
)


def test_resolve_uses_fetcher_then_cache(tmp_path):
    calls: list[int] = []

    def fake(taxid):
        calls.append(taxid)
        return "Homo sapiens", "species"

    resolver = TaxonomyResolver(cache_dir=tmp_path, fetcher=fake)
    taxon = resolver.resolve(9606)
    assert taxon == Taxon(9606, "Homo sapiens", "species")
    assert taxon.is_species
    resolver.resolve(9606)  # second time should be served from cache
    assert calls == [9606]


def test_non_species_rank_flagged(tmp_path):
    resolver = TaxonomyResolver(cache_dir=tmp_path, fetcher=lambda t: ("Mammalia", "class"))
    assert resolver.resolve(40674).is_species is False


def test_cache_persists_across_instances(tmp_path):
    writer = TaxonomyResolver(cache_dir=tmp_path, fetcher=lambda t: ("Mus musculus", "species"))
    writer.resolve(10090)
    writer.flush()

    reader = TaxonomyResolver(cache_dir=tmp_path, allow_network=False)
    assert reader.resolve(10090).name == "Mus musculus"


def test_offline_and_uncached_raises(tmp_path):
    resolver = TaxonomyResolver(cache_dir=tmp_path, allow_network=False)
    with pytest.raises(TaxonomyUnavailable):
        resolver.resolve(9606)


def test_malformed_cache_entry_is_treated_as_miss(tmp_path):
    # Entry missing 'rank' must not crash with KeyError.
    (tmp_path / taxonomy.CACHE_FILENAME).write_text(json.dumps({"9606": {"name": "X"}}))

    offline = TaxonomyResolver(cache_dir=tmp_path, allow_network=False)
    with pytest.raises(TaxonomyUnavailable):  # treated as a miss, not KeyError
        offline.resolve(9606)

    online = TaxonomyResolver(cache_dir=tmp_path, fetcher=lambda t: ("Homo sapiens", "species"))
    assert online.resolve(9606).rank == "species"  # re-fetched cleanly


def test_not_found_propagates(tmp_path):
    def fake(taxid):
        raise TaxonNotFound(f"{taxid} not found")

    resolver = TaxonomyResolver(cache_dir=tmp_path, fetcher=fake)
    with pytest.raises(TaxonNotFound):
        resolver.resolve(99999999)


def test_http_fetch_parses_json(tmp_path):
    def handler(request):
        assert request.url.path.endswith("/9606.json")
        return httpx.Response(200, json={"scientificName": "Homo sapiens", "rank": "species"})

    client = httpx.Client(transport=httpx.MockTransport(handler))
    resolver = TaxonomyResolver(
        cache_dir=tmp_path, base_url="https://example.test/taxonomy", client=client
    )
    taxon = resolver.resolve(9606)
    assert taxon.name == "Homo sapiens"
    assert taxon.rank == "species"


def test_http_fetch_404_is_not_found(tmp_path):
    client = httpx.Client(transport=httpx.MockTransport(lambda req: httpx.Response(404)))
    resolver = TaxonomyResolver(
        cache_dir=tmp_path, base_url="https://example.test/taxonomy", client=client
    )
    with pytest.raises(TaxonNotFound):
        resolver.resolve(123)
