"""Resolve NCBI taxonomy IDs to scientific name + rank via UniProt.

We only need two things from taxonomy: a taxid's scientific name (for display and
output) and its rank (to warn when a taxid is not species-level). That is far too
little to justify downloading the full NCBI taxdump, so instead we make a tiny
request to UniProt's taxonomy endpoint:

    https://rest.uniprot.org/taxonomy/<taxid>.json

Each resolved taxid is cached in a small local JSON file, so repeat runs need no
network and results are reproducible.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import httpx

from . import __version__, paths
from .errors import UserError

DEFAULT_BASE_URL = "https://rest.uniprot.org/taxonomy"
CACHE_FILENAME = "taxonomy_cache.json"
_USER_AGENT = f"frankensearch/{__version__}"
_TIMEOUT = 30.0


class TaxonNotFound(UserError):
    """The taxid does not exist in the taxonomy (a user mistake)."""


class TaxonomyUnavailable(UserError):
    """Taxonomy could not be resolved right now (offline, not cached, server error)."""


@dataclass(frozen=True)
class Taxon:
    taxid: int
    name: str
    rank: str

    @property
    def is_species(self) -> bool:
        return self.rank.lower() == "species"


class TaxonomyResolver:
    """Resolve taxids to :class:`Taxon`, backed by a local JSON cache."""

    def __init__(
        self,
        cache_dir: Path | None = None,
        *,
        base_url: str = DEFAULT_BASE_URL,
        allow_network: bool = True,
        fetcher=None,
        client: httpx.Client | None = None,
        timeout: float = _TIMEOUT,
    ) -> None:
        self.cache_dir = cache_dir or paths.taxonomy_dir()
        self.cache_file = self.cache_dir / CACHE_FILENAME
        self.base_url = base_url.rstrip("/")
        self.allow_network = allow_network
        self.timeout = timeout
        self._fetcher = fetcher
        self._client = client
        self._owns_client = client is None
        self._cache: dict[str, dict[str, str]] = self._load_cache()
        self._dirty = False

    # -- public API -------------------------------------------------------- #
    @property
    def cache_size(self) -> int:
        return len(self._cache)

    def resolve(self, taxid: int) -> Taxon:
        cached = self._cache.get(str(taxid))
        if isinstance(cached, dict) and "name" in cached and "rank" in cached:
            return Taxon(taxid=taxid, name=cached["name"], rank=cached["rank"])
        # A missing/malformed entry is treated as a cache miss (re-fetched below).

        if not self.allow_network:
            raise TaxonomyUnavailable(
                f"Taxid {taxid} is not in the local cache and network lookups are disabled.",
                hint="Run 'frankensearch setup' with internet access to resolve and cache it.",
            )

        name, rank = self._fetch(taxid)
        self._cache[str(taxid)] = {"name": name, "rank": rank}
        self._dirty = True
        return Taxon(taxid=taxid, name=name, rank=rank)

    def resolve_many(self, taxids) -> list[Taxon]:
        result = [self.resolve(taxid) for taxid in taxids]
        self.flush()
        return result

    def flush(self) -> None:
        if not self._dirty:
            return
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self.cache_file.write_text(json.dumps(self._cache, indent=2, sort_keys=True))
            self._dirty = False
        except OSError:
            # A failure to persist the cache should never break the actual run.
            pass

    def close(self) -> None:
        if self._owns_client and self._client is not None:
            self._client.close()
            self._client = None

    def __enter__(self) -> TaxonomyResolver:
        return self

    def __exit__(self, *exc) -> None:
        self.flush()
        self.close()

    # -- internals --------------------------------------------------------- #
    def _load_cache(self) -> dict[str, dict[str, str]]:
        if self.cache_file.exists():
            try:
                data = json.loads(self.cache_file.read_text())
                if isinstance(data, dict):
                    return data
            except (json.JSONDecodeError, OSError):
                pass
        return {}

    def _get_client(self) -> httpx.Client:
        if self._client is None:
            self._client = httpx.Client(timeout=self.timeout, headers={"User-Agent": _USER_AGENT})
        return self._client

    def _fetch(self, taxid: int) -> tuple[str, str]:
        if self._fetcher is not None:
            return self._fetcher(taxid)

        url = f"{self.base_url}/{taxid}.json"
        try:
            response = self._get_client().get(url)
        except httpx.RequestError as exc:
            raise TaxonomyUnavailable(
                f"Could not reach UniProt to look up taxid {taxid}.",
                hint="Check your internet connection and try again.",
            ) from exc

        if response.status_code == 404:
            raise TaxonNotFound(
                f"Taxonomy ID {taxid} was not found.",
                hint="Verify the taxid at https://www.ncbi.nlm.nih.gov/taxonomy",
            )
        if response.status_code != 200:
            raise TaxonomyUnavailable(
                f"UniProt returned status {response.status_code} for taxid {taxid}.",
                hint="Try again in a moment.",
            )

        try:
            data = response.json()
        except ValueError as exc:
            raise TaxonomyUnavailable(
                f"UniProt returned an unreadable response for taxid {taxid}."
            ) from exc

        name = data.get("scientificName")
        if not name:
            raise TaxonomyUnavailable(
                f"UniProt response for taxid {taxid} did not include a scientific name."
            )
        return name, data.get("rank", "no rank")
