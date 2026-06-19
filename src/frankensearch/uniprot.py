"""Download species proteomes from UniProt for building BLAST databases.

Three "sets" are supported per species:

* ``reference``  -- the species' reference proteome (looked up via the proteomes
  API, then downloaded by its proteome ID).
* ``swissprot``  -- reviewed (Swiss-Prot) entries only; small and high quality.
* ``all``        -- every UniProtKB entry for the organism.

All downloads use the ``uniprotkb/stream`` endpoint with ``compressed=true``;
httpx transparently decompresses the gzip response.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import httpx

from . import __version__
from .errors import UserError


@dataclass
class DownloadResult:
    sequence_count: int
    release: str | None = None
    release_date: str | None = None

USER_AGENT = f"frankensearch/{__version__}"
UNIPROTKB_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
PROTEOMES_SEARCH = "https://rest.uniprot.org/proteomes/search"

PROTEOME_SETS = ("reference", "swissprot", "all")

# Preference order when choosing among a taxon's proteomes (lower = better).
_PROTEOME_TYPE_PRIORITY = {
    "Reference proteome": 0,
    "Representative proteome": 1,
    "Other proteome": 2,
}


def make_client(timeout: float = 300.0) -> httpx.Client:
    return httpx.Client(timeout=timeout, headers={"User-Agent": USER_AGENT}, follow_redirects=True)


def build_query(taxid: int, proteome_set: str, *, client: httpx.Client) -> str:
    """Return the UniProtKB query string for a taxid + set."""
    if proteome_set == "swissprot":
        return f"(organism_id:{taxid}) AND (reviewed:true)"
    if proteome_set == "all":
        return f"(organism_id:{taxid})"
    if proteome_set == "reference":
        return f"(proteome:{find_reference_proteome(taxid, client=client)})"
    raise UserError(
        f"Unknown proteome set '{proteome_set}'.",
        hint=f"Choose one of: {', '.join(PROTEOME_SETS)}.",
    )


def find_reference_proteome(taxid: int, *, client: httpx.Client) -> str:
    """Find the best proteome ID (UPxxxxxxxxx) for a taxid."""
    try:
        response = client.get(
            PROTEOMES_SEARCH,
            params={"query": f"organism_id:{taxid}", "format": "json", "size": "25"},
        )
    except httpx.RequestError as exc:
        raise UserError(
            f"Could not reach UniProt to find a reference proteome for taxid {taxid}.",
            hint="Check your internet connection and try again.",
        ) from exc

    if response.status_code != 200:
        raise UserError(
            f"UniProt proteome lookup failed (status {response.status_code}) for taxid {taxid}."
        )

    candidates: list[tuple[int, str]] = []
    for item in response.json().get("results", []):
        priority = _PROTEOME_TYPE_PRIORITY.get(item.get("proteomeType", ""))
        upid = item.get("id")
        if priority is not None and upid:
            candidates.append((priority, upid))

    if not candidates:
        raise UserError(
            f"No usable reference proteome was found for taxid {taxid}.",
            hint="Try '--proteome-set swissprot' (reviewed only) or '--proteome-set all'.",
        )
    candidates.sort()
    return candidates[0][1]


def download_fasta(query: str, dest: Path, *, client: httpx.Client) -> DownloadResult:
    """Stream a FASTA download to ``dest``; return count + UniProt release info.

    Note: we deliberately do NOT pass ``compressed=true``. httpx negotiates gzip
    via ``Accept-Encoding`` and transparently decompresses the stream; the explicit
    ``compressed=true`` payload is *not* decoded by ``iter_bytes()`` and would write
    raw gzip to disk.
    """
    params = {"query": query, "format": "fasta"}
    sequence_count = 0
    seen_data = False
    release = release_date = None
    try:
        with client.stream("GET", UNIPROTKB_STREAM, params=params) as response:
            if response.status_code != 200:
                response.read()
                raise UserError(
                    f"UniProt download failed (status {response.status_code}).",
                    hint="Try again in a moment, or use a different --proteome-set.",
                )
            release = response.headers.get("x-uniprot-release")
            release_date = response.headers.get("x-uniprot-release-date")
            with open(dest, "wb") as handle:
                for chunk in response.iter_bytes():
                    if not chunk:
                        continue
                    if not seen_data:
                        if not chunk.lstrip().startswith(b">"):
                            raise UserError(
                                "UniProt did not return FASTA data.",
                                hint="Try again in a moment, or use a different --proteome-set.",
                            )
                        seen_data = True
                    handle.write(chunk)
                    sequence_count += chunk.count(b">")
    except httpx.RequestError as exc:
        raise UserError(
            "The proteome download was interrupted.",
            hint="Check your internet connection and re-run setup.",
        ) from exc
    return DownloadResult(sequence_count, release, release_date)
