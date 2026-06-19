"""Build and track per-species local BLAST databases.

Each species gets its own directory under the database root::

    <db_dir>/<taxid>/<taxid>.{pin,phr,psq,...}   # the BLAST database
    <db_dir>/<taxid>/metadata.json               # provenance for `doctor` / reuse

A database is considered "built" when both its metadata file and the BLAST index
files exist, so re-running `setup` skips work unless ``--force`` is given.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import paths, uniprot
from .blast import find_tool, tool_version
from .errors import DependencyError, UserError
from .taxonomy import Taxon

# Presence of either marks a built protein database (BLAST v5 writes both).
_DB_MARKER_EXTENSIONS = (".pin", ".pdb")


@dataclass
class DbMetadata:
    taxid: int
    scientific_name: str
    rank: str
    proteome_set: str
    query: str
    sequence_count: int
    built_at: str
    makeblastdb_version: str | None = None
    uniprot_release: str | None = None
    uniprot_release_date: str | None = None


def species_dir(taxid: int, db_dir: Path | None = None) -> Path:
    return (db_dir or paths.database_dir()) / str(taxid)


def db_prefix(taxid: int, db_dir: Path | None = None) -> Path:
    return species_dir(taxid, db_dir) / str(taxid)


def metadata_file(taxid: int, db_dir: Path | None = None) -> Path:
    return species_dir(taxid, db_dir) / "metadata.json"


def is_built(taxid: int, db_dir: Path | None = None) -> bool:
    if not metadata_file(taxid, db_dir).exists():
        return False
    prefix = db_prefix(taxid, db_dir)
    return any(Path(f"{prefix}{ext}").exists() for ext in _DB_MARKER_EXTENSIONS)


def is_current(taxid: int, proteome_set: str, db_dir: Path | None = None) -> bool:
    """True if a built database for this taxid exists AND was built with proteome_set."""
    if not is_built(taxid, db_dir):
        return False
    meta = load_metadata(taxid, db_dir)
    return meta is not None and meta.proteome_set == proteome_set


def load_metadata(taxid: int, db_dir: Path | None = None) -> DbMetadata | None:
    path = metadata_file(taxid, db_dir)
    if not path.exists():
        return None
    try:
        return DbMetadata(**json.loads(path.read_text()))
    except (json.JSONDecodeError, OSError, TypeError):
        return None


def list_built(db_dir: Path | None = None) -> list[DbMetadata]:
    root = db_dir or paths.database_dir()
    if not root.exists():
        return []
    metas: list[DbMetadata] = []
    for child in sorted(root.iterdir()):
        if child.is_dir() and child.name.isdigit():
            meta = load_metadata(int(child.name), db_dir)
            if meta is not None:
                metas.append(meta)
    return metas


def build(
    taxon: Taxon,
    proteome_set: str,
    *,
    db_dir: Path | None = None,
    force: bool = False,
    client,
) -> DbMetadata:
    """Download the proteome and build a BLAST database; reuse an existing one unless ``force``."""
    if not force and is_current(taxon.taxid, proteome_set, db_dir):
        existing = load_metadata(taxon.taxid, db_dir)
        if existing is not None:
            return existing

    makeblastdb = find_tool("makeblastdb")
    if makeblastdb is None:
        raise DependencyError(
            "makeblastdb (BLAST+) was not found on your PATH.",
            hint="Install BLAST+ or run 'conda activate frankensearch'.",
        )

    target_dir = species_dir(taxon.taxid, db_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    query = uniprot.build_query(taxon.taxid, proteome_set, client=client)

    handle = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=target_dir)
    fasta_path = Path(handle.name)
    handle.close()
    try:
        download = uniprot.download_fasta(query, fasta_path, client=client)
        if download.sequence_count == 0:
            raise UserError(
                f"UniProt returned no sequences for taxid {taxon.taxid} ({proteome_set}).",
                hint="Check the taxid, or try a different --proteome-set.",
            )
        _run_makeblastdb(
            makeblastdb, fasta_path, db_prefix(taxon.taxid, db_dir), taxon, proteome_set
        )
    finally:
        fasta_path.unlink(missing_ok=True)

    meta = DbMetadata(
        taxid=taxon.taxid,
        scientific_name=taxon.name,
        rank=taxon.rank,
        proteome_set=proteome_set,
        query=query,
        sequence_count=download.sequence_count,
        built_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        makeblastdb_version=tool_version("makeblastdb"),
        uniprot_release=download.release,
        uniprot_release_date=download.release_date,
    )
    metadata_file(taxon.taxid, db_dir).write_text(json.dumps(asdict(meta), indent=2))
    return meta


def _run_makeblastdb(
    executable: str, fasta: Path, prefix: Path, taxon: Taxon, proteome_set: str
) -> None:
    cmd = [
        executable,
        "-in", str(fasta),
        "-dbtype", "prot",
        "-blastdb_version", "5",
        "-title", f"{taxon.name} (taxid {taxon.taxid}) [{proteome_set}]",
        "-out", str(prefix),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except OSError as exc:
        raise UserError(f"Could not run makeblastdb: {exc}") from exc
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout).strip()
        raise UserError(
            f"makeblastdb failed while building the database for taxid {taxon.taxid}.",
            hint=detail[:400] if detail else None,
        )
