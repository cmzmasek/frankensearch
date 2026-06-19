"""Run blastp per taxid and rank hits by identity (not E-value).

For each taxid we run one blastp call against that species' local database (all
query sequences at once), parse the tabular output, compute both identity ratios,
and keep the top-N hits per (query, taxid) ranked by the chosen ratio.
"""

from __future__ import annotations

import subprocess
import tempfile
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from . import database, scoring
from .blast import find_tool
from .errors import DependencyError, UserError
from .inputs import Query
from .taxonomy import Taxon

# Recognised NCBI/UniProt FASTA id prefixes of the form "<tag>|<accession>|...".
_DB_TAGS = {"sp", "tr", "gb", "dbj", "emb", "ref", "pdb", "pir", "prf", "tpg", "tpe", "tpd"}

# Remote NCBI runs: per-search timeout and a polite gap between submissions.
_REMOTE_TIMEOUT = 600
_REMOTE_DELAY = 3.0

# Tabular fields requested from blastp, in order. stitle is last (it may contain
# spaces, never tabs) so the split stays unambiguous.
_OUTFMT_FIELDS = [
    "qseqid", "sseqid", "nident", "length", "qlen",
    "qstart", "qend", "sstart", "send", "bitscore", "evalue",
    "qseq", "sseq", "stitle",
]


@dataclass
class SearchParams:
    matrix: str = "identity"
    ungapped: bool = False
    evalue: float = 200000.0
    word_size: int = 2
    max_target_seqs: int = 5000
    # How to rank hits: "identity-alignment", "identity-query", or "alignment-length".
    rank_by: str = "identity-alignment"
    num_hits: int = 10
    remote: bool = False
    seg: bool = False  # low-complexity (SEG) filtering; off by default


@dataclass
class Hit:
    query_id: str
    query_len: int
    taxid: int
    species: str
    subject_id: str
    accession: str
    target_name: str
    nident: int
    align_len: int
    bitscore: float
    evalue: float
    qstart: int
    qend: int
    sstart: int
    send: int
    qseq: str
    sseq: str

    @property
    def identity_over_alignment(self) -> float:
        return self.nident / self.align_len if self.align_len else 0.0

    @property
    def identity_over_query(self) -> float:
        return self.nident / self.query_len if self.query_len else 0.0

    def sort_key(self, rank_by: str) -> tuple:
        """Descending sort key for the chosen ranking mode (use reverse=True).

        Ties fall back to bit score, then to the (lower) E-value. For the two
        identity modes we first prefer the longer alignment, so a high-identity
        match over more residues outranks the same identity over fewer; in
        ``alignment-length`` mode that fallback is implied by the primary key.
        """
        if rank_by == "alignment-length":
            return (self.align_len, self.bitscore, -self.evalue)
        primary = (
            self.identity_over_query
            if rank_by == "identity-query"
            else self.identity_over_alignment
        )
        return (primary, self.align_len, self.bitscore, -self.evalue)

    @property
    def match_line(self) -> str:
        return match_string(self.qseq, self.sseq)

    @property
    def alignment_text(self) -> str:
        return render_alignment(self.qseq, self.sseq, self.qstart, self.sstart)


def run_search(
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    *,
    db_dir: Path | None = None,
    on_warning: Callable[[str], None] | None = None,
) -> list[Hit]:
    """Search all queries against every taxon's database; return ranked hits."""
    if find_tool("blastp") is None:
        raise DependencyError(
            "blastp (BLAST+) was not found on your PATH.",
            hint="Install BLAST+ or run 'conda activate frankensearch'.",
        )

    if not params.remote:
        missing = [t for t in taxa if not database.is_built(t.taxid, db_dir)]
        if missing:
            listing = ", ".join(f"{t.taxid} ({t.name})" for t in missing)
            raise UserError(
                f"No local database for: {listing}.",
                hint="Build it first:  frankensearch setup --taxids "
                + ",".join(str(t.taxid) for t in missing),
            )

    # Use safe synthetic IDs in the FASTA so query names with spaces/odd characters
    # can't be mangled by BLAST; map back to the real Query afterwards.
    id_map = {f"q{i}": query for i, query in enumerate(queries)}

    grouped: dict[tuple[str, int], list[Hit]] = {}
    with tempfile.TemporaryDirectory() as raw:
        query_file = Path(raw) / "queries.fasta"
        query_file.write_text("".join(f">{key}\n{q.sequence}\n" for key, q in id_map.items()))
        for index, taxon in enumerate(taxa):
            if params.remote and index > 0:
                time.sleep(_REMOTE_DELAY)  # be polite to NCBI between submissions
            for parts in _run_blastp(query_file, taxon, params, db_dir):
                query = id_map[parts[0]]
                hit = _build_hit(parts, query, taxon)
                grouped.setdefault((query.id, taxon.taxid), []).append(hit)

    results: list[Hit] = []
    capped = 0
    for query in queries:
        for taxon in taxa:
            hits = grouped.get((query.id, taxon.taxid), [])
            if len({h.subject_id for h in hits}) >= params.max_target_seqs:
                capped += 1
            hits.sort(key=lambda h: h.sort_key(params.rank_by), reverse=True)
            results.extend(hits[: params.num_hits])

    if capped and on_warning is not None:
        on_warning(
            f"{capped} query/species search(es) hit the --max-target-seqs limit "
            f"({params.max_target_seqs}); some high-identity hits may be missing. "
            "Re-run with a higher --max-target-seqs to be sure."
        )
    return results


def blastp_command(query_file: Path, taxon: Taxon, params: SearchParams, db_dir) -> list[str]:
    """Assemble the full blastp command for one taxon (local or remote)."""
    common = [
        "-comp_based_stats", "0",
        "-evalue", str(params.evalue),
        "-word_size", str(params.word_size),
        "-max_target_seqs", str(params.max_target_seqs),
    ]
    if params.seg:
        common += ["-seg", "yes"]  # off (blastp default) when not requested
    outfmt = ["-outfmt", "6 " + " ".join(_OUTFMT_FIELDS)]

    if params.remote:
        scoring_args, _ = scoring.remote_blast_args(params.matrix, ungapped=params.ungapped)
        return [
            "blastp", "-remote", "-db", "nr",
            "-query", str(query_file),
            "-entrez_query", f"txid{taxon.taxid}[ORGN]",
            *common, *scoring_args, *outfmt,
        ]
    return [
        "blastp",
        "-query", str(query_file),
        "-db", str(database.db_prefix(taxon.taxid, db_dir)),
        *common,
        *scoring.blast_args(params.matrix, ungapped=params.ungapped),
        *outfmt,
    ]


def _run_blastp(query_file: Path, taxon: Taxon, params: SearchParams, db_dir) -> list[list[str]]:
    cmd = blastp_command(query_file, taxon, params, db_dir)
    timeout: float | None = _REMOTE_TIMEOUT if params.remote else None
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise UserError(
            f"Remote BLAST timed out for taxid {taxon.taxid} ({taxon.name}).",
            hint="NCBI may be busy. Try again, use fewer taxids, or build local databases.",
        ) from exc
    except OSError as exc:
        raise DependencyError(f"Could not run blastp: {exc}") from exc
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout).strip()
        raise UserError(
            f"blastp failed for taxid {taxon.taxid} ({taxon.name}).",
            hint=detail[:300] if detail else None,
        )

    rows: list[list[str]] = []
    n = len(_OUTFMT_FIELDS)
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < n:
            continue
        if len(parts) > n:  # stitle (last) may, in theory, hold extra tabs
            parts = parts[: n - 1] + ["\t".join(parts[n - 1:])]
        rows.append(parts)
    return rows


def _build_hit(parts: list[str], query: Query, taxon: Taxon) -> Hit:
    (
        _qseqid, sseqid, nident, length, qlen,
        qstart, qend, sstart, send, bitscore, evalue,
        qseq, sseq, stitle,
    ) = parts
    accession, target_name = parse_subject(sseqid, stitle)
    return Hit(
        query_id=query.id,
        query_len=int(qlen),
        taxid=taxon.taxid,
        species=taxon.name,
        subject_id=sseqid,
        accession=accession,
        target_name=target_name,
        nident=int(nident),
        align_len=int(length),
        bitscore=float(bitscore),
        evalue=float(evalue),
        qstart=int(qstart),
        qend=int(qend),
        sstart=int(sstart),
        send=int(send),
        qseq=qseq,
        sseq=sseq,
    )


def parse_subject(sseqid: str, stitle: str) -> tuple[str, str]:
    """Extract (accession, protein name) from a UniProt-style subject id + title."""
    accession = sseqid
    parts = sseqid.split("|")
    if len(parts) >= 2 and parts[0] in _DB_TAGS and parts[1]:
        accession = parts[1]

    name = ""
    if stitle:
        after_id = stitle.split(None, 1)
        description = after_id[1] if len(after_id) > 1 else ""
        name = description.split(" OS=")[0].strip()
    return accession, name or accession


def match_string(qseq: str, sseq: str) -> str:
    """The "match" line: the query residue where query and subject are identical
    and a dot (``.``) at every non-identical column (incl. gaps), e.g. ``MKL.EV``.

    ``qseq`` and ``sseq`` are aligned strings and must be the same length (BLAST
    pairwise output always is); ``strict=True`` turns any mismatch into a loud
    error rather than a silently truncated match line.
    """
    return "".join(
        q if q == s and q != "-" else "." for q, s in zip(qseq, sseq, strict=True)
    )


def render_alignment(qseq: str, sseq: str, qstart: int, sstart: int, width: int = 60) -> str:
    """Render a BLAST-website-style pairwise alignment (query / match / subject)."""
    full_match = match_string(qseq, sseq)
    lines: list[str] = []
    qpos, spos = qstart, sstart
    for i in range(0, len(qseq), width):
        qchunk, schunk = qseq[i : i + width], sseq[i : i + width]
        midline = full_match[i : i + width]
        q_nongap = sum(1 for c in qchunk if c != "-")
        s_nongap = sum(1 for c in schunk if c != "-")
        q_end = qpos + q_nongap - 1 if q_nongap else qpos
        s_end = spos + s_nongap - 1 if s_nongap else spos
        lines.append(f"Query  {qpos:>6}  {qchunk}  {q_end}")
        lines.append(f"Match  {'':>6}  {midline}")
        lines.append(f"Sbjct  {spos:>6}  {schunk}  {s_end}")
        lines.append("")
        qpos += q_nongap
        spos += s_nongap
    return "\n".join(lines).rstrip()
