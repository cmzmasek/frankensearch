"""Run blastp per taxid and rank hits by identity (not E-value).

For each taxid we run one blastp call against that species' local database (all
query sequences at once), parse the tabular output, compute both identity ratios,
and keep the top-N hits per (query, taxid) ranked by the chosen ratio.
"""

from __future__ import annotations

import heapq
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
    # Default "identity-query": within a query group it orders by identical-residue
    # count (query_len is constant), which is monotonic with chance surprise.
    rank_by: str = "identity-query"
    num_hits: int = 10
    remote: bool = False
    seg: bool = False  # low-complexity (SEG) filtering; off by default
    # Optional threshold on the active rank_by metric for the _filtered output
    # files. Identity modes: a fraction 0-1; alignment-length: a residue count.
    filter_by: float | None = None
    # When False (--no-alignments), no pairwise alignments are written at all:
    # no _alignments.txt, and the _top1/_filtered .txt become table-only. The
    # compact match column/line is unaffected.
    include_alignments: bool = True
    # When True (--output-query), the full query sequence is added to the outputs:
    # a query_sequence column in the .tsv files and a QUERY SEQUENCES section in
    # the .txt files (the query name is always reported regardless).
    output_query: bool = False
    # When True (--exec-summary), also write a plain-language _executive_summary.txt
    # rolling up, per construct (seq) x species, the strongest match and top-5
    # distinct target proteins. Requires junction-style query ids; the CLI gates on
    # that and clears this flag (with a warning) when the labels don't match.
    exec_summary: bool = False


@dataclass
class SearchResults:
    """What a search produced.

    * ``hits``  -- the top-N (``num_hits``) hits per (query, taxid), the main result.
    * ``top1``  -- the single best hit per (query, taxid), plus any hits tied with
      it on the full ranking (a true dead heat); never cut by ``num_hits``.
    """

    hits: list[Hit]
    top1: list[Hit]
    # Number of (query, taxid) groups that had more hits than ``num_hits`` (so the
    # reported files show only a capped view) -- surfaced as a warning + summary note.
    truncated_groups: int = 0


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
    query_seq: str = ""  # the full query sequence (for --output-query)

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

        NOTE: ``_sort_key_from_parts`` recomputes this exact tuple straight from the
        raw blastp fields (to rank a row before building a Hit). Keep the two in sync.
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


def rank_metric(hit: Hit, rank_by: str) -> float:
    """The hit's value for the active ``--rank-by`` metric, in the same units that
    ``--filter-by`` uses: a fraction (0-1) for the identity modes, or the residue
    count for ``alignment-length``."""
    if rank_by == "alignment-length":
        return float(hit.align_len)
    if rank_by == "identity-query":
        return hit.identity_over_query
    return hit.identity_over_alignment  # identity-alignment (default)


def top1_of_group(sorted_hits: list[Hit], rank_by: str) -> list[Hit]:
    """From hits already sorted best-first, return the rank-1 hit plus any others
    the ranking genuinely cannot separate from it -- i.e. hits with the identical
    full sort key (the ``--rank-by`` metric, alignment length, bit score, and
    E-value). Usually exactly one hit; more only on a true dead heat."""
    if not sorted_hits:
        return []
    best_key = sorted_hits[0].sort_key(rank_by)
    tied: list[Hit] = []
    for hit in sorted_hits:
        if hit.sort_key(rank_by) == best_key:
            tied.append(hit)
        else:
            break  # sorted descending: once the key drops, nothing later ties
    return tied


def _sort_key_from_parts(parts: list[str], rank_by: str) -> tuple:
    """The ``Hit.sort_key`` tuple computed straight from a raw tabular row, so a hit
    can be ranked before the (comparatively costly) ``Hit`` is built. MUST stay
    bit-identical to ``Hit.sort_key`` -- same fields, same arithmetic. Indices follow
    ``_OUTFMT_FIELDS``: nident=2, length=3, qlen=4, bitscore=9, evalue=10."""
    nident = int(parts[2])
    align_len = int(parts[3])
    bitscore = float(parts[9])
    neg_evalue = -float(parts[10])
    if rank_by == "alignment-length":
        return (align_len, bitscore, neg_evalue)
    if rank_by == "identity-query":
        qlen = int(parts[4])
        primary = nident / qlen if qlen else 0.0
    else:  # identity-alignment
        primary = nident / align_len if align_len else 0.0
    return (primary, align_len, bitscore, neg_evalue)


class _GroupAccumulator:
    """Streaming top-N + rank-1 tie tracker for one (query, taxon) group.

    ``add`` is fed every row BLAST returns for the group but builds a ``Hit`` (and
    retains it) only for rows that can affect a kept result -- the ``num_hits`` best
    (a bounded min-heap keyed by the full sort key) and the set of hits sharing the
    single best key (the genuine rank-1 dead heat, for _top1). Rows that can't make
    either are ranked from their raw fields and dropped without ever building a Hit.
    It also counts the total rows (for the -n truncation warning) and detects when the
    distinct-subject count reaches ``max_target_seqs`` (freeing the tracking set once
    it does, since from then on only the flag matters). This keeps peak memory at the
    kept hits, not the tens of millions of HSPs a heavy short-query run can emit."""

    __slots__ = (
        "_num_hits", "_max_targets", "_heap", "_seq", "_subjects",
        "best_key", "tie", "total", "capped",
    )

    def __init__(self, num_hits: int, max_target_seqs: int):
        self._num_hits = num_hits
        self._max_targets = max_target_seqs
        self._heap: list[tuple] = []  # min-heap of (sort_key, -seq, hit); len <= num_hits
        self._seq = 0
        self._subjects: set[str] | None = set()
        self.best_key: tuple | None = None
        self.tie: list[Hit] = []
        self.total = 0
        self.capped = False

    def add(self, key: tuple, parts: list[str], query: Query, taxon: Taxon) -> None:
        """Offer one tabular row (pre-ranked ``key``) to the group; build and keep a
        ``Hit`` only if it belongs in the top-N heap or the rank-1 tie set."""
        self.total += 1
        if self._subjects is not None:
            self._subjects.add(parts[1])  # sseqid
            if len(self._subjects) >= self._max_targets:
                self.capped = True
                self._subjects = None  # only the flag is needed from here on

        new_best = self.best_key is None or key > self.best_key
        ties_best = not new_best and key == self.best_key
        enters_heap = self._num_hits > 0 and (
            len(self._heap) < self._num_hits or key > self._heap[0][0]
        )
        if not (new_best or ties_best or enters_heap):
            return  # row can't affect any kept result -> never build the Hit

        hit = _build_hit(parts, query, taxon)
        if new_best:
            self.best_key = key
            self.tie = [hit]
        elif ties_best:
            self.tie.append(hit)
        if enters_heap:
            # -seq is the heap tie-break: among equal keys the min entry (the one
            # evicted) is the LATEST-seen, so the heap keeps the earliest -- matching a
            # stable descending sort + [:num_hits]. It also makes entries unique, so
            # the hit itself is never compared. (top_hits restores order via -e[1].)
            entry = (key, -self._seq, hit)
            self._seq += 1
            if len(self._heap) < self._num_hits:
                heapq.heappush(self._heap, entry)
            else:
                heapq.heapreplace(self._heap, entry)

    def top_hits(self, rank_by: str) -> list[Hit]:
        """The retained hits, ranked best-first (output order preserved among ties)."""
        # e[1] is -seq; sort by seq ascending (= -e[1]) to recover BLAST output order,
        # then a stable descending sort by the rank key keeps ties in that order.
        hits = [entry[2] for entry in sorted(self._heap, key=lambda e: -e[1])]
        hits.sort(key=lambda h: h.sort_key(rank_by), reverse=True)
        return hits


def run_search(
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    *,
    db_dir: Path | None = None,
    on_warning: Callable[[str], None] | None = None,
) -> SearchResults:
    """Search all queries against every taxon's database; return ranked results."""
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
    rank_by = params.rank_by
    num_hits = params.num_hits

    # Accumulate per (query, taxid) while streaming blastp output, keeping only the
    # trimmed top-N and the rank-1 tie set per group -- never the full hit list. A
    # heavy run (many short queries x a high --max-target-seqs) can emit tens of
    # millions of HSPs; materialising them all before trimming is what made the
    # post-BLAST phase slow and memory-hungry.
    topn: dict[tuple[str, int], list[Hit]] = {}
    tied: dict[tuple[str, int], list[Hit]] = {}
    capped = 0
    truncated = 0
    with tempfile.TemporaryDirectory() as raw:
        query_file = Path(raw) / "queries.fasta"
        query_file.write_text("".join(f">{key}\n{q.sequence}\n" for key, q in id_map.items()))
        for index, taxon in enumerate(taxa):
            if params.remote and index > 0:
                time.sleep(_REMOTE_DELAY)  # be polite to NCBI between submissions
            # One accumulator per query for this taxon; dropped (with its transient
            # subject sets) before the next taxon, so peak memory stays at one taxon.
            groups: dict[str, _GroupAccumulator] = {}
            for parts in _run_blastp(query_file, taxon, params, db_dir):
                query = id_map[parts[0]]
                acc = groups.get(query.id)
                if acc is None:
                    acc = groups[query.id] = _GroupAccumulator(num_hits, params.max_target_seqs)
                acc.add(_sort_key_from_parts(parts, rank_by), parts, query, taxon)
            for query_id, acc in groups.items():
                topn[(query_id, taxon.taxid)] = acc.top_hits(rank_by)
                tied[(query_id, taxon.taxid)] = acc.tie
                if acc.total > num_hits:
                    truncated += 1
                if acc.capped:
                    capped += 1

    # Emit query-major (then taxon), matching the original output row ordering.
    results: list[Hit] = []
    top1: list[Hit] = []
    for query in queries:
        for taxon in taxa:
            key = (query.id, taxon.taxid)
            results.extend(topn.get(key, []))
            top1.extend(tied.get(key, []))

    if capped and on_warning is not None:
        on_warning(
            f"{capped} query/species search(es) hit the --max-target-seqs limit "
            f"({params.max_target_seqs}); some high-identity hits may be missing. "
            "Re-run with a higher --max-target-seqs to be sure."
        )
    if truncated and on_warning is not None:
        on_warning(
            f"{truncated} query/species group(s) had more than {num_hits} hit(s); "
            f"only the top {num_hits} are reported. Re-run with a higher -n/--num-hits "
            "to see them all."
        )
    return SearchResults(results, top1, truncated)


def blastp_command(
    query_file: Path, taxon: Taxon, params: SearchParams, db_dir, out_path: Path | None = None
) -> list[str]:
    """Assemble the full blastp command for one taxon (local or remote).

    With ``out_path`` set, blastp writes its tabular output to that file (``-out``)
    instead of stdout, so a large result set is streamed from disk rather than
    buffered in memory."""
    common = [
        "-comp_based_stats", "0",
        "-evalue", str(params.evalue),
        "-word_size", str(params.word_size),
        "-max_target_seqs", str(params.max_target_seqs),
    ]
    if out_path is not None:
        common += ["-out", str(out_path)]
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


def _run_blastp(query_file: Path, taxon: Taxon, params: SearchParams, db_dir):
    """Run blastp for one taxon and yield its tabular rows (already tab-split).

    blastp writes to a temp file (``-out``) which we stream line by line, so a huge
    result set is never held in memory as one string; yielding lets ``run_search``
    trim each group as it goes rather than collecting every HSP first."""
    timeout: float | None = _REMOTE_TIMEOUT if params.remote else None
    n = len(_OUTFMT_FIELDS)
    with tempfile.TemporaryDirectory() as tmp:
        out_path = Path(tmp) / "hits.tsv"
        cmd = blastp_command(query_file, taxon, params, db_dir, out_path)
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
        with open(out_path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < n:
                    continue
                if len(parts) > n:  # stitle (last) may, in theory, hold extra tabs
                    parts = parts[: n - 1] + ["\t".join(parts[n - 1:])]
                yield parts


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
        query_seq=query.sequence,
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
