"""Write search results to a machine-readable .tsv and a human-readable .txt.

The .txt has two sections: a fixed-width table of every hit (one row per hit,
mirroring the console view plus query/subject coordinates) followed by the
multi-line BLAST-style pairwise alignments grouped by query then species. Both
files carry the alignment; the .tsv encodes it as a single field (newlines
escaped as ``\\n``).
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from . import scoring
from .database import DbMetadata
from .inputs import Query
from .search import Hit, SearchParams
from .taxonomy import Taxon

# Literature references included in the methods-grade summary.
_REFERENCES = [
    "Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and "
    "applications. BMC Bioinformatics. 2009;10:421.",
    "Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment "
    "search tool. J Mol Biol. 1990;215(3):403-410.",
    "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase. "
    "Nucleic Acids Res. (see https://www.uniprot.org for the current citation).",
]

TSV_COLUMNS = [
    "query_id",
    "query_len",
    "taxid",
    "scientific_name",
    "target_accession",
    "target_name",
    "subject_id",
    "nident",
    "alignment_length",
    "identity_ratio_alignment",
    "identity_ratio_query",
    "bit_score",
    "evalue",
    "q_start",
    "q_end",
    "s_start",
    "s_end",
    "alignment",
]


_RANK_LABELS = {
    "identity-alignment": "identity over alignment length",
    "identity-query": "identity over query length",
    "alignment-length": "alignment length",
}

# ASCII so the fixed-width .txt table aligns in every terminal (a Unicode arrow is
# East-Asian-ambiguous width and would shift the column in some locales).
RANK_MARK = "<"  # appended to the header of the column the hits are ranked by


@dataclass(frozen=True)
class HitColumn:
    """One column of the per-hit results table, shared by the console and .txt views."""

    header: str
    value: Callable[[Hit], str]
    numeric: bool = False  # right-aligned in the .txt table / console
    rank_flag: str | None = None  # the rank_by value that flags this column
    max_width: int | None = None  # truncate (.txt) / fold (console) long text
    in_console: bool = True  # shown in the compact console preview


# Single source of truth for the hit table. The console preview (cli._print_results)
# uses the in_console subset; the .txt table (_hit_table) uses them all.
HIT_COLUMNS: list[HitColumn] = [
    HitColumn("Query", lambda h: h.query_id),
    HitColumn("Species", lambda h: h.species),
    HitColumn("Accession", lambda h: h.accession),
    HitColumn("Target", lambda h: h.target_name, max_width=40),
    HitColumn("%id/aln", lambda h: f"{h.identity_over_alignment * 100:.1f}",
              numeric=True, rank_flag="identity-alignment"),
    HitColumn("%id/qry", lambda h: f"{h.identity_over_query * 100:.1f}",
              numeric=True, rank_flag="identity-query"),
    HitColumn("Aln len", lambda h: str(h.align_len),
              numeric=True, rank_flag="alignment-length"),
    HitColumn("Bits", lambda h: f"{h.bitscore:.0f}", numeric=True),
    HitColumn("E-value", lambda h: f"{h.evalue:.1e}", numeric=True),
    HitColumn("Qstart", lambda h: str(h.qstart), numeric=True, in_console=False),
    HitColumn("Qend", lambda h: str(h.qend), numeric=True, in_console=False),
    HitColumn("Sstart", lambda h: str(h.sstart), numeric=True, in_console=False),
    HitColumn("Send", lambda h: str(h.send), numeric=True, in_console=False),
    # Free-form, .txt-only; last so its variable width never shifts other columns.
    HitColumn("Match", lambda h: h.match_line, in_console=False),
]


def hit_column_header(col: HitColumn, rank_by: str) -> str:
    """Header text, with the ranking marker appended when this is the ranked column."""
    return col.header + (f" {RANK_MARK}" if col.rank_flag == rank_by else "")


def _scoring_summary(params: SearchParams) -> tuple[str, str, str, str]:
    """Return (backend, matrix_label, gaps, ranked_by) describing the effective run."""
    ranked_by = _RANK_LABELS.get(params.rank_by, _RANK_LABELS["identity-alignment"])
    if params.remote:
        backend = "NCBI remote (nr)"
        gaps = "ungapped" if params.ungapped else "NCBI defaults"
        matrix_label = params.matrix
        if params.matrix == "identity":
            matrix_label = "PAM30 (remote fallback from IDENTITY)"
    else:
        backend = "local databases"
        gaps = scoring.gap_description(params.matrix, ungapped=params.ungapped)
        matrix_label = params.matrix
    return backend, matrix_label, gaps, ranked_by


def output_paths(out_prefix: Path) -> tuple[Path, Path, Path]:
    """The three core output files for a run: (.tsv, .txt, .summary.md)."""
    base = out_prefix
    return (
        base.parent / f"{base.name}.tsv",
        base.parent / f"{base.name}.txt",
        base.parent / f"{base.name}.summary.md",
    )


def top1_output_paths(out_prefix: Path) -> tuple[Path, Path]:
    """The two best-hit-only files for a run: (_top1.tsv, _top1.txt)."""
    base = out_prefix
    return (
        base.parent / f"{base.name}_top1.tsv",
        base.parent / f"{base.name}_top1.txt",
    )


def write_outputs(
    hits: list[Hit],
    out_prefix: Path,
    *,
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
    command: str,
    frankensearch_version: str,
    blast_versions: dict[str, str | None],
    db_metadata: dict[int, DbMetadata],
    top1_hits: list[Hit],
) -> tuple[Path, Path, Path, Path, Path]:
    """Write the five output files; return their paths.

    Core: .tsv, .txt, .summary.md. Plus two best-hit-only views built from
    ``top1_hits`` (the single best hit per query/species, all ties included):
    _top1.tsv and _top1.txt.
    """
    tsv_path, txt_path, summary_path = output_paths(out_prefix)
    top1_tsv_path, top1_txt_path = top1_output_paths(out_prefix)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(hits, tsv_path)
    write_txt(
        hits,
        txt_path,
        queries=queries,
        taxa=taxa,
        params=params,
        input_path=input_path,
        db_metadata=db_metadata,
    )
    write_summary(
        summary_path,
        hits=hits,
        queries=queries,
        taxa=taxa,
        params=params,
        input_path=input_path,
        command=command,
        frankensearch_version=frankensearch_version,
        blast_versions=blast_versions,
        db_metadata=db_metadata,
    )
    write_tsv(top1_hits, top1_tsv_path)
    write_txt(
        top1_hits,
        top1_txt_path,
        queries=queries,
        taxa=taxa,
        params=params,
        input_path=input_path,
        db_metadata=db_metadata,
        top1=True,
    )
    return tsv_path, txt_path, summary_path, top1_tsv_path, top1_txt_path


def write_tsv(hits: list[Hit], path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for hit in hits:
            writer.writerow(
                [
                    hit.query_id,
                    hit.query_len,
                    hit.taxid,
                    hit.species,
                    hit.accession,
                    hit.target_name,
                    hit.subject_id,
                    hit.nident,
                    hit.align_len,
                    f"{hit.identity_over_alignment:.4f}",
                    f"{hit.identity_over_query:.4f}",
                    f"{hit.bitscore:.1f}",
                    f"{hit.evalue:.2e}",
                    hit.qstart,
                    hit.qend,
                    hit.sstart,
                    hit.send,
                    hit.alignment_text.replace("\n", "\\n"),
                ]
            )


def _database_lines(taxa: list[Taxon], db_metadata: dict[int, DbMetadata]) -> list[str]:
    """One header line per queried local database: name, taxid, count, set.

    Continuation lines are indented to align under the first entry, matching the
    12-column label padding of the surrounding header block.
    """
    entries: list[str] = []
    for taxon in taxa:
        meta = db_metadata.get(taxon.taxid)
        if meta is not None:
            entries.append(
                f"{meta.scientific_name} (taxid {meta.taxid}, "
                f"{meta.sequence_count:,} sequences, {meta.proteome_set})"
            )
        else:
            entries.append(f"{taxon.name} (taxid {taxon.taxid})")
    if not entries:
        return []
    indent = " " * 12
    return [f"Databases:  {entries[0]}"] + [f"{indent}{entry}" for entry in entries[1:]]


def write_txt(
    hits: list[Hit],
    path: Path,
    *,
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
    db_metadata: dict[int, DbMetadata] | None = None,
    top1: bool = False,
) -> None:
    by_group: dict[tuple[str, int], list[Hit]] = {}
    for hit in hits:
        by_group.setdefault((hit.query_id, hit.taxid), []).append(hit)

    backend, matrix_label, gaps, ranked_by = _scoring_summary(params)
    title = (
        "FRANKENSEARCH results — best hit per query/species (ties included)"
        if top1
        else "FRANKENSEARCH results"
    )
    selection = (
        "Selection:  best hit per (query, species); tied rank-1 hits all shown"
        if top1
        else f"Top hits per (query, species): {params.num_hits}"
    )
    lines: list[str] = [
        title,
        f"Input:      {input_path}",
        f"Generated:  {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Backend:    {backend}",
    ]
    if not params.remote:
        lines += _database_lines(taxa, db_metadata or {})
    lines += [
        f"Matrix:     {matrix_label}    Gaps: {gaps}",
        f"Ranked by:  {ranked_by}",
        selection,
        "Note:       E-value is shown for reference only; results are NOT filtered by E-value.",
    ]
    if params.remote:
        lines.append(
            "Note:       Remote nr is non-redundant; a hit's listed organism may "
            "differ from the queried taxid."
        )

    # Section 1: one table of every hit (grouped by query then species in row order).
    lines += ["", "=" * 80, "HITS", "=" * 80, ""]
    lines += _hit_table(hits, params.rank_by) if hits else ["(no hits)"]

    # Section 2: the BLAST-style pairwise alignments, grouped by query then species.
    lines += ["", "=" * 80, "ALIGNMENTS", "=" * 80]
    for query in queries:
        lines.append("")
        lines.append("-" * 80)
        lines.append(f"Query: {query.id}   (length {len(query.sequence)})")
        lines.append("-" * 80)
        for taxon in taxa:
            lines.append("")
            lines.append(f"--- {taxon.name}  (taxid {taxon.taxid}) ---")
            group = by_group.get((query.id, taxon.taxid), [])
            if not group:
                lines.append("    (no hits)")
                continue
            for rank, hit in enumerate(group, start=1):
                lines.extend(_hit_block(rank, hit))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _truncate(text: str, width: int) -> str:
    # ASCII "..." keeps the fixed-width .txt table aligned in every terminal.
    return text if len(text) <= width else text[: max(0, width - 3)] + "..."


def cell_text(col: HitColumn, hit: Hit) -> str:
    """The rendered string for one cell (truncated if the column caps its width)."""
    text = col.value(hit)
    return _truncate(text, col.max_width) if col.max_width is not None else text


def _hit_table(hits: list[Hit], rank_by: str) -> list[str]:
    """A fixed-width text table, one row per hit, driven by HIT_COLUMNS."""
    headers = [hit_column_header(col, rank_by) for col in HIT_COLUMNS]
    rows = [[cell_text(col, hit) for col in HIT_COLUMNS] for hit in hits]
    return _format_table(headers, rows, numeric=[col.numeric for col in HIT_COLUMNS])


def _format_table(headers: list[str], rows: list[list[str]], *, numeric: list[bool]) -> list[str]:
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt(cells: list[str]) -> str:
        # No rstrip: a full rectangular grid, so the rule underlines every column
        # (incl. the wide trailing Match column) rather than stopping short.
        return "  ".join(
            cell.rjust(widths[i]) if numeric[i] else cell.ljust(widths[i])
            for i, cell in enumerate(cells)
        )

    out = [fmt(headers), "  ".join("-" * w for w in widths)]
    out += [fmt(row) for row in rows]
    return out


def _hit_block(rank: int, hit: Hit) -> list[str]:
    block = [
        "",
        f"  [{rank}] {hit.accession}  {hit.target_name}",
        f"      identity (alignment): {hit.identity_over_alignment * 100:5.1f}%  "
        f"({hit.nident}/{hit.align_len})",
        f"      identity (query):     {hit.identity_over_query * 100:5.1f}%  "
        f"({hit.nident}/{hit.query_len})",
        f"      bit score: {hit.bitscore:.1f}   E-value: {hit.evalue:.2e}   "
        f"query {hit.qstart}-{hit.qend} / subject {hit.sstart}-{hit.send}",
        "",
    ]
    block.extend("      " + line for line in hit.alignment_text.splitlines())
    return block


def _sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return "(unavailable)"


def _clean_version(value: str | None) -> str | None:
    """Turn 'blastp: 2.17.0+' into '2.17.0+' for clean display."""
    if not value:
        return None
    return value.split(": ", 1)[1] if ": " in value else value


def _release_note(db_metadata: dict[int, DbMetadata]) -> str:
    """A ' (UniProt release X)' fragment when all databases share one release."""
    releases = sorted({m.uniprot_release for m in db_metadata.values() if m.uniprot_release})
    return f" (UniProt release {releases[0]})" if len(releases) == 1 else ""


def write_summary(
    path: Path,
    *,
    hits: list[Hit],
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
    command: str,
    frankensearch_version: str,
    blast_versions: dict[str, str | None],
    db_metadata: dict[int, DbMetadata],
) -> None:
    """Write a methods-grade Markdown summary for reproduction / publication."""
    backend, matrix_label, gaps, ranked_by = _scoring_summary(params)
    species = ", ".join(f"{t.name} (taxid {t.taxid})" for t in taxa)
    blastp_ver = _clean_version(blast_versions.get("blastp")) or "unknown"
    makeblastdb_ver = _clean_version(blast_versions.get("makeblastdb")) or "n/a (remote run)"

    lines = [
        "# FRANKENSEARCH analysis summary",
        "",
        f"- **Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- **FRANKENSEARCH version:** {frankensearch_version}",
        "",
        "## Command",
        "",
        "```",
        command,
        "```",
        "",
        "## Inputs",
        "",
        f"- Query file: `{input_path}` (format auto-detected)",
        f"- Queries searched: {len(queries)}",
        f"- Input SHA-256: `{_sha256(input_path)}`",
        f"- Hits reported: {len(hits)}",
        "",
        "## Target databases",
        "",
    ]

    if params.remote:
        lines.append(
            "Searched NCBI **nr** remotely, restricted per taxon via "
            "`txid<id>[ORGN]`. Note: nr is non-redundant, so a hit's listed "
            "organism may differ from the queried taxid."
        )
        lines.append("")
        lines.append(f"- Species queried: {species}")
    else:
        lines.append("Local UniProt databases (one BLAST database per species):")
        lines.append("")
        lines.append(
            "| Taxid | Species | Set | Sequences | UniProt release | UniProt query | Built |"
        )
        lines.append("|---|---|---|---|---|---|---|")
        for taxon in taxa:
            meta = db_metadata.get(taxon.taxid)
            if meta is not None:
                release = meta.uniprot_release or "?"
                if meta.uniprot_release and meta.uniprot_release_date:
                    release = f"{meta.uniprot_release} ({meta.uniprot_release_date})"
                lines.append(
                    f"| {meta.taxid} | {meta.scientific_name} | {meta.proteome_set} | "
                    f"{meta.sequence_count} | {release} | `{meta.query}` | {meta.built_at} |"
                )
            else:
                lines.append(f"| {taxon.taxid} | {taxon.name} | ? | ? | ? | ? | ? |")

    lines += [
        "",
        "## Search parameters",
        "",
        f"- Backend: {backend}",
        f"- Scoring matrix: {matrix_label}",
        f"- Gap costs: {gaps}",
        f"- Low-complexity (SEG) filter: {'on' if params.seg else 'off'}",
        "- Composition-based statistics: disabled (`-comp_based_stats 0`)",
        f"- Word size: {params.word_size}",
        f"- Max target sequences: {params.max_target_seqs}",
        f"- E-value cutoff: {params.evalue:g} "
        "(intentionally high — results are NOT filtered by E-value)",
        f"- Ranking: by {ranked_by} "
        "(both identity ratios and the alignment length are reported)",
        f"- Hits kept: top {params.num_hits} per (query, species)",
        "",
        "## Software versions",
        "",
        f"- FRANKENSEARCH {frankensearch_version}",
        f"- NCBI BLAST+ (blastp): {blastp_ver}",
        f"- makeblastdb: {makeblastdb_ver}",
        "- Sequence data: UniProt (https://www.uniprot.org)",
        "",
        "## References",
        "",
    ]
    lines += [f"- {ref}" for ref in _REFERENCES]
    lines += [
        "",
        "## Suggested methods text",
        "",
        "> "
        + (
            f"Short fusion (“franken”) protein sequences (n={len(queries)}) were "
            f"searched against {species} using FRANKENSEARCH v{frankensearch_version}, which "
            f"wraps NCBI BLAST+ (blastp {blastp_ver}). To detect chance-level similarity "
            "rather than homology, searches used "
            + (
                "the PAM30 matrix"
                if matrix_label.startswith("PAM30")
                else f"the {matrix_label} scoring matrix"
            )
            + " with composition-based statistics disabled and no E-value cutoff; hits were "
            f"ranked by {ranked_by}. The top {params.num_hits} hits per "
            "query per species were retained. "
            + (
                "Targets were NCBI nr restricted to each species' taxonomy ID."
                if params.remote
                else "Targets were the corresponding UniProt proteome for each species"
                + _release_note(db_metadata)
                + "; see the database table above for set, query, and download date."
            )
        ),
        "",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
