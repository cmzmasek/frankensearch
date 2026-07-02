"""Write search results to a machine-readable .tsv and a human-readable .txt.

The main .tsv and .txt are kept compact (one row per hit; the .tsv carries a short
``match`` column, the .txt a fixed-width table) so they stay manageable for large,
many-query runs. The bulky multi-line BLAST-style pairwise alignments are written
to a separate ``_alignments.txt``. The small curated views (``_top1``, ``_filtered``)
keep their alignments inline.
"""

from __future__ import annotations

import csv
import hashlib
import math
import re
import statistics
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from . import scoring
from .database import DbMetadata
from .inputs import Query
from .search import Hit, SearchParams, rank_metric
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
    "match",
]


_RANK_LABELS = {
    "identity-alignment": "identity over alignment length",
    "identity-query": "identity over query length",
    "alignment-length": "alignment length",
}

# The "Selection:" header line for the best-hit views (_top1.txt and _top1_alignments.txt).
_TOP1_SELECTION = "Selection:  best hit per (query, species); tied rank-1 hits all shown"

# Background probability that two residues drawn from typical proteins are identical,
# p = sum(f_i^2) over the Robinson & Robinson (1991) amino-acid frequencies. Used as
# the null for the expected-by-chance match length k* = ln(M*Q) / ln(1/p).
_CHANCE_MATCH_PROB = 0.0598


def chance_match_length(db_residues: int, query_len: int, p: float = _CHANCE_MATCH_PROB) -> int:
    """The exact-match length expected ~once by chance: k* = ln(M*Q) / ln(1/p), rounded.

    A match longer than k* is above the chance background; shorter exact matches are
    expected for short queries against a whole proteome.
    """
    if db_residues <= 0 or query_len <= 0:
        return 0
    return round(math.log(db_residues * query_len) / math.log(1.0 / p))

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


def alignments_output_path(out_prefix: Path) -> Path:
    """The pairwise-alignments file: _alignments.txt (kept out of the main .txt)."""
    return out_prefix.parent / f"{out_prefix.name}_alignments.txt"


def top1_alignments_output_path(out_prefix: Path) -> Path:
    """The best-hit pairwise-alignments file: _top1_alignments.txt."""
    return out_prefix.parent / f"{out_prefix.name}_top1_alignments.txt"


def filtered_output_paths(out_prefix: Path, filter_by: float) -> tuple[Path, Path]:
    """The two threshold-filtered files: (_filtered_by_<x>.tsv, _filtered_by_<x>.txt).

    ``<x>`` is the threshold formatted with ``g`` (e.g. 0.8, 12), so distinct
    thresholds get distinct filenames and don't overwrite each other.
    """
    base = out_prefix
    tag = f"{filter_by:g}"
    return (
        base.parent / f"{base.name}_filtered_by_{tag}.tsv",
        base.parent / f"{base.name}_filtered_by_{tag}.txt",
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
    truncated_groups: int = 0,
) -> tuple[Path, ...]:
    """Write the output files; return every path written.

    Always: the compact (table-only) main .tsv/.txt and best-hit-only _top1.tsv/.txt
    (the single best hit per query/species, all ties included), .summary.md, plus the
    bulky pairwise alignments split into _alignments.txt and _top1_alignments.txt.
    With ``--no-alignments`` neither alignments file is written. When
    ``params.filter_by`` is set, also writes _filtered_by_<x>.tsv and
    _filtered_by_<x>.txt (passing hits only, with every query/species that has none
    reported as "no hit above threshold").
    """
    include_alignments = params.include_alignments
    output_query = params.output_query
    tsv_path, txt_path, summary_path = output_paths(out_prefix)
    alignments_path = alignments_output_path(out_prefix)
    top1_tsv_path, top1_txt_path = top1_output_paths(out_prefix)
    top1_alignments_path = top1_alignments_output_path(out_prefix)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(hits, tsv_path, output_query=output_query)
    # The main .txt is always table-only; it points to _alignments.txt unless that
    # file is skipped via --no-alignments.
    write_txt(
        hits,
        txt_path,
        queries=queries,
        taxa=taxa,
        params=params,
        input_path=input_path,
        db_metadata=db_metadata,
        include_alignments=False,
        alignments_file=alignments_path.name if include_alignments else None,
    )
    written = [tsv_path, txt_path]
    if include_alignments:
        write_alignments_txt(
            hits,
            alignments_path,
            queries=queries,
            taxa=taxa,
            params=params,
            input_path=input_path,
            db_metadata=db_metadata,
        )
        written.append(alignments_path)
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
        truncated_groups=truncated_groups,
    )
    # The _top1.txt is table-only too; its alignments go to _top1_alignments.txt.
    write_tsv(top1_hits, top1_tsv_path, output_query=output_query)
    write_txt(
        top1_hits,
        top1_txt_path,
        queries=queries,
        taxa=taxa,
        params=params,
        input_path=input_path,
        db_metadata=db_metadata,
        top1=True,
        include_alignments=False,
        alignments_file=top1_alignments_path.name if include_alignments else None,
    )
    written += [summary_path, top1_tsv_path, top1_txt_path]
    if include_alignments:
        write_alignments_txt(
            top1_hits,
            top1_alignments_path,
            queries=queries,
            taxa=taxa,
            params=params,
            input_path=input_path,
            db_metadata=db_metadata,
            title="FRANKENSEARCH results — pairwise alignments (best hit per query/species)",
            info_lines=[_TOP1_SELECTION],
        )
        written.append(top1_alignments_path)

    if params.filter_by is not None:
        filt_tsv_path, filt_txt_path = filtered_output_paths(out_prefix, params.filter_by)
        passing_by_group, no_hit = _passing_by_group(
            hits, queries, taxa, params.rank_by, params.filter_by
        )
        write_filtered_tsv(
            filt_tsv_path,
            queries=queries,
            taxa=taxa,
            passing_by_group=passing_by_group,
            output_query=output_query,
        )
        write_filtered_txt(
            filt_txt_path,
            queries=queries,
            taxa=taxa,
            params=params,
            input_path=input_path,
            passing_by_group=passing_by_group,
            no_hit=no_hit,
            db_metadata=db_metadata,
            include_alignments=include_alignments,
        )
        written += [filt_tsv_path, filt_txt_path]

    if params.exec_summary:
        exec_path = exec_summary_output_path(out_prefix)
        write_exec_summary(
            exec_path,
            hits,
            queries=queries,
            taxa=taxa,
            params=params,
            input_path=input_path,
        )
        written.append(exec_path)

    return tuple(written)


def tsv_columns(output_query: bool) -> list[str]:
    """The .tsv header. With ``output_query`` a ``query_sequence`` column is
    inserted right after ``query_len``."""
    if not output_query:
        return TSV_COLUMNS
    cols = list(TSV_COLUMNS)
    cols.insert(2, "query_sequence")
    return cols


def _hit_tsv_row(hit: Hit, include_query_seq: bool = False) -> list:
    """One .tsv data row for a hit, in ``tsv_columns()`` order. With
    ``include_query_seq`` the full query sequence is inserted after ``query_len``."""
    row = [
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
        hit.match_line,
    ]
    if include_query_seq:
        row.insert(2, hit.query_seq)
    return row


def write_tsv(hits: list[Hit], path: Path, *, output_query: bool = False) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(tsv_columns(output_query))
        for hit in hits:
            writer.writerow(_hit_tsv_row(hit, output_query))


def _database_lines(
    taxa: list[Taxon],
    db_metadata: dict[int, DbMetadata],
    queries: list[Query] | None = None,
) -> list[str]:
    """One header line per queried local database: name, taxid, count, set, and -- when
    the DB residue count is known -- the chance-match length k* (per species; computed
    at the median query length, since k* barely varies across queries).

    Continuation lines are indented to align under the first entry, matching the
    12-column label padding of the surrounding header block.
    """
    median_q = statistics.median(len(q.sequence) for q in queries) if queries else 0
    entries: list[str] = []
    for taxon in taxa:
        meta = db_metadata.get(taxon.taxid)
        if meta is None:
            entries.append(f"{taxon.name} (taxid {taxon.taxid})")
            continue
        entry = (
            f"{meta.scientific_name} (taxid {meta.taxid}, "
            f"{meta.sequence_count:,} sequences, {meta.proteome_set}"
        )
        if meta.residue_count and median_q:
            k = chance_match_length(meta.residue_count, int(median_q))
            entry += f"; {meta.residue_count:,} residues, chance-match length k*≈{k}"
        entries.append(entry + ")")
    if not entries:
        return []
    indent = " " * 12
    return [f"Databases:  {entries[0]}"] + [f"{indent}{entry}" for entry in entries[1:]]


def _txt_header(
    params: SearchParams,
    input_path: Path,
    taxa: list[Taxon],
    db_metadata: dict[int, DbMetadata] | None,
    *,
    title: str,
    info_lines: list[str],
    queries: list[Query] | None = None,
) -> list[str]:
    """The shared .txt header: title, run metadata, queried databases (with the
    chance-match length k*), scoring and ranking, the caller's ``info_lines``, and
    the standard notes."""
    backend, matrix_label, gaps, ranked_by = _scoring_summary(params)
    lines: list[str] = [
        title,
        f"Input:      {input_path}",
        f"Generated:  {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Backend:    {backend}",
    ]
    if not params.remote:
        lines += _database_lines(taxa, db_metadata or {}, queries)
    lines += [
        f"Matrix:     {matrix_label}    Gaps: {gaps}",
        f"Ranked by:  {ranked_by}",
        *info_lines,
        "Note:       E-value is reported for reference only; it is never used as a filter.",
    ]
    if params.remote:
        lines.append(
            "Note:       Remote nr is non-redundant; a hit's listed organism may "
            "differ from the queried taxid."
        )
    return lines


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
    include_alignments: bool = True,
    alignments_file: str | None = None,
) -> None:
    by_group: dict[tuple[str, int], list[Hit]] = {}
    for hit in hits:
        by_group.setdefault((hit.query_id, hit.taxid), []).append(hit)

    if top1:
        title = "FRANKENSEARCH results — best hit per query/species (ties included)"
        info_lines = [_TOP1_SELECTION]
    else:
        title = "FRANKENSEARCH results"
        info_lines = [f"Top hits per (query, species): {params.num_hits}"]
    lines = _txt_header(
        params, input_path, taxa, db_metadata, title=title, info_lines=info_lines, queries=queries
    )

    # Section 1: one table of every hit (grouped by query then species in row order).
    lines += ["", "=" * 80, "HITS", "=" * 80, ""]
    lines += _hit_table(hits, params.rank_by, params.output_query) if hits else ["(no hits)"]

    # Section 2: the BLAST-style pairwise alignments. The bulky main file omits
    # them; it points to the companion _alignments.txt unless that file is also
    # being skipped (--no-alignments).
    if include_alignments:
        lines += ["", "=" * 80, "ALIGNMENTS", "=" * 80]
        lines += _alignment_blocks(queries, taxa, by_group, "(no hits)")
    elif alignments_file:
        lines.append("")
        lines.append(f"Pairwise alignments are in the companion {alignments_file}.")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_alignments_txt(
    hits: list[Hit],
    path: Path,
    *,
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
    db_metadata: dict[int, DbMetadata] | None = None,
    title: str = "FRANKENSEARCH results — pairwise alignments",
    info_lines: list[str] | None = None,
) -> None:
    """The pairwise alignments for the given hits, in their own file so the
    table .txt stays compact for large, many-query runs."""
    by_group: dict[tuple[str, int], list[Hit]] = {}
    for hit in hits:
        by_group.setdefault((hit.query_id, hit.taxid), []).append(hit)

    lines = _txt_header(
        params,
        input_path,
        taxa,
        db_metadata,
        title=title,
        info_lines=info_lines or [f"Top hits per (query, species): {params.num_hits}"],
        queries=queries,
    )
    lines += ["", "=" * 80, "ALIGNMENTS", "=" * 80]
    lines += _alignment_blocks(
        queries, taxa, by_group, "(no hits)", show_query_seq=params.output_query
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _alignment_blocks(
    queries: list[Query],
    taxa: list[Taxon],
    group_lookup: dict[tuple[str, int], list[Hit]],
    empty_label: str,
    show_query_seq: bool = False,
) -> list[str]:
    """The pairwise-alignment blocks, grouped by query then species. Each
    (query, species) with no hits shows ``empty_label`` instead. With
    ``show_query_seq`` the full query sequence is printed under each query header
    (used by the table-less _alignments.txt, where there is no Query-Seq column)."""
    lines: list[str] = []
    for query in queries:
        lines.append("")
        lines.append("-" * 80)
        lines.append(f"Query: {query.id}   (length {len(query.sequence)})")
        if show_query_seq:
            lines.append(query.sequence)
        lines.append("-" * 80)
        for taxon in taxa:
            lines.append("")
            lines.append(f"--- {taxon.name}  (taxid {taxon.taxid}) ---")
            group = group_lookup.get((query.id, taxon.taxid), [])
            if not group:
                lines.append(f"    {empty_label}")
                continue
            for rank, hit in enumerate(group, start=1):
                lines.extend(_hit_block(rank, hit))
    return lines


def _passing_by_group(
    hits: list[Hit],
    queries: list[Query],
    taxa: list[Taxon],
    rank_by: str,
    threshold: float,
) -> tuple[dict[tuple[str, int], list[Hit]], list[tuple[Query, Taxon]]]:
    """Split hits by (query, taxid), keeping only those whose ``rank_by`` metric
    is ``>= threshold``. Returns (passing_by_group, no_hit_combos), where
    no_hit_combos lists every (query, taxon) with no passing hit -- the whole
    point of the filter: surfacing queries that match nothing."""
    by_group: dict[tuple[str, int], list[Hit]] = {}
    for hit in hits:
        if rank_metric(hit, rank_by) >= threshold:
            by_group.setdefault((hit.query_id, hit.taxid), []).append(hit)
    no_hit = [(q, t) for q in queries for t in taxa if not by_group.get((q.id, t.taxid))]
    return by_group, no_hit


def write_filtered_tsv(
    path: Path,
    *,
    queries: list[Query],
    taxa: list[Taxon],
    passing_by_group: dict[tuple[str, int], list[Hit]],
    output_query: bool = False,
) -> None:
    """The .tsv of hits that passed the filter (``passing_by_group``), with a
    leading ``status`` column. Each (query, species) with no passing hit gets one
    row flagged ``no_hit_above_threshold`` (hit fields blank)."""
    blanks = [""] * (len(TSV_COLUMNS) - 4)  # query_id/len/taxid/species are filled
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["status", *tsv_columns(output_query)])
        for query in queries:
            for taxon in taxa:
                passing = passing_by_group.get((query.id, taxon.taxid), [])
                if passing:
                    for hit in passing:
                        writer.writerow(["hit", *_hit_tsv_row(hit, output_query)])
                else:
                    prefix = [query.id, len(query.sequence)]
                    if output_query:
                        prefix.append(query.sequence)
                    prefix += [taxon.taxid, taxon.name]
                    writer.writerow(["no_hit_above_threshold", *prefix, *blanks])


def write_filtered_txt(
    path: Path,
    *,
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
    passing_by_group: dict[tuple[str, int], list[Hit]],
    no_hit: list[tuple[Query, Taxon]],
    db_metadata: dict[int, DbMetadata] | None = None,
    include_alignments: bool = True,
) -> None:
    """The .txt of hits that passed the filter, with an explicit section listing
    every (query, species) in ``no_hit`` -- those with no hit above the threshold."""
    _, _, _, ranked_by = _scoring_summary(params)
    passing_hits = [
        h for q in queries for t in taxa for h in passing_by_group.get((q.id, t.taxid), [])
    ]

    info_lines = [
        f"Top hits per (query, species): {params.num_hits}",
        f"Filter:     keep hits with {ranked_by} >= {params.filter_by:g}",
    ]
    lines = _txt_header(
        params,
        input_path,
        taxa,
        db_metadata,
        title="FRANKENSEARCH results — filtered",
        info_lines=info_lines,
        queries=queries,
    )

    # Section 1: the hits that pass the threshold.
    lines += ["", "=" * 80, "HITS ABOVE THRESHOLD", "=" * 80, ""]
    lines += (
        _hit_table(passing_hits, params.rank_by, params.output_query)
        if passing_hits
        else ["(none)"]
    )

    # Section 2: the point of the filter -- query/species combos with no hit.
    lines += ["", "=" * 80, "QUERIES WITH NO HIT ABOVE THRESHOLD", "=" * 80, ""]
    if no_hit:
        for query, taxon in no_hit:
            lines.append(f"{query.id}   {taxon.name} (taxid {taxon.taxid})")
    else:
        lines.append("(every query has a hit above the threshold in every species)")

    # Section 3: pairwise alignments for the passing hits (skipped with --no-alignments).
    # Queries with no passing hit in any species are omitted here entirely -- they are
    # already listed in the no-hit section above.
    if include_alignments:
        lines += ["", "=" * 80, "ALIGNMENTS", "=" * 80]
        queries_with_hits = [
            q for q in queries if any(passing_by_group.get((q.id, t.taxid)) for t in taxa)
        ]
        lines += _alignment_blocks(
            queries_with_hits, taxa, passing_by_group, "(no hit above threshold)"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _truncate(text: str, width: int) -> str:
    # ASCII "..." keeps the fixed-width .txt table aligned in every terminal.
    return text if len(text) <= width else text[: max(0, width - 3)] + "..."


def cell_text(col: HitColumn, hit: Hit) -> str:
    """The rendered string for one cell (truncated if the column caps its width)."""
    text = col.value(hit)
    return _truncate(text, col.max_width) if col.max_width is not None else text


def _hit_table(hits: list[Hit], rank_by: str, output_query: bool = False) -> list[str]:
    """A fixed-width text table, one row per hit, driven by HIT_COLUMNS. With
    ``output_query`` a full-sequence ``Query-Seq`` column is inserted after Query."""
    columns = list(HIT_COLUMNS)
    if output_query:
        columns.insert(1, HitColumn("Query-Seq", lambda h: h.query_seq, in_console=False))
    headers = [hit_column_header(col, rank_by) for col in columns]
    rows = [[cell_text(col, hit) for col in columns] for hit in hits]
    return _format_table(headers, rows, numeric=[col.numeric for col in columns])


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
    truncated_groups: int = 0,
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

    chance_metas = [
        db_metadata[taxon.taxid]
        for taxon in taxa
        if db_metadata.get(taxon.taxid) and db_metadata[taxon.taxid].residue_count
    ]
    if chance_metas and queries:
        median_q = int(statistics.median(len(q.sequence) for q in queries))
        lines += [
            "",
            "## Chance reference",
            "",
            "Expected-once-by-chance exact-match length k\\* = ln(M·Q) / ln(1/p), with "
            "background match probability p = 0.0598 (Σf_i² over the Robinson & Robinson "
            "(1991) amino-acid frequencies), M = total database residues, and Q = query "
            f"length (median {median_q} aa here; k\\* varies by <1 residue across queries). "
            "Matches longer than k\\* are above the chance background; shorter exact matches "
            "are expected for short queries.",
            "",
            "| Taxid | Species | DB residues (M) | k\\* |",
            "|---|---|---|---|",
        ]
        for meta in chance_metas:
            k = chance_match_length(meta.residue_count, median_q)
            lines.append(
                f"| {meta.taxid} | {meta.scientific_name} | {meta.residue_count:,} | {k} |"
            )

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
        *(
            [
                f"- Note: {truncated_groups} (query, species) group(s) had more hits than "
                f"the -n cap ({params.num_hits}); only the top {params.num_hits} are in "
                "these files. Re-run with a higher -n/--num-hits to capture them all."
            ]
            if truncated_groups
            else []
        ),
        *(
            [
                f"- Filter (`_filtered_by_{params.filter_by:g}` files): keep hits with "
                f"{ranked_by} >= {params.filter_by:g}; query/species with none are "
                'reported as "no hit above threshold"'
            ]
            if params.filter_by is not None
            else []
        ),
        *(
            [
                "- Executive summary (`_executive_summary.txt`): per-construct strongest "
                "match + up to 5 distinct target proteins per species (requires "
                "junction-style query IDs)"
            ]
            if params.exec_summary
            else []
        ),
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
                f"Hits with {ranked_by} below {params.filter_by:g} were filtered out in "
                "a separate report, and queries with no hit at or above that threshold "
                "were recorded. "
                if params.filter_by is not None
                else ""
            )
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


# ---------------------------------------------------------------------------
# Executive summary
#
# A plain-language roll-up for a non-technical reader (e.g. a time-poor
# executive). Only meaningful when query ids carry the junction structure that
# extract_junctions.py emits: "<source_file>|<seq>|motif_<N>|<start>_<end>".
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JunctionLabel:
    """The four fields of an ``extract_junctions.py`` query id, e.g.
    ``sample.tsv|seq2117|motif_10|936_974``."""

    source: str
    seq: str
    motif: str
    start: int
    end: int


_JUNCTION_SPAN = re.compile(r"^(\d+)_(\d+)$")

# The executive summary standardises on a 100-column width (wider than the classic
# 80 so the alignment blocks stay readable). Rules span this width and the lines we
# render stay within it: long protein and species names are truncated to fit.
# Accessions are the one exception -- they are shown in full (never truncated), so a
# hit with an unusually long accession can run slightly over.
_EXEC_WIDTH = 100

# Residues per alignment line. An alignment row is
#   [6-space indent]"Query  "[6-col coord]"  "[residues]"  "[end coord]
# so the overhead is 6+7+6+2+2 = 23 plus the trailing coordinate. Budgeting 7 for
# that coordinate (no real protein exceeds ~35 kaa = 5 digits) leaves 100-30 = 70,
# so a typical junction-peptide alignment fits on a single line instead of spilling
# a few residues onto a second.
_EXEC_ALN_WIDTH = _EXEC_WIDTH - 30


def parse_junction_label(query_id: str) -> JunctionLabel | None:
    """Parse an ``extract_junctions.py`` query id, or return ``None`` when it does
    not match the 4-field ``file|seq|motif_N|start_end`` shape -- so callers can
    skip the executive summary rather than emit nonsense."""
    parts = query_id.split("|")
    if len(parts) != 4:
        return None
    source, seq, motif, span = parts
    span_match = _JUNCTION_SPAN.match(span)
    if not (source and seq and motif.startswith("motif") and span_match):
        return None
    return JunctionLabel(source, seq, motif, int(span_match.group(1)), int(span_match.group(2)))


def labels_are_junctions(queries: list[Query]) -> bool:
    """True only when every query id parses as a junction label -- the precondition
    for producing an executive summary."""
    return bool(queries) and all(parse_junction_label(q.id) is not None for q in queries)


def exec_summary_output_path(out_prefix: Path) -> Path:
    """The executive-summary file: _executive_summary.txt."""
    return out_prefix.parent / f"{out_prefix.name}_executive_summary.txt"


def _fit_species(prefix: str, taxon: Taxon, gap: str = " ") -> str:
    """A species line ``<prefix><name><gap>(taxid N)`` with the name truncated so the
    whole line stays within the 100-column width even for long scientific names."""
    suffix = f"{gap}(taxid {taxon.taxid})"
    name = _truncate(taxon.name, _EXEC_WIDTH - len(prefix) - len(suffix))
    return f"{prefix}{name}{suffix}"


def _exec_sort_key(hit: Hit) -> tuple:
    """Executive-summary ranking: % identity over query length, then longer
    alignment, then lower E-value, then earlier subject start. Deliberately
    differs from the main ranking (no bit score; subject start as final tiebreak)."""
    return (hit.identity_over_query, hit.align_len, -hit.evalue, -hit.sstart)


def _full_query_alignment(hit: Hit, width: int = _EXEC_ALN_WIDTH) -> list[str]:
    """A pairwise alignment drawn against the COMPLETE query sequence: residues of
    the construct outside the matched region (the fusion's non-matching flanks) are
    still printed on the Query line with a blank subject beneath, so the reader sees
    how much of the whole construct the match covers -- not just the HSP."""
    full_q = hit.query_seq or hit.qseq.replace("-", "")
    left = full_q[: hit.qstart - 1]
    right = full_q[hit.qend :]
    bars = "".join(
        "|" if q == s and q != "-" else " " for q, s in zip(hit.qseq, hit.sseq, strict=True)
    )
    q_row = left + hit.qseq + right
    m_row = " " * len(left) + bars + " " * len(right)
    s_row = " " * len(left) + hit.sseq + " " * len(right)

    lines: list[str] = []
    qpos, spos = 1, hit.sstart
    for i in range(0, len(q_row), width):
        qseg, mseg, sseg = q_row[i : i + width], m_row[i : i + width], s_row[i : i + width]
        q_res = sum(1 for c in qseg if c != "-")  # q_row has no spaces
        s_res = sum(1 for c in sseg if c not in "- ")  # subject residues (flank = space)
        q_end = qpos + q_res - 1 if q_res else qpos
        s_lo = f"{spos:>6}" if s_res else " " * 6
        s_hi = f"  {spos + s_res - 1}" if s_res else ""
        lines.append(f"Query  {qpos:>6}  {qseg}  {q_end}".rstrip())
        lines.append(f"       {'':>6}  {mseg}".rstrip())
        lines.append(f"Sbjct  {s_lo}  {sseg}{s_hi}".rstrip())
        lines.append("")
        qpos += q_res
        spos += s_res
    while lines and lines[-1] == "":
        lines.pop()
    return lines


def _exec_best_match(hit: Hit) -> list[str]:
    """The headline block for one (construct, species): the strongest match with a
    full-query alignment."""
    label = parse_junction_label(hit.query_id)
    motif = label.motif if label else "?"
    span = f"{label.start}-{label.end}" if label else "?"
    block = [
        "",
        "    STRONGEST MATCH",
        f"      Protein : {_truncate(hit.target_name, 68)}  ({hit.accession})",
        f"      Identity: {hit.identity_over_query * 100:.0f}% of the query is identical to "
        f"this protein ({hit.nident} of {hit.query_len} residues)",
        f"      Junction: {motif}, residues {span} of the source sequence",
        f"      Details : alignment length {hit.align_len}, E-value {hit.evalue:.1e}, "
        f"query {hit.qstart}-{hit.qend}, subject {hit.sstart}-{hit.send}",
        "",
    ]
    block += ["      " + line for line in _full_query_alignment(hit)]
    return block


def _exec_other_targets(reps: list[Hit]) -> list[str]:
    """The compact ranked list of up to five distinct proteins the construct
    resembles. The E-value column lets the reader tell biologically meaningful hits
    (low E-value) from chance matches (E-value near 1). Row 1 is the headline,
    flagged with ``*`` (shown in full above)."""
    lines = [
        "",
        "    TOP DISTINCT PROTEINS THIS CONSTRUCT RESEMBLES  (ranked by % of query identical)",
        f"      {'%qry':>5}  {'E-value':>8}  {'Protein':<40}  {'Accession':<10}  Junction",
    ]
    for rank, hit in enumerate(reps, start=1):
        label = parse_junction_label(hit.query_id)
        motif = label.motif if label else "?"
        pct = f"{hit.identity_over_query * 100:.1f}"
        evalue = f"{hit.evalue:.1e}"
        name = _truncate(hit.target_name, 40)
        mark = " *" if rank == 1 else ""
        lines.append(f"      {pct:>5}  {evalue:>8}  {name:<40}  {hit.accession:<10}  {motif}{mark}")
    lines.append("      * strongest match, shown in full above")
    return lines


def _exec_construct_block(
    source: str,
    seq: str,
    taxa: list[Taxon],
    by_seq_species: dict[tuple[str, str, int], list[Hit]],
) -> list[str]:
    """The full block for one construct: a per-species headline + top-5 list."""
    rule = "=" * _EXEC_WIDTH
    lines = ["", rule, f"CONSTRUCT  {seq}", f"source file: {source}", rule]
    for taxon in taxa:
        lines += ["", _fit_species("  vs ", taxon, gap="  ")]
        group = by_seq_species.get((source, seq, taxon.taxid), [])
        if not group:
            lines.append("     (no hits in this species)")
            continue
        ranked = sorted(group, key=_exec_sort_key, reverse=True)
        reps: list[Hit] = []
        seen: set[str] = set()
        for hit in ranked:  # collapse by target name -> up to 5 distinct proteins
            if hit.target_name in seen:
                continue
            seen.add(hit.target_name)
            reps.append(hit)
            if len(reps) == 5:
                break
        lines += _exec_best_match(reps[0])
        lines += _exec_other_targets(reps)
    return lines


def _exec_header(params: SearchParams, input_path: Path, taxa: list[Taxon]) -> list[str]:
    backend, matrix_label, _, _ = _scoring_summary(params)
    lines = [
        "FRANKENSEARCH — EXECUTIVE SUMMARY",
        "",
        'What this is: each input construct ("seq") is an artificial fusion peptide,',
        "cut into overlapping junction windows (motifs). For every species searched,",
        "this report shows the one protein each construct most resembles by chance",
        "(its strongest match, with a full-length alignment), then up to five distinct",
        'proteins it resembles. Matches are ranked by "% of query identical" — the',
        "fraction of the construct's residues identical to the target protein. This is",
        "a chance-similarity screen, NOT evidence of homology or biological function.",
        "",
        f"Input:      {input_path}",
        f"Generated:  {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Backend:    {backend}",
    ]
    for i, taxon in enumerate(taxa):  # compact species list (one per line)
        lines.append(_fit_species("Species:    " if i == 0 else " " * 12, taxon))
    lines += [
        f"Matrix:     {matrix_label}",
        "Ranked by:  % of query identical (ties: longer alignment, then lower E-value,",
        "            then earlier subject position)",
    ]
    return lines


def write_exec_summary(
    path: Path,
    hits: list[Hit],
    *,
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
) -> None:
    """Write the plain-language executive summary: per construct (seq) x species, the
    single strongest chance-similarity match (with a full-query alignment) plus up to
    five distinct target proteins it resembles. Requires junction-style query ids
    (see ``parse_junction_label``); callers gate on ``labels_are_junctions``."""
    by_seq_species: dict[tuple[str, str, int], list[Hit]] = {}
    for hit in hits:
        label = parse_junction_label(hit.query_id)
        if label is None:
            continue  # non-junction ids are ignored (caller should have gated already)
        by_seq_species.setdefault((label.source, label.seq, hit.taxid), []).append(hit)

    # Construct order = first appearance in the query list; keep only constructs that
    # produced at least one hit in at least one species.
    seq_order: list[tuple[str, str]] = []
    seen_seq: set[tuple[str, str]] = set()
    for query in queries:
        label = parse_junction_label(query.id)
        if label is None:
            continue
        key = (label.source, label.seq)
        if key in seen_seq:
            continue
        seen_seq.add(key)
        if any((label.source, label.seq, t.taxid) in by_seq_species for t in taxa):
            seq_order.append(key)

    lines = _exec_header(params, input_path, taxa)
    for source, seq in seq_order:
        lines += _exec_construct_block(source, seq, taxa, by_seq_species)
    if not seq_order:
        lines += ["", "(no constructs had any hits to summarise)"]
    else:
        lines += ["", "=" * _EXEC_WIDTH, f"{len(seq_order)} construct(s) summarised."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
