"""Write search results to a machine-readable .tsv and a human-readable .txt.

Both files carry the alignment: the .tsv encodes it as a single field (newlines
escaped as ``\\n``), the .txt shows the multi-line BLAST-style view.
"""

from __future__ import annotations

import csv
import hashlib
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


def _scoring_summary(params: SearchParams) -> tuple[str, str, str, str]:
    """Return (backend, matrix_label, gaps, ranked_by) describing the effective run."""
    ranked_by = "query length" if params.identity_denominator == "query" else "alignment length"
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
    """The three output files for a run: (.tsv, .txt, .summary.md)."""
    base = out_prefix
    return (
        base.parent / f"{base.name}.tsv",
        base.parent / f"{base.name}.txt",
        base.parent / f"{base.name}.summary.md",
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
) -> tuple[Path, Path, Path]:
    """Write the .tsv, .txt, and summary.md; return their paths."""
    tsv_path, txt_path, summary_path = output_paths(out_prefix)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(hits, tsv_path)
    write_txt(hits, txt_path, queries=queries, taxa=taxa, params=params, input_path=input_path)
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
    return tsv_path, txt_path, summary_path


def write_tsv(hits: list[Hit], path: Path) -> None:
    with open(path, "w", newline="") as handle:
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


def write_txt(
    hits: list[Hit],
    path: Path,
    *,
    queries: list[Query],
    taxa: list[Taxon],
    params: SearchParams,
    input_path: Path,
) -> None:
    by_group: dict[tuple[str, int], list[Hit]] = {}
    for hit in hits:
        by_group.setdefault((hit.query_id, hit.taxid), []).append(hit)

    backend, matrix_label, gaps, ranked_by = _scoring_summary(params)
    lines: list[str] = [
        "FRANKENSEARCH results",
        f"Input:      {input_path}",
        f"Generated:  {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Backend:    {backend}",
        f"Matrix:     {matrix_label}    Gaps: {gaps}",
        f"Ranked by:  identity over {ranked_by}",
        f"Top hits per (query, species): {params.num_hits}",
        "Note:       E-value is shown for reference only; results are NOT filtered by E-value.",
    ]
    if params.remote:
        lines.append(
            "Note:       Remote nr is non-redundant; a hit's listed organism may "
            "differ from the queried taxid."
        )

    for query in queries:
        lines.append("")
        lines.append("=" * 80)
        lines.append(f"Query: {query.id}   (length {len(query.sequence)})")
        lines.append("=" * 80)
        for taxon in taxa:
            lines.append("")
            lines.append(f"--- {taxon.name}  (taxid {taxon.taxid}) ---")
            group = by_group.get((query.id, taxon.taxid), [])
            if not group:
                lines.append("    (no hits)")
                continue
            for rank, hit in enumerate(group, start=1):
                lines.extend(_hit_block(rank, hit))

    path.write_text("\n".join(lines) + "\n")


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
        f"- Ranking: by identity over {ranked_by} (both ratios are reported)",
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
            f"ranked by percent identity over {ranked_by}. The top {params.num_hits} hits per "
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
    path.write_text("\n".join(lines) + "\n")
