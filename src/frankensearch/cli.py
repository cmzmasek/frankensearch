"""Command-line interface for FRANKENSEARCH."""

from __future__ import annotations

import os
import platform
import shlex
import sys
from enum import Enum
from pathlib import Path

import typer
from rich.table import Table

from . import __version__, blast, database, inputs, paths, scoring, uniprot
from .errors import FrankensearchError, UserError
from .inputs import Query
from .output import (
    HIT_COLUMNS,
    alignments_output_path,
    cell_text,
    filtered_output_paths,
    hit_column_header,
    output_paths,
    top1_output_paths,
    write_outputs,
)
from .runtime import STATE
from .search import Hit, SearchParams, run_search
from .taxonomy import Taxon, TaxonomyResolver, TaxonomyUnavailable
from .ui import console, error_panel, info, warn

app = typer.Typer(
    name="frankensearch",
    help=(
        "Search short fusion ('franken') protein sequences against "
        "species-restricted protein databases.\n\n"
        "This is [bold]not[/] a homology search: results are ranked by identity, "
        "not filtered by E-value — the goal is to find proteins that are similar "
        "by chance, not homologs.\n\n"
        "Typical workflow:\n"
        "  1. [cyan]frankensearch doctor[/]              check your setup\n"
        "  2. [cyan]frankensearch setup --taxids ...[/]  build local databases (once per species)\n"
        "  3. [cyan]frankensearch search ...[/]          run a search"
    ),
    no_args_is_help=True,
    add_completion=False,
    rich_markup_mode="rich",
)


class Matrix(str, Enum):
    identity = "identity"
    pam30 = "pam30"
    blosum45 = "blosum45"
    blosum62 = "blosum62"


class RankBy(str, Enum):
    identity_alignment = "identity-alignment"
    identity_query = "identity-query"
    alignment_length = "alignment-length"


class ProteomeSet(str, Enum):
    reference = "reference"
    swissprot = "swissprot"
    all = "all"


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"frankensearch {__version__}")
        raise typer.Exit()


def _parse_taxids(raw: str) -> list[int]:
    """Parse a comma/space separated taxid string into a list of ints."""
    taxids: list[int] = []
    for piece in raw.replace(",", " ").split():
        try:
            taxids.append(int(piece))
        except ValueError:
            raise UserError(
                f"'{piece}' is not a valid NCBI taxonomy ID.",
                hint="Provide numeric, comma-separated taxids, e.g. --taxids 9606,10029",
            ) from None
    if not taxids:
        raise UserError(
            "No taxonomy IDs were provided.",
            hint="Add --taxids with one or more species taxids, e.g. --taxids 9606",
        )
    return taxids


@app.callback()
def _main(
    version: bool = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the version and exit.",
    ),
    debug: bool = typer.Option(
        False,
        "--debug",
        help="Show full Python tracebacks on error (for bug reports).",
    ),
) -> None:
    STATE.debug = debug


@app.command()
def search(
    input_file: Path = typer.Argument(
        ..., metavar="INPUT", help="Query sequences: FASTA, TSV, or CSV (format auto-detected)."
    ),
    taxids: str = typer.Option(
        ..., "--taxids", "-t", help="Comma-separated species taxids to search, e.g. 9606,10029."
    ),
    num_hits: int = typer.Option(
        10, "-n", "--num-hits", min=1, help="Top hits to report per query, per species."
    ),
    rank_by: RankBy = typer.Option(
        RankBy.identity_alignment,
        "--rank-by",
        help=(
            "How to rank hits within each (query, species) group: "
            "'identity-alignment' (identical residues ÷ alignment length, default), "
            "'identity-query' (identical residues ÷ query length), or "
            "'alignment-length' (longest alignment first). Both identity ratios and the "
            "alignment length are always reported regardless."
        ),
    ),
    filter_by: float = typer.Option(
        None,
        "--filter-by",
        help="Threshold on the --rank-by metric. If set, also writes "
        "<prefix>_filtered_by_<x>.tsv/.txt with only hits at or above it; every "
        "query/species with none is listed as 'no hit above threshold'. A fraction "
        "0-1 for the identity modes, or a residue count for alignment-length.",
    ),
    matrix: Matrix = typer.Option(
        Matrix.identity,
        "--matrix",
        help="Scoring matrix. 'identity' is BLAST's built-in pure-identity matrix.",
    ),
    no_alignments: bool = typer.Option(
        False,
        "--no-alignments",
        help="Skip all pairwise alignment output: no _alignments.txt, and the "
        "_top1/_filtered .txt become table-only. The compact match column is kept.",
    ),
    output_query: bool = typer.Option(
        False,
        "--output-query",
        help="Also report the full query sequence (in addition to its name): a "
        "query_sequence column in the .tsv files and a QUERY SEQUENCES section in "
        "the .txt files.",
    ),
    ungapped: bool = typer.Option(False, "--ungapped", help="Only ungapped alignments."),
    seg: bool = typer.Option(
        False, "--seg/--no-seg", help="Low-complexity (SEG) filtering (off by default)."
    ),
    evalue: float = typer.Option(
        200000.0,
        "--evalue",
        help="E-value cap; intentionally high so hits are NOT filtered by E-value.",
    ),
    word_size: int = typer.Option(2, "--word-size", min=2, help="BLAST word size."),
    max_target_seqs: int = typer.Option(
        5000,
        "--max-target-seqs",
        min=1,
        help="Hits to retrieve from BLAST before local re-ranking.",
    ),
    remote: bool = typer.Option(
        False,
        "--remote",
        help="Search NCBI remotely instead of local DBs (falls back to a built-in matrix).",
    ),
    output: Path = typer.Option(
        None,
        "-o",
        "--output",
        help="Output prefix; writes <prefix>.tsv/.txt/.summary.md plus best-hit-only "
        "<prefix>_top1.tsv/.txt. Defaults to the input name.",
    ),
    out_dir: Path = typer.Option(
        None,
        "--out-dir",
        help="Directory to write outputs into (basename comes from -o or the input file).",
    ),
    force: bool = typer.Option(
        False, "--force", help="Overwrite existing output files instead of stopping."
    ),
    db_dir: Path = typer.Option(
        None, "--db-dir", help="Directory containing local BLAST databases."
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show the plan without running anything."
    ),
) -> None:
    """Search query sequences against one or more species databases."""
    taxid_list = _parse_taxids(taxids)
    if filter_by is not None:
        if rank_by is RankBy.alignment_length:
            if filter_by < 0:
                raise UserError(
                    "--filter-by must be 0 or greater for --rank-by alignment-length.",
                    hint="It is a residue count, e.g. --filter-by 12.",
                )
        elif not 0.0 <= filter_by <= 1.0:
            raise UserError(
                "--filter-by must be between 0 and 1 for identity ranking.",
                hint="It is a fraction, e.g. --filter-by 0.8 for 80% identity.",
            )
    base = output if output is not None else Path(input_file.stem)
    out_prefix = (out_dir / base.name) if out_dir is not None else base

    parsed = inputs.parse(input_file)
    for note in parsed.notes:
        info(f"[dim]{note}[/]")
    for message in parsed.warnings:
        warn(message)

    targets = _resolve_species(taxid_list, db_dir=db_dir, prefer_metadata=not remote)

    matrix_desc = matrix.value
    if matrix is Matrix.identity:
        matrix_desc += " (built-in identity scoring)"

    table = Table(title="Search plan", show_header=False, title_style="bold cyan")
    table.add_column(style="bold")
    table.add_column()
    table.add_row("Input", f"{input_file} (format: {parsed.fmt.value})")
    table.add_row("Queries", str(len(parsed.records)))
    table.add_row("Taxids", ", ".join(str(t) for t in taxid_list))
    table.add_row("Top hits per query/species", str(num_hits))
    table.add_row("Rank by", rank_by.value)
    if filter_by is not None:
        table.add_row("Filter", f"keep {rank_by.value} >= {filter_by:g}")
    if no_alignments:
        table.add_row("Alignments", "off (--no-alignments)")
    if output_query:
        table.add_row("Query sequence", "included in outputs")
    table.add_row("Matrix", matrix_desc)
    table.add_row("Gaps", scoring.gap_description(matrix.value, ungapped=ungapped))
    table.add_row("SEG filter", "on" if seg else "off")
    table.add_row("E-value cap", f"{evalue:g} (not used for filtering)")
    table.add_row("Word size", str(word_size))
    table.add_row("Max target seqs", str(max_target_seqs))
    table.add_row("Backend", "NCBI remote" if remote else "local databases")
    tsv_path, txt_path, summary_path = output_paths(out_prefix)
    alignments_path = alignments_output_path(out_prefix)
    top1_tsv_path, top1_txt_path = top1_output_paths(out_prefix)
    all_output_paths = [tsv_path, txt_path]
    if not no_alignments:
        all_output_paths.append(alignments_path)
    all_output_paths += [summary_path, top1_tsv_path, top1_txt_path]
    if filter_by is not None:
        all_output_paths += list(filtered_output_paths(out_prefix, filter_by))
    all_output_paths = tuple(all_output_paths)
    table.add_row("Outputs", "\n".join(str(p) for p in all_output_paths))
    console.print(table)
    if targets:
        _print_targets(targets)

    if dry_run:
        _print_query_preview(parsed.records)
        console.print("\n[green]Dry run:[/] no search was executed.")
        return

    if not targets:
        raise UserError(
            "No taxids could be resolved, so there is nothing to search.",
            hint="Check your internet connection and the taxids.",
        )

    _ensure_writable_prefix(out_prefix)
    _check_no_overwrite(all_output_paths, force)

    if remote:
        _, matrix_warning = scoring.remote_blast_args(matrix.value, ungapped=ungapped)
        if matrix_warning:
            warn(matrix_warning)

    params = SearchParams(
        matrix=matrix.value,
        ungapped=ungapped,
        evalue=evalue,
        word_size=word_size,
        max_target_seqs=max_target_seqs,
        rank_by=rank_by.value,
        num_hits=num_hits,
        remote=remote,
        seg=seg,
        filter_by=filter_by,
        include_alignments=not no_alignments,
        output_query=output_query,
    )
    status = (
        "Searching NCBI remotely (this can take a while)…"
        if remote
        else "Running BLAST searches…"
    )
    search_warnings: list[str] = []
    with console.status(status):
        results = run_search(
            parsed.records, targets, params, db_dir=db_dir, on_warning=search_warnings.append
        )
    hits = results.hits
    for message in search_warnings:
        warn(message)

    db_metadata = {}
    if not remote:
        for taxon in targets:
            meta = database.load_metadata(taxon.taxid, db_dir)
            if meta is not None:
                db_metadata[taxon.taxid] = meta
    command = "frankensearch " + " ".join(shlex.quote(arg) for arg in sys.argv[1:])

    _print_results(hits, rank_by.value)
    write_outputs(
        hits,
        out_prefix,
        queries=parsed.records,
        taxa=targets,
        params=params,
        input_path=input_file,
        command=command,
        frankensearch_version=__version__,
        blast_versions={
            "blastp": blast.tool_version("blastp"),
            "makeblastdb": blast.tool_version("makeblastdb"),
        },
        db_metadata=db_metadata,
        top1_hits=results.top1,
    )
    files_block = "\n".join(f"  {p}" for p in all_output_paths)
    console.print(
        f"\n[green]Done.[/] {len(hits)} hit(s) written to:\n{files_block}\n"
        f"  [dim]({len(results.top1)} best-hit row(s) in the _top1 files, ties included)[/]"
    )


def _resolve_species(
    taxid_list: list[int], *, db_dir: Path | None = None, prefer_metadata: bool = False
) -> list[Taxon]:
    """Resolve taxids to species, warning (not failing) when a name can't be fetched.

    When ``prefer_metadata`` is set (local search), name/rank are read from an
    existing database's metadata.json first — no network needed if the DB is built.
    UniProt is consulted only for taxids without local metadata. A genuinely wrong
    taxid (TaxonNotFound) still propagates as a clean fatal error; only transient
    unavailability (offline / not cached) degrades to a warning.
    """
    taxa: list[Taxon] = []
    with TaxonomyResolver() as resolver:
        for taxid in taxid_list:
            taxon: Taxon | None = None
            if prefer_metadata:
                meta = database.load_metadata(taxid, db_dir)
                if meta is not None:
                    taxon = Taxon(taxid=taxid, name=meta.scientific_name, rank=meta.rank)
            if taxon is None:
                try:
                    taxon = resolver.resolve(taxid)
                except TaxonomyUnavailable as exc:
                    warn(str(exc))
                    continue
            taxa.append(taxon)
            if not taxon.is_species:
                warn(
                    f"taxid {taxon.taxid} ({taxon.name}) has rank '{taxon.rank}', not species — "
                    "species-level taxids are recommended so one taxon doesn't dominate."
                )
    return taxa


def _print_targets(taxa: list[Taxon]) -> None:
    table = Table(title="Target species", title_style="bold cyan")
    table.add_column("Taxid", justify="right")
    table.add_column("Scientific name")
    table.add_column("Rank")
    for taxon in taxa:
        rank = taxon.rank if taxon.is_species else f"[yellow]{taxon.rank}[/]"
        table.add_row(str(taxon.taxid), taxon.name, rank)
    console.print(table)


def _print_query_preview(records: list[Query], limit: int = 25) -> None:
    table = Table(title="Parsed queries", title_style="bold cyan")
    table.add_column("#", justify="right")
    table.add_column("Query ID")
    table.add_column("Length", justify="right")
    for index, query in enumerate(records[:limit], start=1):
        table.add_row(str(index), query.id, str(len(query.sequence)))
    console.print(table)
    if len(records) > limit:
        console.print(f"... and {len(records) - limit} more")


def _ensure_writable_prefix(out_prefix: Path) -> None:
    """Fail fast (before a possibly long search) if outputs can't be written."""
    parent = out_prefix.parent if str(out_prefix.parent) else Path(".")
    try:
        parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise UserError(
            f"Cannot create the output directory {parent}: {exc}",
            hint="Pick a different location with -o/--output.",
        ) from exc
    if not os.access(parent, os.W_OK):
        raise UserError(
            f"The output directory is not writable: {parent}",
            hint="Pick a writable location with -o/--output.",
        )


def _check_no_overwrite(paths: tuple[Path, ...], force: bool) -> None:
    """Refuse to clobber existing output files unless --force was given."""
    if force:
        return
    existing = [p for p in paths if p.exists()]
    if existing:
        listing = ", ".join(str(p) for p in existing)
        raise UserError(
            f"Output file(s) already exist: {listing}.",
            hint="Re-run with --force to overwrite, or change -o/--output or --out-dir.",
        )


def _print_results(hits: list[Hit], rank_by: str, limit: int = 20) -> None:
    if not hits:
        console.print("\n[yellow]No hits found.[/]")
        return
    table = Table(title=f"Top hits (ranked by {rank_by})", title_style="bold cyan")
    columns = [col for col in HIT_COLUMNS if col.in_console]
    for col in columns:
        header = hit_column_header(col, rank_by)
        if col.max_width is not None:  # long text: truncate to one line, like the .txt table
            table.add_column(header, no_wrap=True, max_width=col.max_width)
        else:
            table.add_column(header, justify="right" if col.numeric else "left")
    for hit in hits[:limit]:
        table.add_row(*(cell_text(col, hit) for col in columns))
    console.print(table)
    if len(hits) > limit:
        console.print(f"[dim]... and {len(hits) - limit} more (see the .txt / .tsv files).[/]")


@app.command()
def setup(
    taxids: str = typer.Option(
        ..., "--taxids", "-t", help="Species taxids to build databases for, e.g. 9606,10029."
    ),
    proteome_set: ProteomeSet = typer.Option(
        ProteomeSet.reference,
        "--proteome-set",
        help=(
            "UniProt set to download per species. "
            "'reference' (default): the species' reference proteome — one protein "
            "per gene, mixing reviewed (Swiss-Prot) and unreviewed (TrEMBL) entries; "
            "the recommended, complete-but-non-redundant search space. "
            "'swissprot': reviewed entries only — small and high quality, but may be "
            "sparse or empty for non-model organisms. "
            "'all': every UniProtKB entry — largest and most redundant (isoforms, "
            "fragments, strains)."
        ),
    ),
    db_dir: Path = typer.Option(None, "--db-dir", help="Where to store the built databases."),
    force: bool = typer.Option(False, "--force", help="Rebuild even if a database already exists."),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show what would be built without downloading."
    ),
) -> None:
    """Download UniProt proteomes and build local BLAST databases (run once per species)."""
    taxid_list = _parse_taxids(taxids)
    targets = _resolve_species(taxid_list)
    if not targets:
        raise UserError(
            "None of the requested taxids could be resolved.",
            hint="Check your internet connection and the taxids, then try again.",
        )

    resolved_db_dir = db_dir or paths.database_dir()

    table = Table(title="Setup plan", show_header=False, title_style="bold cyan")
    table.add_column(style="bold")
    table.add_column()
    table.add_row("Taxids", ", ".join(str(t) for t in taxid_list))
    table.add_row("Proteome set", proteome_set.value)
    table.add_row("Database dir", str(resolved_db_dir))
    table.add_row("Force rebuild", "yes" if force else "no")
    console.print(table)
    _print_targets(targets)

    if dry_run:
        for taxon in targets:
            already = (not force) and database.is_current(taxon.taxid, proteome_set.value, db_dir)
            verb = "already built — would skip" if already else "would download + build"
            info(f"  • taxid {taxon.taxid} ({taxon.name}): {verb}")
        console.print("\n[green]Dry run:[/] nothing was downloaded or built.")
        return

    built, skipped = 0, 0
    with uniprot.make_client() as client:
        for taxon in targets:
            if not force and database.is_current(taxon.taxid, proteome_set.value, db_dir):
                info(
                    f"taxid {taxon.taxid} ({taxon.name}) already built "
                    f"[{proteome_set.value}] — skipping (use --force to rebuild)."
                )
                skipped += 1
                continue
            with console.status(f"Building database for {taxon.name} (taxid {taxon.taxid})…"):
                meta = database.build(
                    taxon, proteome_set.value, db_dir=db_dir, force=force, client=client
                )
            info(
                f"[green]✓[/] taxid {taxon.taxid} ({taxon.name}): "
                f"{meta.sequence_count:,} sequences."
            )
            built += 1

    console.print(
        f"\n[green]Setup complete.[/] Built {built}, skipped {skipped}. "
        f"Databases in {resolved_db_dir}"
    )


@app.command()
def doctor(
    taxids: str = typer.Option(
        None, "--taxids", "-t", help="Also check that local databases exist for these taxids."
    ),
    db_dir: Path = typer.Option(
        None, "--db-dir", help="Directory containing local BLAST databases."
    ),
) -> None:
    """Check that everything FRANKENSEARCH needs is installed and working."""
    table = Table(title="frankensearch doctor", title_style="bold cyan")
    table.add_column("Check")
    table.add_column("Status")
    table.add_column("Detail", overflow="fold")

    all_ok = True
    for tool in blast.REQUIRED_TOOLS:
        version = blast.tool_version(tool)
        if version:
            table.add_row(tool, "[green]OK[/]", version)
        else:
            table.add_row(tool, "[red]missing[/]", "not found on PATH")
            all_ok = False

    matrix_ok, matrix_detail = blast.selftest_matrix("identity")
    table.add_row(
        "identity matrix", "[green]OK[/]" if matrix_ok else "[red]error[/]", matrix_detail
    )
    if not matrix_ok and blast.find_tool("blastp") is not None:
        all_ok = False

    table.add_row("python", "[green]OK[/]", platform.python_version())
    table.add_row("frankensearch", "[green]OK[/]", __version__)

    cache_probe = TaxonomyResolver()
    cache_status = "[green]OK[/]" if cache_probe.cache_size else "[yellow]empty[/]"
    table.add_row(
        "taxonomy cache",
        cache_status,
        f"{cache_probe.cache_size} taxid(s) cached at {cache_probe.cache_file}",
    )

    built_dbs = database.list_built(db_dir)
    table.add_row(
        "databases",
        "[green]OK[/]" if built_dbs else "[yellow]none[/]",
        f"{len(built_dbs)} built in {db_dir or paths.database_dir()}",
    )

    if taxids:
        with TaxonomyResolver() as resolver:
            for taxid in _parse_taxids(taxids):
                try:
                    taxon = resolver.resolve(taxid)
                except FrankensearchError as exc:
                    table.add_row(f"taxid {taxid}", "[red]error[/]", str(exc))
                    continue
                status = "[green]OK[/]" if taxon.is_species else "[yellow]not species[/]"
                db_note = (
                    "database built"
                    if database.is_built(taxid, db_dir)
                    else "[yellow]no database (run setup)[/]"
                )
                table.add_row(
                    f"taxid {taxid}", status, f"{taxon.name} (rank: {taxon.rank}); {db_note}"
                )

    console.print(table)

    if not all_ok:
        console.print()
        error_panel(
            "Setup incomplete",
            "BLAST+ tools were not found on your PATH.",
            hint="conda activate frankensearch   (or install: conda install -c bioconda blast)",
        )
        raise typer.Exit(1)

    console.print("\n[green]All required tools are available.[/]")


@app.command()
def databases(
    db_dir: Path = typer.Option(
        None, "--db-dir", help="Directory containing local BLAST databases."
    ),
) -> None:
    """List the local species databases that have been built."""
    root = db_dir or paths.database_dir()
    metas = database.list_built(db_dir)
    if not metas:
        console.print(f"No databases have been built yet in {root}.")
        console.print("Build one with:  [cyan]frankensearch setup --taxids 9606[/]")
        return

    table = Table(title="Local databases", title_style="bold cyan")
    table.add_column("Taxid", justify="right")
    table.add_column("Species")
    table.add_column("Set")
    table.add_column("Sequences", justify="right")
    table.add_column("UniProt release")
    table.add_column("Built")
    for meta in metas:
        table.add_row(
            str(meta.taxid),
            meta.scientific_name,
            meta.proteome_set,
            f"{meta.sequence_count:,}",
            meta.uniprot_release or "?",
            meta.built_at,
        )
    console.print(table)
    console.print(f"[dim]Location: {root}[/]")


def main() -> None:
    """Entry point: run the CLI with friendly error handling."""
    try:
        app()
    except FrankensearchError as exc:
        error_panel("Error", str(exc), getattr(exc, "hint", None))
        raise SystemExit(2) from None
    except Exception as exc:  # noqa: BLE001 -- last-resort guard for unexpected bugs
        if STATE.debug:
            raise
        error_panel(
            "Unexpected error",
            f"{type(exc).__name__}: {exc}",
            hint="This is probably a bug. Re-run with --debug to see the full traceback.",
        )
        raise SystemExit(1) from None
