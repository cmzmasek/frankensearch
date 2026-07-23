# FRANKENSEARCH

Chance-similarity search of short fusion ("franken") protein sequences against
species-restricted protein databases.

Unlike an ordinary BLAST search, FRANKENSEARCH is **not** looking for homologs.
Its inputs are artificial fusion proteins, so the usual E-value/homology
statistics are the wrong lens. Results are therefore **ranked by identity and are
*not* filtered by E-value** — the goal is to surface proteins that are similar by
chance, including short, high-identity matches that homology search would discard.

---

## Installation

FRANKENSEARCH drives the NCBI **BLAST+** tools (`blastp`, `makeblastdb`), so those
must be installed. The easiest way is a dedicated conda environment:

```sh
conda create -n frankensearch -c conda-forge -c bioconda python=3.11 blast
conda activate frankensearch
pip install -e ".[dev]"
```

Check everything is ready:

```sh
frankensearch doctor
```

`doctor` verifies that BLAST+ is found, that the identity matrix loads, and reports
the taxonomy cache and which databases you have built.

---

## Quick start

```sh
# 1. Build a local database for each species you want to search (once per species)
frankensearch setup --taxids 9606            # human (UniProt reference proteome)

# 2. Search your franken proteins against it
frankensearch search myproteins.fasta --taxids 9606 -n 10

# 3. Read the results
#    myproteins.txt             -> human-readable hit table
#    myproteins.tsv             -> for downstream analysis
#    myproteins_alignments.txt  -> the full pairwise alignments
#    myproteins_top1.txt/.tsv   -> just the best hit per query/species
```

There is an example input at `examples/franken_demo.fasta`.

---

## Inputs

**Query sequences** — short amino-acid sequences (typically < 200 residues). The
format is auto-detected:

| Format | Shape |
|--------|-------|
| FASTA  | `>name` lines followed by sequence |
| TSV    | two columns: name `<tab>` sequence |
| CSV    | two columns: name `,` sequence |

Bad records (invalid characters, empty sequences, etc.) are skipped with a plain
warning rather than aborting the run.

**Taxids** — one or more NCBI Taxonomy IDs (`--taxids 9606,10090`). They should be
**species-level** (e.g. human, mouse), not higher clades; FRANKENSEARCH warns if a
taxid is not species-rank. Each taxid is searched separately so no single taxon
drowns out another.

---

## Key options (`frankensearch search --help`)

| Option | Meaning |
|--------|---------|
| `-n, --num-hits` | Top hits to report **per (query, species)** (default 10). A warning + a `.summary.md` note flag any group that had more hits than the cap. |
| `--rank-by {identity-query,identity-alignment,alignment-length}` | How to **rank** hits within each (query, species) group: identity ÷ query length (**default**), identity ÷ alignment length, or longest alignment first. Both ratios and the alignment length are always reported. (Default is `identity-query` because, within a query, it orders by identical-residue count — i.e. by how surprising the match is under the chance model.) |
| `--filter-by VALUE` | Threshold on the `--rank-by` metric. Writes two extra files keeping only hits at or above it, and lists every query/species with **no** passing hit. A fraction `0–1` for the identity modes, or a residue count for `alignment-length`. See [Output](#output). |
| `--matrix {identity,pam30,blosum45,blosum62}` | Scoring matrix; `identity` (built-in pure-identity) is the default. |
| `--no-alignments` | Skip all pairwise alignment output (no `_alignments.txt`; `_top1`/`_filtered` `.txt` become table-only). The compact `match` column is kept. |
| `--output-query` | Also report the full query sequence (in addition to its name): a `query_sequence` column in the `.tsv` files and a `Query-Seq` column in the `.txt` HITS tables. |
| `--exec-summary` | Also write `<prefix>_executive_summary.txt` (a plain-language, per-construct roll-up — each construct's strongest match plus up to five distinct proteins it resembles, per species, with a full-query alignment) and `<prefix>_executive_summary_matrix.tsv` (a construct × species matrix of highest % of query identical). Needs junction-style query IDs (as produced by `extract_junctions.py`); skipped with a warning otherwise. See [Executive summary](#executive-summary). |
| `--ungapped` | Ungapped alignments only. |
| `--remote` | Search NCBI remotely instead of using local databases (see below). |
| `-o, --output` | Output path prefix (defaults to the input file's name). |
| `--dry-run` | Show the plan (parsed queries, resolved species) without searching. |
| `--debug` | Show full tracebacks (otherwise errors are concise, friendly messages). |

> **`-n` only limits the report, not the search.** BLAST retrieves up to
> `--max-target-seqs` hits (default 5000) regardless of `-n`; those are ranked
> locally, and `-n` just keeps the top of that ranked list. So the **best hit is
> always kept for any `n ≥ 1`** — a small `-n` only drops lower-ranked extras (which
> the truncation warning flags). The single best hit per (query, species) is also
> always written to the `_top1` files, independent of `-n`. The knob that affects what
> BLAST *finds* is `--max-target-seqs`, not `-n`.

---

## Output

Seven files are written per run (nine with `--filter-by`, plus two with
`--exec-summary`), sharing the `-o/--output` prefix. The table outputs
(`.txt`/`.tsv`) are kept **compact** (one row per hit) so they stay manageable for
large, many-query runs; the bulky pairwise alignments go to their own files (skip
them all with `--no-alignments`).

- **`.txt`** — human-readable: a fixed-width table of every hit (query, species,
  target, **both** identity ratios, alignment length, bit score, E-value,
  query/subject start–end coordinates, and a final **match** column). The match
  column shows the query residue where the two sequences are identical and a dot
  (`.`) where they differ (e.g. `MKL.EV`). It points to the `_alignments.txt` file
  for the full alignments. With `--output-query` the HITS table gains a `Query-Seq`
  column (right after `Query`). The header reports, per database, a **chance-match
  length `k*`**: matches **longer** than `k*` stand out from the random-chance
  background, while shorter exact matches are *expected* for short queries against a
  whole proteome (see `.summary.md` for the formula). This is interpretive context,
  not a significance filter.
- **`.tsv`** — one row per hit for downstream processing: query ID, queried taxid +
  species, target accession + name, **both** identity ratios (over alignment length
  and over query length), bit score, E-value, alignment coordinates, and a compact
  **`match`** column (same dotted string as the `.txt`). With `--output-query` a
  `query_sequence` column is added.
- **`_alignments.txt`** — the BLAST-style pairwise alignments for every hit, grouped
  by query then species (kept out of the `.txt` so it stays small). Ignore or delete
  it if you don't need the alignments, or pass `--no-alignments` to skip all
  alignment output (no `_alignments.txt`/`_top1_alignments.txt` are written).
- **`_top1.txt` / `_top1.tsv`** — the same table views filtered to just the **single
  best hit per (query, species)**, with the alignments split off to
  `_top1_alignments.txt`. More than one row appears only on a genuine tie (hits the
  ranking cannot separate at all). Independent of `-n`, so a tie is never dropped.
- **`_filtered_by_<x>.txt` / `_filtered_by_<x>.tsv`** — *only when `--filter-by <x>`
  is given.* Keeps just the hits whose `--rank-by` metric is ≥ `<x>`, and — the key
  point — **explicitly lists every query/species that has no hit above the
  threshold** (the `.txt` has a dedicated section; the `.tsv` adds a leading
  `status` column with `hit` / `no_hit_above_threshold`). Built so you can find
  queries that match *nothing* without post-filtering in Excel.
- **`.summary.md`** — a methods-grade record of the run (command, input checksum,
  per-species database provenance, full effective parameters, software versions,
  references, a **chance-reference** section with the `k*` formula + `p` + per-species
  `M`/`k*`, and a ready-to-paste Methods paragraph).
- **`_executive_summary.txt`** — *only with `--exec-summary`.* A plain-language
  roll-up for a non-specialist reader (see below).
- **`_executive_summary_matrix.tsv`** — *only with `--exec-summary`.* A
  construct × species matrix: one row per construct, one column per species, each
  cell the highest *% of query identical* that construct reached against that
  species (`0` when nothing hit). A spreadsheet-friendly pivot of the executive
  summary — ideal for spotting a construct that matches its positive-control
  species but nothing else.

> E-value is reported for reference only; it is never used to filter results.

### Executive summary

`--exec-summary` writes `<prefix>_executive_summary.txt`, a self-contained,
plain-language report. It opens with an **overview table** — per species, how many
constructs have at least one junction reaching each *% of query identical* level
(`>=50 / >=70 / >=90 / =100%`), out of all constructs searched. Then, for **each
construct (`seq`) in each species**, it shows:

- the **strongest match** — the single protein that construct most resembles by
  chance (ranked by *% of query identical*, then longer alignment, lower E-value,
  earlier subject position), drawn as a **full-query alignment** so the flanks of
  the fusion that *don't* match are visible, not just the aligned core; and
- **up to five distinct proteins** it resembles, collapsed by target so one protein
  isn't listed five times, each with its **E-value** — so a biologically meaningful
  hit (low E-value) is easy to tell from a chance match (E-value near 1).

With `--filter-by`, a second overview table counts constructs that clear the cutoff
(on the active `--rank-by` metric), and every shown match at or above it is flagged
(`>=` in the list, a tag on the strongest-match heading). The report is laid out at
a fixed 100-column width. It requires junction-style query
IDs of the form
`source|seq|motif_N|start_end` — exactly what `extract_junctions.py` produces. If a
run's IDs don't match, the file is skipped with a warning (the rest of the run is
unaffected).

Alongside it, `--exec-summary` also writes **`<prefix>_executive_summary_matrix.tsv`**
— the same information as a compact **construct × species matrix**: one row per
construct, one column per species, each cell the highest *% of query identical* any
of that construct's junctions reached against that species (`0` for no hit). Open it
in Excel to scan for the design goal at a glance — high in the intended
(positive-control) species, near-zero elsewhere.

---

## Matrix choice changes how junctions align

Because a query is an artificial fusion, its two halves often come from
different positions of the *same* target protein — farther apart there than in
the query. Reporting the fusion as **one** alignment then requires opening a gap
across the junction, and whether BLAST does that depends on the **gap cost**,
which is tied to `--matrix`:

- `--matrix identity` is forced to gap **open 15 / extend 2** (the only gapped
  costs BLAST+ allows for the IDENTITY matrix).
- the other matrices use blastp's cheaper per-matrix default (e.g. PAM30 =
  **open 9 / extend 1**).

Cheap gaps let BLAST **stitch** both flanks into a single alignment; expensive
gaps make it cheaper to report **two** separate HSPs, one per flank. Since
`identity-query` divides identical residues by query length *for the best single
HSP*, splitting the match roughly halves that ratio — so the same query can
clear a `--filter-by 0.7` threshold under PAM30 but not under `identity`:

```
query  GQVFGLYKNTCVGS GGSGGGGSGG CTERLKLFAAETLK   (two SARS fragments + GS linker)

identity  -> 2 HSPs:  14/38 = 0.37  and  14/38 = 0.37   (below 0.7)
pam30     -> 1 HSP:   29/38 = 0.76                        (above 0.7)
```

Neither is "more correct" — both find the same identical residues, just packaged
differently. If you do **not** want alignments bridging the synthetic linker
(cleaner for chance detection), keep `identity` or add `--ungapped`; if you
**do** want the whole fusion scored as one span, use a non-identity matrix. This
is also why `--remote`, which falls back to PAM30, can report more or longer
single-HSP hits than a local identity run.

---

## Local vs. remote

- **Local (default, recommended).** `setup` downloads each species' UniProt proteome
  and builds a BLAST database under `~/.frankensearch/blastdb/`. Reproducible,
  offline after setup, and uses the pure-identity matrix.
  - `--proteome-set {reference,swissprot,all}` chooses what to download:
    - `reference` (default) — the species' **reference proteome**: one protein per
      gene, mixing reviewed (Swiss-Prot) and unreviewed (TrEMBL) entries. Complete
      but non-redundant; the recommended search space.
    - `swissprot` — **reviewed entries only**. Small and high quality, but can be
      sparse or empty for non-model organisms.
    - `all` — **every UniProtKB entry** for the organism. Largest and most
      redundant (isoforms, fragments, strains).

    `reference` already includes most Swiss-Prot sequences (it uses the reviewed
    entry per gene where one exists), so it is not a strict superset of `swissprot`
    but overlaps it heavily.
- **Remote (`--remote`).** Searches NCBI's `nr` remotely (no local database needed),
  restricting to each taxid. Convenient for one-offs, but:
  - NCBI's remote service has no IDENTITY matrix, so it falls back to **PAM30**
    (with a warning).
  - `nr` is non-redundant, so a hit's *listed organism* may differ from the queried
    taxid (the output notes this).
  - It is slower and subject to NCBI's load.

See what you have built:

```sh
frankensearch databases
```

---

## Where data lives

Everything is stored under `~/.frankensearch/` (taxonomy cache + BLAST databases).
Override the location with the `FRANKENSEARCH_HOME` environment variable.

---

## Companion script: `extract_junctions.py`

`extract_junctions.py` (in the project root) is a small standalone helper for
preparing **junction peptides** to search. Given a motif (e.g. a linker), it
extracts every occurrence of that motif plus a fixed number of flanking residues
from your sequences and writes them as a FASTA — ready to feed to FRANKENSEARCH or
paste into NCBI blastp.

It reads a *folder* of input files (FASTA, or TSV/CSV with column 1 = name and
column 2 = sequence; a header row and any extra columns are ignored), chosen by
file extension, and writes one `<name>_<format>_extracted.fasta` per input file
(the format is part of the name so same-named inputs of different formats don't
overwrite each other).

```sh
python3 extract_junctions.py -i input_folder -o output_folder -m GGSGGGGSGG -f 14
```

| Option | Meaning |
|--------|---------|
| `-i, --input` | Folder of FASTA/TSV/CSV input files. |
| `-o, --output` | Output folder (created if needed). |
| `-m, --motif` | Motif to center each peptide on (matched literally; default `GGSGGGGSGG`). |
| `-f, --flank` | Flanking residues kept on each side of the motif (default 14). |

Each output record is one motif occurrence with its flanks, named
`>input_file|sequence_name|motif_N|start_end`. Rows whose sequence contains
non-amino-acid characters are skipped with a warning.

---

## Troubleshooting

- **"BLAST+ tools were not found"** — `conda activate frankensearch` (or install
  BLAST+). Run `frankensearch doctor` to confirm.
- **"No local database for ..."** — build it: `frankensearch setup --taxids <id>`.
- **A taxid won't resolve** — check it at <https://www.ncbi.nlm.nih.gov/taxonomy>;
  taxonomy lookups need internet on first use (results are then cached).
- For a full traceback when reporting a bug, re-run with `--debug`.

---

## License

FRANKENSEARCH is free software, licensed under the **GNU General Public License
v3.0 or later** (GPL-3.0-or-later). See the [`LICENSE`](LICENSE) file for the
full text.
