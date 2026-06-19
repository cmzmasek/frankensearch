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
#    myproteins.txt  -> human-readable, with alignments
#    myproteins.tsv  -> for downstream analysis
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
| `-n, --num-hits` | Top hits to report **per (query, species)** (default 10). |
| `--rank-by {identity-alignment,identity-query,alignment-length}` | How to **rank** hits within each (query, species) group: identity ÷ alignment length (default), identity ÷ query length, or longest alignment first. Both ratios and the alignment length are always reported. |
| `--matrix {identity,pam30,blosum45,blosum62}` | Scoring matrix; `identity` (built-in pure-identity) is the default. |
| `--ungapped` | Ungapped alignments only. |
| `--remote` | Search NCBI remotely instead of using local databases (see below). |
| `-o, --output` | Output path prefix (defaults to the input file's name). |
| `--dry-run` | Show the plan (parsed queries, resolved species) without searching. |
| `--debug` | Show full tracebacks (otherwise errors are concise, friendly messages). |

---

## Output

Three files are written per run, sharing the `-o/--output` prefix:

- **`.txt`** — human-readable, in two sections: (1) a table of every hit (query,
  species, target, **both** identity ratios, alignment length, bit score, E-value,
  query/subject start–end coordinates, and a final **match** column), then (2) the
  BLAST-style pairwise alignments grouped by query then species. The match column /
  middle "match" line shows the query residue where the two sequences are identical
  and a dot (`.`) where they differ (e.g. `MKL.EV`).
- **`.tsv`** — one row per hit for downstream processing, with columns including the
  query ID, queried taxid + species, target accession + name, **both** identity
  ratios (over alignment length and over query length), bit score, E-value,
  alignment coordinates, and the alignment itself (as a single field with newlines
  escaped as `\n`).
- **`.summary.md`** — a methods-grade record of the run (command, input checksum,
  per-species database provenance, full effective parameters, software versions,
  references, and a ready-to-paste Methods paragraph).

> E-value is reported for reference only; it is never used to filter results.

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
