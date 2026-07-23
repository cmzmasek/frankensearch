#!/usr/bin/env python3
"""Validate FRANKENSEARCH's chance-match length k* against a random-query null.

WHAT k* CLAIMS
--------------
FRANKENSEARCH reports, per database, the exact-match length expected ~once by
chance between one query and the whole DB:

    k* = ln(M * Q) / ln(1 / p)          (rounded to an integer in the tool)

    M = database residues, Q = query length,
    p = P(two residues are identical) = Sum f_i^2 = 0.0598
        over the Robinson & Robinson (1991) amino-acid frequencies.

Derivation (the thing we are testing): between a length-Q query and a length-M
database there are ~Q*M ungapped alignment offsets. Under an i.i.d. model each
offset yields an exact match of length >= k with probability p^k, so the expected
number of chance matches of length >= k is

    lambda(k) = Q * M * p^k .

Setting lambda(k*) = 1 gives k* above. Two crisp, falsifiable predictions follow
for L = the *longest exact match* of a random query against the DB:

  (1)  P(L >= k)      ~=  1 - exp(-Q * M * p^k)                (full CDF)
  (2)  P(L >= k*)     ~=  1 - 1/e ~= 0.632   at the real-valued k*
       (k* is the "expected-count-1" point, i.e. the ~63rd percentile of L,
        NOT the mean -- a specific number the experiment can check)

Across databases of very different size, the *observed* longest chance match
should track k*(M) = ln(Q*M)/ln(1/p): a straight line vs ln(M) with slope
1/ln(1/p). That sweep is the figure for the paper.

WHAT THIS SCRIPT DOES
---------------------
Generates random "queries" from the Robinson & Robinson composition, computes each
query's longest exact substring match against a database, and compares the observed
distribution of L to predictions (1)-(2).

Two database models:
  * synthetic (default): a random string at the same composition -> a clean test
    of the *formula* (match probability is exactly Sum f_i^2 by construction).
  * real proteome (--proteome file.fasta): queries are drawn from the proteome's
    OWN composition and searched against it (so p = Sum g_i^2 of that proteome).
    Real sequence has repeats / low-complexity / homologous families that the
    i.i.d. model ignores; any excess of long matches over prediction (1) is a
    reportable finding -- it shows k* is a *conservative* chance baseline and by
    how much.

USAGE
-----
  # single synthetic point (fast demo):
  python validate_kstar.py --db-size 200000 --query-len 150 --n-queries 500

  # the paper figure: sweep DB size over orders of magnitude
  python validate_kstar.py --sweep 1e4,1e5,1e6,1e7 --n-queries 1000 --plot sweep.png

  # test on a real proteome (subsampled to --db-size residues)
  python validate_kstar.py --proteome human.fasta --db-size 2000000 --n-queries 500

  # emit the random queries as FASTA to run end-to-end through frankensearch
  python validate_kstar.py --db-size 200000 --emit-fasta random_queries.fasta

Every run also writes a TSV of the per-k observed-vs-predicted table (--out).
Only matplotlib is optional (for --plot); everything else is stdlib.
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
import sys

# Robinson & Robinson (1991) amino-acid background frequencies (fractions).
# Sum f_i^2 over these ~= 0.0587; the tool pins p = 0.0598. The tiny difference
# reflects the exact frequency table and does not affect the conclusion -- the
# script computes p from whatever composition it actually samples from, so the
# null is always internally consistent with its own k*.
RR_FREQ: dict[str, float] = {
    "A": 0.0780, "R": 0.0512, "N": 0.0448, "D": 0.0536, "C": 0.0193,
    "Q": 0.0426, "E": 0.0632, "G": 0.0737, "H": 0.0229, "I": 0.0514,
    "L": 0.0901, "K": 0.0574, "M": 0.0224, "F": 0.0385, "P": 0.0520,
    "S": 0.0712, "T": 0.0584, "W": 0.0133, "Y": 0.0322, "V": 0.0673,
}

TOOL_P = 0.0598  # frankensearch's pinned background match probability
STANDARD_AA = set(RR_FREQ)
SENTINEL = "\x00"  # joins proteins so no exact match spans a protein boundary


# --------------------------------------------------------------------------- #
# Composition and theory
# --------------------------------------------------------------------------- #
def normalized(freq: dict[str, float]) -> dict[str, float]:
    total = sum(freq.values())
    return {k: v / total for k, v in freq.items()}


def match_probability(freq: dict[str, float]) -> float:
    """p = Sum f_i^2: probability two independently drawn residues are identical."""
    return sum(f * f for f in normalized(freq).values())


def k_star(m_residues: float, query_len: int, p: float, *, rounded: bool = True) -> float:
    """k* = ln(M*Q) / ln(1/p); rounded matches the tool's reported integer."""
    k = math.log(m_residues * query_len) / math.log(1.0 / p)
    return round(k) if rounded else k


def predicted_ge(k: float, m_residues: float, query_len: int, p: float) -> float:
    """P(L >= k) ~= 1 - exp(-Q*M*p^k) under the i.i.d. null."""
    return 1.0 - math.exp(-query_len * m_residues * p**k)


def predicted_cdf(k: int, m_residues: float, query_len: int, p: float) -> float:
    """P(L <= k) ~= exp(-Q*M*p^(k+1))."""
    return math.exp(-query_len * m_residues * p ** (k + 1))


def predicted_mean(m_residues: float, query_len: int, p: float, kmax: int = 60) -> float:
    """E[L] = Sum_k k * (F(k) - F(k-1)) from the predicted CDF."""
    prev = 0.0
    mean = 0.0
    for k in range(0, kmax + 1):
        cdf = predicted_cdf(k, m_residues, query_len, p)
        mean += k * (cdf - prev)
        prev = cdf
    return mean


# --------------------------------------------------------------------------- #
# Sequence generation and the longest-exact-match core
# --------------------------------------------------------------------------- #
def random_seq(n: int, letters: list[str], weights: list[float], rng: random.Random) -> str:
    return "".join(rng.choices(letters, weights=weights, k=n))


def longest_exact_match(query: str, db_joined: str) -> int:
    """Length of the longest substring of `query` that occurs in `db_joined`.

    Binary search on length L: "some length-L substring of the query occurs in the
    DB" is monotone in L (a substring of a match is itself a match), so we can
    bisect. `sub in db_joined` uses CPython's fast C substring search; the sentinel
    join guarantees a hit never spans a protein boundary (queries contain no
    sentinel). Memory-lean -- holds only the DB string, so it scales to real
    proteomes.
    """
    q = len(query)

    def exists(length: int) -> bool:
        if length == 0:
            return True
        return any(query[i:i + length] in db_joined for i in range(q - length + 1))

    lo, hi = 0, q
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if exists(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


# --------------------------------------------------------------------------- #
# Proteome loading
# --------------------------------------------------------------------------- #
def load_proteome(path: str, db_size: int, rng: random.Random) -> list[str]:
    """Read FASTA, keep standard-AA runs, subsample to ~db_size residues.

    k* depends on M only through ln(M), so subsampling a huge proteome to a couple
    of million residues shifts k* by well under one residue while keeping the run
    fast. Non-standard symbols (X/U/B/Z/*) split a sequence so we never build exact
    matches through an ambiguity code.
    """
    proteins: list[str] = []
    current: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    proteins.append("".join(current))
                    current = []
                continue
            current.append(line.strip().upper())
    if current:
        proteins.append("".join(current))

    # Split on any non-standard residue so matches can't run through it.
    clean: list[str] = []
    for prot in proteins:
        run: list[str] = []
        for ch in prot:
            if ch in STANDARD_AA:
                run.append(ch)
            elif run:
                clean.append("".join(run))
                run = []
        if run:
            clean.append("".join(run))

    rng.shuffle(clean)
    picked: list[str] = []
    total = 0
    for prot in clean:
        picked.append(prot)
        total += len(prot)
        if total >= db_size:
            break
    if not picked:
        sys.exit(f"error: no usable standard-AA sequence found in {path}")
    return picked


def composition(proteins: list[str]) -> dict[str, float]:
    counts = {aa: 0 for aa in RR_FREQ}
    for prot in proteins:
        for ch in prot:
            if ch in counts:
                counts[ch] += 1
    total = sum(counts.values()) or 1
    return {aa: c / total for aa, c in counts.items()}


# --------------------------------------------------------------------------- #
# One experimental point
# --------------------------------------------------------------------------- #
class Point:
    def __init__(self, m_residues: int, query_len: int, p: float, lengths: list[int]):
        self.m = m_residues
        self.q = query_len
        self.p = p
        self.lengths = lengths
        self.k_star = k_star(m_residues, query_len, p, rounded=True)
        self.k_star_real = k_star(m_residues, query_len, p, rounded=False)
        self.mean = statistics.mean(lengths)
        self.median = statistics.median(lengths)
        self.frac_ge_kstar = sum(1 for x in lengths if x >= self.k_star) / len(lengths)


def run_point(
    m_residues: int,
    query_len: int,
    n_queries: int,
    freq: dict[str, float],
    rng: random.Random,
    db_string: str | None = None,
) -> tuple[Point, list[str]]:
    """Build (or reuse) a DB, sample queries, measure longest exact match each."""
    freq_n = normalized(freq)
    letters = list(freq_n)
    weights = [freq_n[a] for a in letters]
    p = match_probability(freq)

    if db_string is None:
        db_string = random_seq(m_residues, letters, weights, rng)

    queries = [random_seq(query_len, letters, weights, rng) for _ in range(n_queries)]
    lengths = [longest_exact_match(q, db_string) for q in queries]
    return Point(m_residues, query_len, p, lengths), queries


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def report_point(pt: Point) -> list[str]:
    """Human-readable observed-vs-predicted table for one (M, Q) point."""
    rows: list[str] = []
    rows.append(f"  M (DB residues) : {pt.m:,}")
    rows.append(f"  Q (query length): {pt.q}")
    rows.append(f"  p (Sum f_i^2)   : {pt.p:.4f}   (tool pins {TOOL_P})")
    rows.append(f"  k* (tool round) : {pt.k_star}   (real-valued {pt.k_star_real:.3f})")
    rows.append("")
    rows.append(f"  observed mean L : {pt.mean:.3f}   predicted {predicted_mean(pt.m, pt.q, pt.p):.3f}")
    rows.append(f"  observed median : {pt.median}")
    rows.append(
        f"  P(L >= k*)      : observed {pt.frac_ge_kstar:.3f}   "
        f"predicted {predicted_ge(pt.k_star, pt.m, pt.q, pt.p):.3f}   "
        f"(anchor at real k*: {1 - 1 / math.e:.3f})"
    )
    rows.append("")
    rows.append("  k     obs P(L>=k)   pred P(L>=k)   obs count   <- k* row marked *")
    n = len(pt.lengths)
    lo = max(1, pt.k_star - 4)
    hi = pt.k_star + 6
    for k in range(lo, hi + 1):
        obs = sum(1 for x in pt.lengths if x >= k) / n
        pred = predicted_ge(k, pt.m, pt.q, pt.p)
        mark = " *" if k == pt.k_star else ""
        cnt = sum(1 for x in pt.lengths if x == k)
        rows.append(f"  {k:<4}  {obs:>10.3f}   {pred:>12.3f}   {cnt:>9}{mark}")
    return rows


def write_tsv(path: str, points: list[Point]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "M\tQ\tp\tk_star\tk_star_real\tobs_mean_L\tobs_median_L\t"
            "obs_frac_ge_kstar\tpred_mean_L\tpred_frac_ge_kstar\n"
        )
        for pt in points:
            fh.write(
                f"{pt.m}\t{pt.q}\t{pt.p:.5f}\t{pt.k_star}\t{pt.k_star_real:.4f}\t"
                f"{pt.mean:.4f}\t{pt.median}\t{pt.frac_ge_kstar:.4f}\t"
                f"{predicted_mean(pt.m, pt.q, pt.p):.4f}\t"
                f"{predicted_ge(pt.k_star, pt.m, pt.q, pt.p):.4f}\n"
            )


def plot_sweep(points: list[Point], path: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed; skipping --plot)", file=sys.stderr)
        return
    ks = [pt.k_star_real for pt in points]
    obs = [pt.mean for pt in points]
    lo, hi = min(ks) - 1, max(ks) + 1
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot([lo, hi], [lo, hi], "--", color="grey", label="observed = k*")
    ax.scatter(ks, obs, color="C0", zorder=3, label="observed mean longest match")
    for pt in points:
        ax.annotate(f"M={pt.m:,}", (pt.k_star_real, pt.mean),
                    textcoords="offset points", xytext=(6, -4), fontsize=8)
    ax.set_xlabel("predicted k* = ln(M·Q)/ln(1/p)")
    ax.set_ylabel("observed longest exact match (mean over queries)")
    ax.set_title("Random-query null vs k*")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"  wrote plot: {path}")


def write_fasta(path: str, queries: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for i, seq in enumerate(queries, 1):
            fh.write(f">random_query_{i}\n{seq}\n")
    print(f"  wrote {len(queries)} random queries: {path}")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_sizes(spec: str) -> list[int]:
    return [int(float(tok)) for tok in spec.split(",") if tok.strip()]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Validate FRANKENSEARCH k* against a random-query null.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--db-size", type=lambda s: int(float(s)), default=200_000,
                    help="synthetic DB residues, or subsample size for --proteome (default 200000)")
    ap.add_argument("--query-len", type=int, default=150,
                    help="length of each random query (default 150)")
    ap.add_argument("--n-queries", type=int, default=500,
                    help="number of random queries per point (default 500)")
    ap.add_argument("--sweep", type=str, default=None,
                    help="comma-separated DB sizes for the k*(M) figure, e.g. 1e4,1e5,1e6,1e7")
    ap.add_argument("--proteome", type=str, default=None,
                    help="FASTA proteome to search (real-composition test) instead of synthetic")
    ap.add_argument("--seed", type=int, default=1,
                    help="RNG seed for reproducibility (default 1)")
    ap.add_argument("--out", type=str, default="kstar_validation.tsv",
                    help="TSV summary output (default kstar_validation.tsv)")
    ap.add_argument("--plot", type=str, default=None,
                    help="write a PNG figure (needs matplotlib); sweep -> mean L vs k*")
    ap.add_argument("--emit-fasta", type=str, default=None,
                    help="write the random queries to FASTA for an end-to-end frankensearch check")
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)

    # Choose the composition / database model.
    if args.proteome:
        proteins = load_proteome(args.proteome, args.db_size, rng)
        freq = composition(proteins)
        db_string = SENTINEL.join(proteins)
        m_residues = sum(len(p) for p in proteins)
        print(f"proteome: {args.proteome}  ->  {len(proteins):,} sequences, "
              f"{m_residues:,} residues (subsampled to ~{args.db_size:,})")
        if args.sweep:
            print("note: --sweep is ignored with --proteome (single real-DB point).",
                  file=sys.stderr)
    else:
        freq = RR_FREQ
        db_string = None
        m_residues = args.db_size

    points: list[Point] = []
    first_queries: list[str] = []

    if args.sweep and not args.proteome:
        for m in parse_sizes(args.sweep):
            pt, queries = run_point(m, args.query_len, args.n_queries, freq, rng)
            points.append(pt)
            if not first_queries:
                first_queries = queries
            print(f"\n=== synthetic DB, M = {m:,} ===")
            print("\n".join(report_point(pt)))
    else:
        pt, queries = run_point(m_residues, args.query_len, args.n_queries, freq, rng,
                                db_string=db_string)
        points.append(pt)
        first_queries = queries
        header = "real proteome" if args.proteome else "synthetic DB"
        print(f"\n=== {header}, M = {pt.m:,} ===")
        print("\n".join(report_point(pt)))

    write_tsv(args.out, points)
    print(f"\nwrote TSV: {args.out}")
    if args.plot:
        plot_sweep(points, args.plot)
    if args.emit_fasta:
        write_fasta(args.emit_fasta, first_queries)

    # Headline verdict for the default single-point run.
    if len(points) == 1:
        pt = points[0]
        gap = pt.mean - pt.k_star_real
        pred_ge = predicted_ge(pt.k_star, pt.m, pt.q, pt.p)
        dev = pt.frac_ge_kstar - pred_ge
        print(f"\nVERDICT: observed mean longest match {pt.mean:.2f} vs k* {pt.k_star_real:.2f} "
              f"(diff {gap:+.2f} residue); P(L>=k*) observed {pt.frac_ge_kstar:.3f} vs "
              f"i.i.d. {pred_ge:.3f} (dev {dev:+.3f}).")
        if args.proteome:
            if dev > 0.02:
                print("  Excess over the i.i.d. baseline: repeats / low-complexity inflate the"
                      " longest chance match; k* is then a conservative floor.")
            elif dev < -0.02:
                print("  Slight deficit vs i.i.d.: k-mer redundancy in real sequence lowers the"
                      " chance a RANDOM query matches, so k* is faithful (mildly conservative).")
            else:
                print("  Composition-matched chance background matches the i.i.d. k* prediction.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
