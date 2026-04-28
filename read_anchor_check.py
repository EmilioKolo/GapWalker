#!/usr/bin/env python3

"""
read_anchor_check.py

Functions for aligning two short sequences (FASTA) against a single
read extracted from a FASTQ file. Produces:
  - <output_prefix>_alignment.png   : linear map visualization
  - <output_prefix>_alignment.txt   : classic text alignment view
  - <output_prefix>_summary.csv     : per-alignment stats table

Dependencies:
  pip install biopython matplotlib
  samtools (must be on PATH)
"""

import csv
import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Optional

from Bio import Align, SeqIO
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


# Internal data structures

@dataclass
class AlignmentResult:
    label: str
    score: float
    ref_start: int       # 0-based, on the reference
    ref_end: int         # exclusive
    query_start: int     # 0-based, on the query
    query_end: int       # exclusive
    query_len: int
    aligned_ref: str     # gapped reference string
    aligned_query: str   # gapped query string
    n_matches: int
    n_mismatches: int
    n_gaps: int


# FASTQ indexing + extraction (samtools)

def _ensure_index(fastq_path: str) -> None:
    """Create <fastq_path>.fai if it does not already exist."""
    index_path = fastq_path + ".fai"
    if not os.path.exists(index_path):
        result = subprocess.run(
            ["samtools", "fqidx", fastq_path],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"samtools fqidx failed for '{fastq_path}':\n{result.stderr}"
            )


def _extract_read(fastq_path: str, read_id: str) -> str:
    """
    Return the nucleotide sequence for *read_id* from *fastq_path*.
    Assumes the index exists (call _ensure_index first).
    Uses a temp FASTA file so BioPython can parse the output cleanly.
    The temp file is removed before returning.
    """
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".fa")
    os.close(tmp_fd)
    try:
        result = subprocess.run(
            ["samtools", "faidx", fastq_path, read_id],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not result.stdout.strip():
            raise RuntimeError(
                f"samtools faidx could not extract read '{read_id}' "
                f"from '{fastq_path}':\n{result.stderr}"
            )
        with open(tmp_path, "w") as fh:
            fh.write(result.stdout)
        record = SeqIO.read(tmp_path, "fasta")
        return str(record.seq).upper()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


# Alignment

def _build_aligner(
    match_score: float,
    mismatch_score: float,
    gap_open: float,
    gap_extend: float,
) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    return aligner


def _parse_aligned_strings(alignment) -> tuple[str, str]:
    """
    Extract gapped aligned strings (query, target) from a BioPython
    PairwiseAlignment object in a version-safe way.
    """
    fmt = format(alignment, "fasta")
    lines = [l for l in fmt.splitlines() if not l.startswith(">")]
    # fasta format interleaves query/target as two sequences
    if len(lines) >= 2:
        return lines[0].upper(), lines[1].upper()
    # fallback: use the alignment string representation
    parts = str(alignment).split("\n")
    seqs = [p for p in parts if p and not p.startswith(" ")]
    return (seqs[0].upper() if len(seqs) > 0 else "",
            seqs[-1].upper() if len(seqs) > 1 else "")


def _count_alignment_stats(aligned_query: str, aligned_ref: str) -> tuple[int, int, int]:
    """Return (n_matches, n_mismatches, n_gaps) from two gapped strings."""
    matches = mismatches = gaps = 0
    for q, r in zip(aligned_query, aligned_ref):
        if q == "-" or r == "-":
            gaps += 1
        elif q == r:
            matches += 1
        else:
            mismatches += 1
    return matches, mismatches, gaps


def _align_one(
    label: str,
    query_seq: str,
    ref_seq: str,
    aligner: Align.PairwiseAligner,
    print_debug: bool = False,
) -> AlignmentResult:
    alignments:Align.PairwiseAlignments = aligner.align(ref_seq, query_seq)
    best = next(iter(alignments))

    # Print type of object for "best"
    if print_debug:
        print(f"Best alignment type: {type(best)}")

    # Coordinates on reference (target) and query
    ref_coords = best.aligned[0]    # array of (start, end) blocks on ref
    qry_coords = best.aligned[1]    # array of (start, end) blocks on query

    ref_start = int(ref_coords[0][0])
    ref_end   = int(ref_coords[-1][1])
    qry_start = int(qry_coords[0][0])
    qry_end   = int(qry_coords[-1][1])

    aligned_query, aligned_ref = _parse_aligned_strings(best)
    matches, mismatches, gaps = _count_alignment_stats(aligned_query, aligned_ref)

    return AlignmentResult(
        label=label,
        score=best.score,
        ref_start=ref_start,
        ref_end=ref_end,
        query_start=qry_start,
        query_end=qry_end,
        query_len=len(query_seq),
        aligned_ref=aligned_ref,
        aligned_query=aligned_query,
        n_matches=matches,
        n_mismatches=mismatches,
        n_gaps=gaps,
    )


# Outputs

# Visual identity for each sequence — distinct colors, hatching, and labels
_SEQ_STYLES = [
    {
        "color":     "#2171b5",   # blue
        "edge":      "#08519c",
        "hatch":     None,
        "mismatch":  "#fdae6b",   # orange highlight for mismatches
    },
    {
        "color":     "#238b45",   # green
        "edge":      "#005a32",
        "hatch":     "///",
        "mismatch":  "#d94801",   # red-orange highlight for mismatches
    },
]


def _write_text_alignment(results: list[AlignmentResult], output_prefix: str) -> str:
    """Write a classic gapped-string alignment view.  Returns the file path."""
    out_path = output_prefix + "_alignment.txt"
    line_width = 80

    with open(out_path, "w") as fh:
        for res in results:
            fh.write("=" * 80 + "\n")
            fh.write(f"Sequence : {res.label}\n")
            fh.write(f"Score    : {res.score:.1f}\n")
            fh.write(f"Ref span : {res.ref_start + 1}–{res.ref_end}  "
                     f"(length {res.ref_end - res.ref_start} bp)\n")
            fh.write(f"Query span: {res.query_start + 1}–{res.query_end}  "
                     f"(of {res.query_len} bp)\n")
            fh.write(f"Matches  : {res.n_matches}  |  "
                     f"Mismatches: {res.n_mismatches}  |  "
                     f"Gaps: {res.n_gaps}\n")
            fh.write("-" * 80 + "\n\n")

            aq, ar = res.aligned_query, res.aligned_ref

            # Build match string ( | match, X mismatch, space gap)
            match_str = ""
            for q, r in zip(aq, ar):
                if q == "-" or r == "-":
                    match_str += " "
                elif q == r:
                    match_str += "|"
                else:
                    match_str += "X"

            ref_pos   = res.ref_start
            query_pos = res.query_start

            for i in range(0, len(aq), line_width):
                chunk_q = aq[i: i + line_width]
                chunk_m = match_str[i: i + line_width]
                chunk_r = ar[i: i + line_width]

                q_end = query_pos + len(chunk_q.replace("-", ""))
                r_end = ref_pos   + len(chunk_r.replace("-", ""))

                fh.write(f"Query {query_pos + 1:>8}  {chunk_q}  {q_end}\n")
                fh.write(f"           {'':>8}  {chunk_m}\n")
                fh.write(f"Ref   {ref_pos + 1:>8}  {chunk_r}  {r_end}\n\n")

                query_pos = q_end
                ref_pos   = r_end

            fh.write("\n")

    return out_path


def _write_summary_csv(results: list[AlignmentResult], ref_len: int, output_prefix: str) -> str:
    """Write a per-alignment stats CSV.  Returns the file path."""
    out_path = output_prefix + "_summary.csv"
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "sequence", "score",
            "ref_start_1based", "ref_end", "ref_span_bp",
            "query_start_1based", "query_end", "query_len",
            "matches", "mismatches", "gaps",
            "identity_pct", "query_coverage_pct",
        ])
        for res in results:
            aligned_len = res.n_matches + res.n_mismatches + res.n_gaps
            identity    = 100.0 * res.n_matches / aligned_len if aligned_len else 0.0
            coverage    = 100.0 * (res.query_end - res.query_start) / res.query_len \
                          if res.query_len else 0.0
            writer.writerow([
                res.label,
                f"{res.score:.1f}",
                res.ref_start + 1,
                res.ref_end,
                res.ref_end - res.ref_start,
                res.query_start + 1,
                res.query_end,
                res.query_len,
                res.n_matches,
                res.n_mismatches,
                res.n_gaps,
                f"{identity:.2f}",
                f"{coverage:.2f}",
            ])
    return out_path


def _write_plot(
    results: list[AlignmentResult],
    ref_len: int,
    read_id: str,
    output_prefix: str,
) -> str:
    """Draw a linear alignment map and save as PNG. Returns the file path."""
    out_path = output_prefix + "_alignment.png"

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(14, 5),
        gridspec_kw={"height_ratios": [1, 3]},
    )
    fig.patch.set_facecolor("#f8f9fa")

    # ---- Top panel: overview (full reference, both alignments as bars) ----
    ax_top = axes[0]
    ax_top.set_facecolor("#f0f2f5")
    ax_top.set_xlim(0, ref_len)
    ax_top.set_ylim(0, 1)
    ax_top.set_yticks([])
    ax_top.set_xlabel("Reference position (bp)", fontsize=9, color="#444")
    ax_top.set_title(f"Alignment overview — reference: {read_id}", fontsize=10,
                     color="#222", pad=6)

    # Reference backbone
    ax_top.barh(0.5, ref_len, left=0, height=0.18, color="#cccccc",
                linewidth=0, zorder=1)

    for i, res in enumerate(results):
        style = _SEQ_STYLES[i % len(_SEQ_STYLES)]
        ax_top.barh(
            0.5,
            res.ref_end - res.ref_start,
            left=res.ref_start,
            height=0.38,
            color=style["color"],
            edgecolor=style["edge"],
            linewidth=0.6,
            hatch=style["hatch"],
            alpha=0.85,
            label=res.label,
            zorder=2,
        )

    ax_top.legend(
        loc="upper right", fontsize=8, framealpha=0.7,
        handlelength=1.4, handleheight=0.9,
    )
    ax_top.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{int(x):,}"
    ))

    # ---- Bottom panel: per-base detail (zoomed to the aligned region) ----
    ax_bot = axes[1]
    ax_bot.set_facecolor("#f0f2f5")
    ax_bot.set_yticks([])

    # Determine zoom window (union of both alignment spans + 2 % margin)
    span_start = min(r.ref_start for r in results)
    span_end   = max(r.ref_end   for r in results)
    margin     = max(20, int(0.02 * (span_end - span_start)))
    zoom_start = max(0, span_start - margin)
    zoom_end   = min(ref_len, span_end + margin)

    ax_bot.set_xlim(zoom_start, zoom_end)
    ax_bot.set_ylim(0, 1)
    ax_bot.set_xlabel("Reference position (bp)", fontsize=9, color="#444")
    ax_bot.set_title("Zoomed view with mismatch / gap annotation",
                     fontsize=10, color="#222", pad=6)

    # Reference backbone in zoomed window
    ax_bot.barh(0.5, zoom_end - zoom_start, left=zoom_start, height=0.10,
                color="#cccccc", linewidth=0, zorder=1)

    n = len(results)
    bar_height = 0.22
    # Spread the two sequences vertically so they don't overlap
    y_centers = [0.68, 0.32] if n > 1 else [0.5]

    for i, res in enumerate(results):
        style = _SEQ_STYLES[i % len(_SEQ_STYLES)]
        yc    = y_centers[i]

        # Main alignment bar
        ax_bot.barh(
            yc,
            res.ref_end - res.ref_start,
            left=res.ref_start,
            height=bar_height,
            color=style["color"],
            edgecolor=style["edge"],
            linewidth=0.6,
            hatch=style["hatch"],
            alpha=0.80,
            zorder=2,
            label=res.label,
        )

        # Mismatch highlights — iterate over the aligned columns
        ref_pos = res.ref_start
        for q_base, r_base in zip(res.aligned_query, res.aligned_ref):
            if r_base != "-":          # advance ref cursor only on non-ref-gap
                if q_base != "-" and q_base != r_base:
                    ax_bot.barh(
                        yc,
                        1,
                        left=ref_pos,
                        height=bar_height,
                        color=style["mismatch"],
                        linewidth=0,
                        zorder=3,
                        alpha=0.9,
                    )
                ref_pos += 1

        # Gap markers — thin vertical lines at gap positions
        ref_pos = res.ref_start
        for q_base, r_base in zip(res.aligned_query, res.aligned_ref):
            if q_base == "-":
                ax_bot.vlines(
                    ref_pos,
                    yc - bar_height / 2,
                    yc + bar_height / 2,
                    color="white",
                    linewidth=0.8,
                    zorder=4,
                )
            if r_base != "-":
                ref_pos += 1

        # Sequence label on the left of each bar
        ax_bot.text(
            zoom_start + (zoom_end - zoom_start) * 0.005,
            yc + bar_height / 2 + 0.03,
            res.label,
            fontsize=8,
            color=style["edge"],
            va="bottom",
            zorder=5,
        )

        # Stats annotation on the right
        aligned_len = res.n_matches + res.n_mismatches + res.n_gaps
        identity    = 100.0 * res.n_matches / aligned_len if aligned_len else 0.0
        ax_bot.text(
            zoom_end - (zoom_end - zoom_start) * 0.005,
            yc,
            f"{identity:.1f}% id  |  {res.n_mismatches} mm  |  {res.n_gaps} gaps",
            fontsize=7.5,
            color="#444",
            va="center",
            ha="right",
            zorder=5,
        )

    ax_bot.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{int(x):,}"
    ))

    # Shared legend for mismatch / gap symbols
    mismatch_patches = [
        mpatches.Patch(color=_SEQ_STYLES[i]["mismatch"],
                       label=f"{results[i].label} mismatch")
        for i in range(len(results))
    ]
    ax_bot.legend(
        handles=mismatch_patches,
        loc="lower right",
        fontsize=7.5,
        framealpha=0.7,
        title="Highlighted positions",
        title_fontsize=7.5,
    )

    plt.tight_layout(h_pad=1.5)
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return out_path


# Entry points

def check_sequence_presence(
    fasta_query: str,
    fasta_target: str,
    label_query: Optional[str] = None,
    label_target: Optional[str] = None,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> dict:
    """
    Evaluate whether a short query sequence is present within a longer target
    sequence, both supplied as single-record FASTA files.

    Uses local pairwise alignment (Smith-Waterman via BioPython). The query
    is considered "present" when both identity and query coverage exceed the
    supplied thresholds.

    Parameters
    ----------
    fasta_query : str
        Path to the short query sequence (FASTA, single record).
    fasta_target : str
        Path to the longer target sequence (FASTA, single record).
    label_query : str, optional
        Display name for the query. Defaults to the FASTA record ID.
    label_target : str, optional
        Display name for the target. Defaults to the FASTA record ID.
    match_score : float
        Score for a matching base (default 2.0).
    mismatch_score : float
        Penalty for a mismatched base (default -1.0).
    gap_open : float
        Gap-open penalty (default -10.0).
    gap_extend : float
        Gap-extension penalty per base (default -0.5).
    min_identity : float
        Minimum percent identity (0-100) to consider the query present
        (default 90.0).
    min_coverage : float
        Minimum percent of the query that must be aligned (0-100) to
        consider the query present (default 80.0).

    Returns
    -------
    dict with keys:
        "present"          : bool   - True when both thresholds are met
        "score"            : float  - Raw alignment score
        "identity_pct"     : float  - % identical bases in the aligned region
        "coverage_pct"     : float  - % of query bases covered by the alignment
        "target_start"     : int    - 1-based start position on the target
        "target_end"       : int    - 1-based end position on the target (inclusive)
        "query_start"      : int    - 1-based start position on the query
        "query_end"        : int    - 1-based end position on the query (inclusive)
        "n_matches"        : int
        "n_mismatches"     : int
        "n_gaps"           : int
        "query_len"        : int    - Full length of the query sequence
        "target_len"       : int    - Full length of the target sequence
        "label_query"      : str
        "label_target"     : str
        "aligned_query"    : str    - Gapped query string from best alignment
        "aligned_target"   : str    - Gapped target string from best alignment
    """
    rec_query  = SeqIO.read(fasta_query,  "fasta")
    rec_target = SeqIO.read(fasta_target, "fasta")

    seq_query  = str(rec_query.seq).upper()
    seq_target = str(rec_target.seq).upper()

    lbl_query  = label_query  or rec_query.id
    lbl_target = label_target or rec_target.id

    aligner = _build_aligner(match_score, mismatch_score, gap_open, gap_extend)

    # _align_one expects (label, query_seq, ref_seq, aligner) — same convention
    result = _align_one(lbl_query, seq_query, seq_target, aligner)

    aligned_len  = result.n_matches + result.n_mismatches + result.n_gaps
    identity_pct = 100.0 * result.n_matches / aligned_len if aligned_len else 0.0
    coverage_pct = 100.0 * (result.query_end - result.query_start) / result.query_len \
                   if result.query_len else 0.0

    present = (identity_pct >= min_identity) and (coverage_pct >= min_coverage)

    return {
        "present":       present,
        "score":         result.score,
        "identity_pct":  round(identity_pct, 4),
        "coverage_pct":  round(coverage_pct, 4),
        "target_start":  result.ref_start + 1,   # convert to 1-based
        "target_end":    result.ref_end,         # ref_end is already exclusive → last included base
        "query_start":   result.query_start + 1,
        "query_end":     result.query_end,
        "n_matches":     result.n_matches,
        "n_mismatches":  result.n_mismatches,
        "n_gaps":        result.n_gaps,
        "query_len":     result.query_len,
        "target_len":    len(seq_target),
        "label_query":   lbl_query,
        "label_target":  lbl_target,
        "aligned_query":  result.aligned_query,
        "aligned_target": result.aligned_ref,
    }

def run_dual_alignment(
    fasta_a: str,
    fasta_b: str,
    fastq_ref: str,
    read_id: str,
    output_prefix: str,
    label_a: Optional[str] = None,
    label_b: Optional[str] = None,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> dict:
    """
    Align two short sequences (FASTA) against one read from a FASTQ file.

    Parameters
    ----------
    fasta_a : str
        Path to the first short sequence (FASTA, single record).
    fasta_b : str
        Path to the second short sequence (FASTA, single record).
    fastq_ref : str
        Path to the FASTQ file containing the reference read.
        A samtools .fai index is created next to it if absent.
    read_id : str
        Read ID to extract from *fastq_ref* as the alignment target.
    output_prefix : str
        All output files will start with this string.
        E.g. "results/run1" → "results/run1_alignment.png", etc.
    label_a : str, optional
        Display name for sequence A. Defaults to the FASTA record ID.
    label_b : str, optional
        Display name for sequence B. Defaults to the FASTA record ID.
    match_score : float
        Score for a matching base (default 2.0).
    mismatch_score : float
        Penalty for a mismatched base (default -1.0).
    gap_open : float
        Gap-open penalty (default -10.0; high value → fewer gaps).
    gap_extend : float
        Gap-extension penalty per base (default -0.5).

    Returns
    -------
    dict with keys:
        "plot"    : path to the PNG visualization
        "text"    : path to the text alignment file
        "csv"     : path to the summary CSV
        "results" : list of AlignmentResult dataclass instances
    """
    # -- 1. Index + extract reference read ----------------------------------
    _ensure_index(fastq_ref)
    ref_seq = _extract_read(fastq_ref, read_id)
    ref_len = len(ref_seq)

    # -- 2. Load query sequences --------------------------------------------
    rec_a = SeqIO.read(fasta_a, "fasta")
    rec_b = SeqIO.read(fasta_b, "fasta")

    seq_a = str(rec_a.seq).upper()
    seq_b = str(rec_b.seq).upper()

    lbl_a = label_a or rec_a.id
    lbl_b = label_b or rec_b.id

    # -- 3. Build aligner and run -------------------------------------------
    aligner = _build_aligner(match_score, mismatch_score, gap_open, gap_extend)

    result_a = _align_one(lbl_a, seq_a, ref_seq, aligner)
    result_b = _align_one(lbl_b, seq_b, ref_seq, aligner)
    results  = [result_a, result_b]

    # -- 4. Write outputs ---------------------------------------------------
    # Ensure output directory exists if prefix contains a path component
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plot_path = _write_plot(results, ref_len, read_id, output_prefix)
    text_path = _write_text_alignment(results, output_prefix)
    csv_path  = _write_summary_csv(results, ref_len, output_prefix)

    return {
        "plot":    plot_path,
        "text":    text_path,
        "csv":     csv_path,
        "results": results,
    }