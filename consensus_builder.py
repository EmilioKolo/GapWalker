#!/usr/bin/env python3

"""
ONT Long-Read Consensus Builder
================================
Takes 10+ sequences (already trimmed to start AFTER a shared region)
and builds a consensus sequence by:
  1. Force-aligning all reads from position 0 (no soft-clipping at the start)
  2. Running a pairwise-guided progressive MSA (via MAFFT or parasail fallback)
  3. Calling consensus at every column with configurable thresholds
  4. Stopping extension when coverage drops below a minimum

Usage:
    python consensus_builder.py -i reads.fasta -o consensus.fasta [options]
    python consensus_builder.py -i reads.fasta -o consensus.fasta --min-coverage 3 --min-identity 0.4

Dependencies:
    pip install biopython parasail numpy
    Optional (better MSA): mafft (system install: apt/brew/conda)
"""

import argparse
import subprocess
import sys
import os
import tempfile
import csv
import numpy as np
from collections import Counter

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit("BioPython not found. Run: pip install biopython")

try:
    import parasail
    HAS_PARASAIL = True
except ImportError:
    HAS_PARASAIL = False
    print("[WARN] parasail not found (pip install parasail). "
          "Pairwise gap-open/gap-extend scoring will use basic NW.")


# ─── Pairwise alignment ──────────────────────────────────────────────────────

def align_pair_parasail(ref: str, query: str,
                        match=2, mismatch=-4,
                        gap_open=8, gap_extend=4) -> tuple[str, str]:
    """Global NW alignment anchored at position 0 using parasail."""
    matrix = parasail.matrix_create("ACGTN", match, mismatch)
    result = parasail.nw_trace_striped_32(query, ref, gap_open, gap_extend, matrix)
    cigar = result.cigar
    ref_al, qry_al = [], []
    ri, qi = 0, 0
    for item in cigar.seq:
        op = item & 0xF
        length = item >> 4
        op_char = cigar.decode_op(op)  # returns '=', 'X', 'I', 'D', etc.
        if op_char == "=" or op_char == "X":
            ref_al.append(ref[ri:ri+length])
            qry_al.append(query[qi:qi+length])
            ri += length; qi += length
        elif op_char == "D":
            ref_al.append(ref[ri:ri+length])
            qry_al.append("-" * length)
            ri += length
        elif op_char == "I":
            ref_al.append("-" * length)
            qry_al.append(query[qi:qi+length])
            qi += length
    return "".join(ref_al), "".join(qry_al)


def align_pair_biopython(ref: str, query: str) -> tuple[str, str]:
    """Fallback: BioPython PairwiseAligner (global, left-anchored)."""
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -2
    # Force left anchor: no end gap penalty on the right (reads differ in length)
    aligner.target_right_gap_score = 0
    aligner.query_right_gap_score = 0
    alignments = aligner.align(ref, query)
    best = next(iter(alignments))
    return str(best[0]), str(best[1])


def align_pair(ref: str, query: str) -> tuple[str, str]:
    if HAS_PARASAIL:
        return align_pair_parasail(ref, query)
    return align_pair_biopython(ref, query)


# ─── MSA ─────────────────────────────────────────────────────────────────────

def mafft_available() -> bool:
    try:
        subprocess.run(["mafft", "--version"],
                       capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def run_mafft(sequences: list[SeqRecord]) -> list[str]:
    """Run MAFFT with --adjustdirectionaccurately and return aligned strings."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        SeqIO.write(sequences, fh, "fasta")
        tmp_in = fh.name
    tmp_out = tmp_in + ".aln"
    cmd = [
        "mafft",
        "--auto",            # auto strategy (FFT-NS or L-INS-i based on size)
        "--adjustdirection", # handle rare reverse-complement edge cases
        "--thread", "-1",    # use all cores
        tmp_in
    ]
    with open(tmp_out, "w") as out_fh:
        subprocess.run(cmd, stdout=out_fh, stderr=subprocess.DEVNULL, check=True)
    aligned = [str(r.seq).upper() for r in SeqIO.parse(tmp_out, "fasta")]
    os.unlink(tmp_in); os.unlink(tmp_out)
    return aligned


def progressive_pairwise_msa(sequences: list[str]) -> list[str]:
    """
    Progressive MSA without MAFFT:
    Align each sequence against a growing consensus profile, re-inserting
    gap columns into all previous sequences when insertions are found.
    """
    print("  [MSA] Using progressive pairwise alignment (no MAFFT)...")
    aligned = [sequences[0]]  # first sequence is anchor (no leading gaps)

    for i, seq in enumerate(sequences[1:], 1):
        print(f"  [MSA] Aligning sequence {i+1}/{len(sequences)} ...")
        # Build a gap-stripped consensus of current aligned sequences as reference
        ref = gapless_consensus(aligned)
        ref_al, qry_al = align_pair(ref, seq)

        # Map ref alignment columns back to position in the aligned matrix
        new_aligned = insert_gaps_from_alignment(aligned, ref, ref_al)
        # Pad query to same length
        qry_padded = pad_query_to_ref(ref_al, qry_al, len(new_aligned[0]))
        new_aligned.append(qry_padded)
        aligned = new_aligned

    # Final length normalisation
    max_len = max(len(s) for s in aligned)
    aligned = [s.ljust(max_len, "-") for s in aligned]
    return aligned


def gapless_consensus(aligned: list[str]) -> str:
    """Majority-vote consensus ignoring gap columns, no gaps in output."""
    if not aligned:
        return ""
    result = []
    for col in zip(*aligned):
        bases = [b for b in col if b != "-"]
        if bases:
            result.append(Counter(bases).most_common(1)[0][0])
    return "".join(result)


def insert_gaps_from_alignment(aligned: list[str],
                                ref_ungapped: str,
                                ref_aligned: str) -> list[str]:
    """
    When aligning a new sequence to the gapless consensus, insertions in the
    new sequence correspond to new gap columns that must be inserted into ALL
    existing aligned sequences at the correct position.
    """
    # Build a map: position in ref_ungapped → position(s) in ref_aligned
    # We need to know WHERE to insert extra gap columns
    new_cols_before: dict[int, int] = {}  # ref pos → how many new gaps before it
    ref_pos = 0
    for al_char in ref_aligned:
        if al_char == "-":  # insertion in query = gap in ref_aligned
            new_cols_before[ref_pos] = new_cols_before.get(ref_pos, 0) + 1
    # Rebuild each aligned sequence inserting gap columns at insertion sites
    result = []
    for seq in aligned:
        # seq is aligned against ref_ungapped (may already have gaps)
        new_seq = []
        ungapped_pos = 0
        for ch in seq:
            # Insert any new gap columns before this ref position
            extra = new_cols_before.get(ungapped_pos, 0)
            new_seq.append("-" * extra)
            new_seq.append(ch)
            if ch != "-":
                ungapped_pos += 1
        # Handle gaps at the very end
        extra = new_cols_before.get(ungapped_pos, 0)
        new_seq.append("-" * extra)
        result.append("".join(new_seq))
    return result


def pad_query_to_ref(ref_al: str, qry_al: str, target_len: int) -> str:
    """Return qry_al padded/trimmed to match the updated MSA column count."""
    # The query aligned string already has the correct gap structure;
    # just ensure it matches target_len by padding with trailing gaps.
    if len(qry_al) < target_len:
        return qry_al + "-" * (target_len - len(qry_al))
    return qry_al[:target_len]


# ─── Consensus calling ───────────────────────────────────────────────────────

def call_consensus(aligned_seqs: list[str],
                   min_coverage: int = 3,
                   min_base_fraction: float = 0.4,
                   gap_fraction_threshold: float = 0.5) -> tuple[str, list[dict]]:
    """
    Column-by-column consensus calling.

    Parameters
    ----------
    aligned_seqs        : list of gap-padded aligned strings (same length)
    min_coverage        : minimum number of non-gap bases required to call a position
    min_base_fraction   : a base must appear in ≥ this fraction of reads to be called
    gap_fraction_threshold : if gaps make up ≥ this fraction, skip the column

    Returns
    -------
    consensus   : consensus string (stops when coverage drops below threshold)
    stats       : per-position statistics list (useful for QC)
    """
    n_reads = len(aligned_seqs)
    ncols = len(aligned_seqs[0])
    consensus = []
    stats = []

    last_good_col = 0
    read_lengths = [len(s.rstrip("-")) for s in aligned_seqs]

    for col_i in range(ncols):
        col = [s[col_i] for s in aligned_seqs]
        bases = [b for b in col if b != "-"]
        active_reads = sum(1 for i in range(n_reads) if col_i < read_lengths[i])
        coverage = len(bases)
        gap_frac = (active_reads - coverage) / active_reads if active_reads > 0 else 1.0

        if active_reads < min_coverage:
            stats.append({"col": col_i, "coverage": coverage,
                          "base": None, "reason": "low_coverage"})
            break  # stop extension here

        if gap_frac >= gap_fraction_threshold:
            # This column is mostly gaps (likely a minority insertion) — skip it
            stats.append({"col": col_i, "coverage": coverage,
                          "base": "-", "reason": "majority_gap"})
            continue

        counter = Counter(bases)
        top_base, top_count = counter.most_common(1)[0]
        fraction = top_count / coverage

        if fraction < min_base_fraction:
            # No clear winner — use IUPAC ambiguity or N
            top2 = [b for b, _ in counter.most_common(2)]
            iupac = get_iupac(top2)
            stats.append({"col": col_i, "coverage": coverage,
                          "base": iupac, "reason": "ambiguous",
                          "top_base": top_base, "fraction": fraction})
            consensus.append(iupac)
        else:
            stats.append({"col": col_i, "coverage": coverage,
                          "base": top_base, "reason": "ok",
                          "fraction": fraction})
            consensus.append(top_base)

        last_good_col = col_i

    return "".join(consensus), stats


IUPAC_MAP = {
    frozenset("AC"): "M", frozenset("AG"): "R", frozenset("AT"): "W",
    frozenset("CG"): "S", frozenset("CT"): "Y", frozenset("GT"): "K",
    frozenset("ACG"): "V", frozenset("ACT"): "H", frozenset("AGT"): "D",
    frozenset("CGT"): "B", frozenset("ACGT"): "N",
}

def get_iupac(bases: list[str]) -> str:
    return IUPAC_MAP.get(frozenset(b.upper() for b in bases), "N")


def run_consensus_builder(
    input_file: str,
    output_file: str,
    min_coverage: int = 3,
    min_base_fraction: float = 0.4,
    gap_fraction_threshold: float = 0.5,
    no_mafft: bool = False,
    save_msa: str|None = None,
    stats: str|None = None
):
    """
    Programmatic entry point mirroring the CLI of consensus_builder.py.

    Parameters
    ----------
    input_file              : Path to input FASTA (reads trimmed past anchor)
    output_file             : Path to output FASTA for consensus
    min_coverage            : Min reads covering a position to include it
    min_base_fraction       : Min fraction for a base to be called (else IUPAC/N)
    gap_fraction_threshold  : Column gap fraction above which column is skipped
    no_mafft                : Force progressive pairwise MSA even if MAFFT available
    save_msa                : If set, save the MSA to this FASTA file
    stats                   : If set, save per-column stats to this TSV file
    """

    # ── Load sequences ───────────────────────────────────────────────────────
    print(f"\n[1/4] Loading sequences from {input_file} ...")
    records = list(SeqIO.parse(input_file, "fasta"))
    if len(records) < 2:
        raise ValueError("Need at least 2 sequences.")
    sequences = [str(r.seq).upper().replace("U", "T") for r in records]
    lengths = [len(s) for s in sequences]
    print(f"  Loaded {len(sequences)} sequences")
    print(f"  Lengths: min={min(lengths)}, max={max(lengths)}, "
          f"mean={int(np.mean(lengths))}")

    # ── MSA ──────────────────────────────────────────────────────────────────
    print("\n[2/4] Running multiple sequence alignment ...")
    if not no_mafft and mafft_available():
        print("  Using MAFFT (--auto) ...")
        aligned_seqs = run_mafft(records)
    else:
        if not no_mafft:
            print("  MAFFT not found — falling back to progressive pairwise MSA.")
        aligned_seqs = progressive_pairwise_msa(sequences)

    print(f"  MSA complete. Alignment length: {len(aligned_seqs[0])} columns")

    if save_msa:
        msa_records = [
            SeqRecord(Seq(s), id=records[i].id, description="aligned")
            for i, s in enumerate(aligned_seqs)
        ]
        SeqIO.write(msa_records, save_msa, "fasta")
        print(f"  MSA saved to {save_msa}")

    # ── Consensus calling ────────────────────────────────────────────────────
    print("\n[3/4] Calling consensus ...")
    consensus_seq, col_stats = call_consensus(
        aligned_seqs,
        min_coverage=min_coverage,
        min_base_fraction=min_base_fraction,
        gap_fraction_threshold=gap_fraction_threshold,
    )

    called = [s for s in col_stats if s.get("base") not in (None, "-")]
    ambiguous = sum(1 for s in called if s.get("reason") == "ambiguous")
    print(f"  Consensus length : {len(consensus_seq)} bp")
    print(f"  Ambiguous bases  : {ambiguous} ({100*ambiguous/max(len(consensus_seq),1):.1f}%)")

    low_cov = [s for s in col_stats if s["coverage"] < min_coverage]
    if low_cov:
        print(f"  Extension stopped at column {low_cov[0]['col']} "
              f"(coverage={low_cov[0]['coverage']} < {min_coverage})")
    else:
        print("  All columns had sufficient coverage throughout.")

    # Print and save a coverage plot
    cov_plot_file = stats.replace(".tsv", "_coverage_plot.txt") if stats else None
    print_coverage_plot(col_stats, output_file=cov_plot_file)

    # ── Save outputs ─────────────────────────────────────────────────────────
    print(f"\n[4/4] Saving outputs ...")
    out_record = SeqRecord(
        Seq(consensus_seq),
        id="consensus",
        description=(f"len={len(consensus_seq)} "
                     f"min_cov={min_coverage} "
                     f"n_reads={len(sequences)}")
    )
    SeqIO.write([out_record], output_file, "fasta")
    print(f"  Consensus saved to {output_file}")

    if stats:
        with open(stats, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["col", "coverage", "base", "reason", "fraction", "top_base"],
                extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(col_stats)
        print(f"  Per-column stats saved to {stats}")

    print("\nDone ✓\n")
    return consensus_seq, col_stats

# ─── Reporting ───────────────────────────────────────────────────────────────

def print_coverage_plot(stats: list[dict], width: int = 80, output_file=None):
    """
    ASCII coverage plot across the consensus.
    Also saves the figure to a text file if output_file is provided.
    """
    coverages = [s["coverage"] for s in stats if s.get("base") not in (None,)]
    if not coverages:
        return
    max_cov = max(coverages)
    step = max(1, len(coverages) // width)
    if output_file:
        with open(output_file, "w") as f:
            f.write("Coverage across consensus (ASCII plot):\n")
            f.write(f"Max coverage: {max_cov}  Length: {len(coverages)} bp\n\n")
            bar_height = 10
            for row in range(bar_height, 0, -1):
                threshold = (row / bar_height) * max_cov
                line = ""
                for i in range(0, len(coverages), step):
                    chunk = coverages[i:i+step]
                    avg = sum(chunk) / len(chunk)
                    line += "█" if avg >= threshold else " "
                label = f"{int(threshold):>4} |"
                f.write(f"{label}{line}\n")
            f.write(" " * 6 + "─" * (len(coverages) // step) + "\n")
            f.write(f"      0{' ' * (len(coverages) // step - 10)}{len(coverages)} bp\n")
    print("\n  Coverage across consensus (ASCII plot):")
    print(f"  Max coverage: {max_cov}  Length: {len(coverages)} bp\n")
    bar_height = 10
    for row in range(bar_height, 0, -1):
        threshold = (row / bar_height) * max_cov
        line = ""
        for i in range(0, len(coverages), step):
            chunk = coverages[i:i+step]
            avg = sum(chunk) / len(chunk)
            line += "█" if avg >= threshold else " "
        label = f"{int(threshold):>4} |"
        print(f"  {label}{line}")
    print("  " + " " * 6 + "─" * (len(coverages) // step))
    print(f"       0{' ' * (len(coverages) // step - 10)}{len(coverages)} bp")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build a consensus from ONT long reads anchored at position 0.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file with reads (already trimmed past anchor)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA file for consensus sequence")
    parser.add_argument("--min-coverage", type=int, default=3,
                        help="Min reads covering a position to include it")
    parser.add_argument("--min-base-fraction", type=float, default=0.4,
                        help="Min fraction for a base to be called (else IUPAC/N)")
    parser.add_argument("--gap-fraction-threshold", type=float, default=0.5,
                        help="Column gap fraction above which column is skipped (minority insertion)")
    parser.add_argument("--no-mafft", action="store_true",
                        help="Force progressive pairwise MSA even if MAFFT is available")
    parser.add_argument("--save-msa", metavar="FILE",
                        help="Save the MSA to this FASTA file (for inspection in e.g. AliView)")
    parser.add_argument("--stats", metavar="FILE",
                        help="Save per-column stats to TSV file")
    args = parser.parse_args()

    # ── Load sequences ───────────────────────────────────────────────────────
    print(f"\n[1/4] Loading sequences from {args.input} ...")
    records = list(SeqIO.parse(args.input, "fasta"))
    if len(records) < 2:
        sys.exit("Need at least 2 sequences.")
    sequences = [str(r.seq).upper().replace("U", "T") for r in records]
    lengths = [len(s) for s in sequences]
    print(f"  Loaded {len(sequences)} sequences")
    print(f"  Lengths: min={min(lengths)}, max={max(lengths)}, "
          f"mean={int(np.mean(lengths))}")

    # ── MSA ──────────────────────────────────────────────────────────────────
    print("\n[2/4] Running multiple sequence alignment ...")
    if not args.no_mafft and mafft_available():
        print("  Using MAFFT (--auto) ...")
        aligned_seqs = run_mafft(records)
    else:
        if not args.no_mafft:
            print("  MAFFT not found — falling back to progressive pairwise MSA.")
        aligned_seqs = progressive_pairwise_msa(sequences)

    print(f"  MSA complete. Alignment length: {len(aligned_seqs[0])} columns")

    # Optionally save MSA
    if args.save_msa:
        msa_records = [
            SeqRecord(Seq(s), id=records[i].id, description="aligned")
            for i, s in enumerate(aligned_seqs)
        ]
        SeqIO.write(msa_records, args.save_msa, "fasta")
        print(f"  MSA saved to {args.save_msa}")

    # ── Consensus calling ────────────────────────────────────────────────────
    print("\n[3/4] Calling consensus ...")
    consensus_seq, stats = call_consensus(
        aligned_seqs,
        min_coverage=args.min_coverage,
        min_base_fraction=args.min_base_fraction,
        gap_fraction_threshold=args.gap_fraction_threshold,
    )

    called = [s for s in stats if s.get("base") not in (None, "-")]
    ambiguous = sum(1 for s in called if s.get("reason") == "ambiguous")
    print(f"  Consensus length : {len(consensus_seq)} bp")
    print(f"  Ambiguous bases  : {ambiguous} ({100*ambiguous/max(len(consensus_seq),1):.1f}%)")

    # Show where coverage drops
    low_cov = [s for s in stats if s["coverage"] < args.min_coverage]
    if low_cov:
        print(f"  Extension stopped at column {low_cov[0]['col']} "
              f"(coverage={low_cov[0]['coverage']} < {args.min_coverage})")
    else:
        print("  All columns had sufficient coverage throughout.")

    print_coverage_plot(stats)

    # ── Save outputs ─────────────────────────────────────────────────────────
    print(f"\n[4/4] Saving outputs ...")
    out_record = SeqRecord(
        Seq(consensus_seq),
        id="consensus",
        description=(f"len={len(consensus_seq)} "
                     f"min_cov={args.min_coverage} "
                     f"n_reads={len(sequences)}")
    )
    SeqIO.write([out_record], args.output, "fasta")
    print(f"  Consensus saved to {args.output}")

    if args.stats:
        import csv
        with open(args.stats, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["col","coverage","base","reason","fraction","top_base"],
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(stats)
        print(f"  Per-column stats saved to {args.stats}")

    print("\nDone ✓\n")


if __name__ == "__main__":
    main()
