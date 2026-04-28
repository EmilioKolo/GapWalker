#!/usr/bin/env python3

"""
Unified CLI entry point for the anchor-gap assembly pipeline.

Subcommands:
  extract-reads          Step 1  - align reads to anchors and extract shared reads
  align-reads            Step 2  - align every read to both anchors and record scores
  get-top-align          Step 3  - rank alignments and build a consensus sequence
  check-and-iter         Step 4  - check if the other anchor is present; 
                                   produce new anchor or stop
  concatenate-evaluate   Step 5  - concatenate consensus sequences and evaluate coverage
  run-all                        - Run the full pipeline end-to-end

Usage examples
--------------
# Run the full pipeline, stop after at most 5 iterations:
python pipeline_cli.py run-all \
    -1 anchor1.fasta -2 anchor2.fasta \
    -i reads.fastq -o ./results/run1 \
    --max-iter 5 --threads 4

# Run the full pipeline, loop indefinitely until the anchor is found:
python pipeline_cli.py run-all \
    -1 anchor1.fasta -2 anchor2.fasta \
    -i reads.fastq -o ./results/run1 \
    --max-iter -1 --threads 4

# Run only step 1:
python pipeline_cli.py extract-reads \
    -1 anchor1.fasta -2 anchor2.fasta \
    -i reads.fastq -o ./results/run1

# Dry-run the full pipeline (no commands executed):
python pipeline_cli.py run-all --dry-run \
    -1 anchor1.fasta -2 anchor2.fasta \
    -i reads.fastq -o ./results/run1
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from cli_extract_reads import main_pipeline as main1
from cli_align_reads_to_anchors import main_pipeline as main2
from cli_get_top_align import main_pipeline as main3
from cli_check_and_iter import main_pipeline as main4
from cli_concatenate_evaluate import main_pipeline as main5


# ===========================================================================
# Logging setup
# ===========================================================================

def setup_logging(log_file: str | None = None, verbose: bool = False) -> logging.Logger:
    """Configure root logger with console + optional file handler."""
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(level=level, format=fmt, datefmt=datefmt, handlers=handlers)
    return logging.getLogger(__name__)


logger = logging.getLogger(__name__)


# ===========================================================================
# Guard helpers
# ===========================================================================

class PipelineError(RuntimeError):
    """Raised when a pipeline step fails a pre- or post-condition check."""


def require_file_exists(path: str | Path, label: str = "") -> None:
    """Raise PipelineError if *path* does not exist."""
    p = Path(path)
    tag = f" ({label})" if label else ""
    if not p.exists():
        raise PipelineError(f"Required file not found{tag}: {p}")
    logger.debug("File exists check passed%s: %s", tag, p)


def require_file_nonempty(
    path: str | Path, 
    label: str = "", 
    fail_warn: bool = False
) -> None:
    """Raise PipelineError if *path* is missing or has zero bytes."""
    require_file_exists(path, label)
    p = Path(path)
    tag = f" ({label})" if label else ""
    if p.stat().st_size == 0:
        if fail_warn:
            logger.warning("File is empty%s: %s", tag, p)
            return
        else:
            raise PipelineError(f"File is empty{tag}: {p}")
    logger.debug("Non-empty check passed%s: %s", tag, p)


def ensure_output_dir(prefix: str) -> Path:
    """Create parent directory of *prefix* if it does not exist."""
    output_dir = Path(prefix.rsplit("/", 1)[0]) if "/" in prefix else Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


# ===========================================================================
# Per-step pre/post condition checks
# ===========================================================================

def check_inputs_extract_reads(args: argparse.Namespace) -> None:
    """Input guards for step 1 (extract-reads)."""
    require_file_exists(args.anchor1, "anchor1")
    require_file_exists(args.anchor2, "anchor2")
    require_file_nonempty(args.input_reads, "input-reads")


def check_outputs_extract_reads(prefix: str) -> None:
    """Output guards for step 1 - verify shared-read files were created."""
    require_file_nonempty(f"{prefix}_shared.txt", "shared reads", fail_warn=True)
    # _shared_other.txt may legitimately be empty in some runs; only check existence.
    require_file_nonempty(f"{prefix}_shared_other.txt", "shared reads (other)", fail_warn=True)


def check_inputs_align_reads(args: argparse.Namespace) -> None:
    """Input guards for step 2 (align-reads)."""
    require_file_exists(args.anchor1, "anchor1")
    require_file_exists(args.anchor2, "anchor2")
    require_file_nonempty(args.input_reads, "input-reads")


def check_outputs_align_reads(prefix: str) -> None:
    """Output guards for step 2 - mapping CSV must exist and be non-empty."""
    require_file_nonempty(f"{prefix}_read_id_mapping.csv", "read-id mapping")


def check_inputs_get_top_align(args: argparse.Namespace) -> None:
    """Input guards for step 3 (get-top-align)."""
    require_file_exists(args.input_folder, "input-folder")
    require_file_nonempty(args.read_id_mapping_file, "read-id-mapping-file")
    require_file_nonempty(args.fastq_file, "fastq-file")


def check_outputs_get_top_align(prefix: str, top_n: int) -> None:
    """Output guards for step 3 - consensus fasta must exist and be non-empty."""
    require_file_nonempty(
        f"{prefix}_top_{top_n}_alignments_with_seqs.csv",
        "top-alignments CSV",
    )
    require_file_nonempty(f"{prefix}_consensus_seq.fasta", "consensus sequence")


def check_inputs_check_and_iter(args: argparse.Namespace) -> None:
    """Input guards for step 4 (check-and-iter)."""
    require_file_nonempty(args.iterated_consensus, "iterated-consensus")
    require_file_nonempty(args.other_anchor_consensus, "other-anchor-consensus")


def check_outputs_check_and_iter(prefix: str) -> None:
    """Output guards for step 4 - results text file must exist."""
    require_file_exists(f"{prefix}_check_results.txt", "check results")


def check_inputs_concatenate_evaluate(args: argparse.Namespace) -> None:
    """Input guards for step 5 (concatenate-evaluate)."""
    for fasta in args.consensus_sequences:
        require_file_nonempty(fasta, f"consensus sequence: {fasta}")
    require_file_nonempty(args.fastq, "fastq")


def check_outputs_concatenate_evaluate(prefix: str) -> None:
    """Output guards for step 5 - full consensus fasta must exist."""
    require_file_nonempty(f"{prefix}_full_consensus.fasta", "full consensus")


def _extract_query_end(check_file: str) -> int:
    """Read check_file and return the integer value of 'Query end'."""
    with open(check_file) as f:
        for line in f:
            if line.startswith("Query end:"):
                return int(line.split(":", 1)[1].strip())
    raise PipelineError(f"'Query end' not found in {check_file}")


# ===========================================================================
# Step runners  (thin wrappers around the original main* functions)
# ===========================================================================

def run_extract_reads(args: argparse.Namespace, dry_run: bool = False) -> None:
    """Step 1 - extract_reads."""
    logger.info("=== Step 1: extract-reads ===")
    check_inputs_extract_reads(args)
    ensure_output_dir(args.output_prefix)

    if dry_run:
        logger.info("[DRY-RUN] Would run extract-reads with prefix=%s", args.output_prefix)
        return

    # ---- actual work ----
    main1_args = argparse.Namespace(
        anchor1=args.anchor1,
        anchor2=args.anchor2,
        input_reads=args.input_reads,
        output_prefix=args.output_prefix,
        threads=args.threads,
        skip_align=getattr(args, "skip_align", False),
        skip_extract=getattr(args, "skip_extract", False),
    )
    main1(main1_args)
    # ---- post-checks ----
    check_outputs_extract_reads(args.output_prefix)
    logger.info("Step 1 finished successfully.")


def run_align_reads(args: argparse.Namespace, dry_run: bool = False) -> None:
    """Step 2 - align_reads_to_anchors."""
    logger.info("=== Step 2: align-reads ===")
    check_inputs_align_reads(args)
    ensure_output_dir(args.output_prefix)

    if dry_run:
        logger.info("[DRY-RUN] Would run align-reads with prefix=%s", args.output_prefix)
        return

    main2_args = argparse.Namespace(
        anchor1=args.anchor1,
        anchor2=args.anchor2,
        input_reads=args.input_reads,
        output_prefix=args.output_prefix,
    )
    main2(main2_args)
    check_outputs_align_reads(args.output_prefix)
    logger.info("Step 2 finished successfully.")


def run_get_top_align(args: argparse.Namespace, dry_run: bool = False) -> None:
    """Step 3 - get_top_align."""
    logger.info("=== Step 3: get-top-align ===")
    check_inputs_get_top_align(args)
    ensure_output_dir(args.output_prefix)

    if dry_run:
        logger.info("[DRY-RUN] Would run get-top-align with prefix=%s", args.output_prefix)
        return

    main3_args = argparse.Namespace(
        input_folder=args.input_folder,
        top_n=args.top_n,
        read_id_mapping_file=args.read_id_mapping_file,
        fastq_file=args.fastq_file,
        output_prefix=args.output_prefix,
        priority_anchor=args.priority_anchor,
    )
    main3(main3_args)
    check_outputs_get_top_align(args.output_prefix, args.top_n)
    logger.info("Step 3 finished successfully.")


def run_check_and_iter(args: argparse.Namespace, dry_run: bool = False) -> bool:
    """
    Step 4 - check_and_iter.

    Returns
    -------
    True  - other anchor found; pipeline should proceed to step 5.
    False - other anchor not found; caller should loop back to step 1.
    """
    logger.info("=== Step 4: check-and-iter ===")
    check_inputs_check_and_iter(args)
    ensure_output_dir(args.output_prefix)

    if dry_run:
        logger.info("[DRY-RUN] Would run check-and-iter with prefix=%s", args.output_prefix)
        return False  # Assume loop continues in dry-run

    main4_args = argparse.Namespace(
        iterated_consensus=args.iterated_consensus,
        other_anchor_consensus=args.other_anchor_consensus,
        output_prefix=args.output_prefix,
        extension_direction=args.extension_direction,
    )
    found = main4(main4_args)  # main4 should return True when anchor is found
    check_outputs_check_and_iter(args.output_prefix)

    if found:
        logger.info("Other anchor found. Proceeding to concatenate-evaluate.")
    else:
        logger.info("Other anchor not found. Preparing next iteration.")
    return bool(found)


def run_concatenate_evaluate(args: argparse.Namespace, dry_run: bool = False) -> None:
    """Step 5 - concatenate_evaluate."""
    logger.info("=== Step 5: concatenate-evaluate ===")
    check_inputs_concatenate_evaluate(args)
    ensure_output_dir(args.output_prefix)

    if dry_run:
        logger.info("[DRY-RUN] Would run concatenate-evaluate with prefix=%s", args.output_prefix)
        return

    main5_args = argparse.Namespace(
        consensus_sequences=args.consensus_sequences,
        output_prefix=args.output_prefix,
        final_pos=args.final_pos,
        fastq=args.fastq,
        threads=args.threads,
    )
    main5(main5_args)
    check_outputs_concatenate_evaluate(args.output_prefix)
    logger.info("Step 5 finished successfully.")


# ===========================================================================
# run-all orchestrator
# ===========================================================================

def run_all(args: argparse.Namespace) -> None:
    """
    Orchestrate the full pipeline:
      1. extract-reads
      2. align-reads
      3. get-top-align
      loop (up to --max-iter times, or indefinitely if max-iter < 0):
        4. check-and-iter
           ├─ found  → break → 5. concatenate-evaluate
           └─ not found → update anchor → back to 1
    """
    dry_run: bool = args.dry_run
    max_iter: int = args.max_iter
    infinite: bool = max_iter < 0

    if infinite:
        logger.info(
            "Starting full pipeline. max_iter=unlimited  dry_run=%s  "
            "(pass Ctrl-C to abort)",
            dry_run,
        )
    else:
        logger.info("Starting full pipeline. max_iter=%d  dry_run=%s", max_iter, dry_run)

    # Determine extension direction from priority anchor (inverse relationship)
    extension_direction = "right" if args.priority_anchor == "left" else "left"

    # Track consensus sequences produced across iterations (fed into step 5).
    consensus_files: list[str] = []

    # The "current" anchor on the extending side changes after each iteration.
    current_anchor_ext = args.anchor1  # or anchor2, depending on direction

    iteration = 0
    anchor_found = False

    while True:
        iteration += 1

        # Enforce max-iter cap when not running indefinitely
        if not infinite and iteration > max_iter:
            logger.warning(
                "Reached maximum number of iterations (%d) without finding the other anchor. "
                "Pass --max-iter -1 to iterate indefinitely, or review your anchors.",
                max_iter,
            )
            sys.exit(2)

        iter_label = f"{iteration}/∞" if infinite else f"{iteration}/{max_iter}"
        logger.info("--- Iteration %s ---", iter_label)
        iter_dir = os.path.join(args.output_prefix, f"iter{iteration}")

        # Step‑specific output prefixes inside iteration folder
        step1_prefix = os.path.join(iter_dir, "01-extract-reads", "output")
        step2_prefix = os.path.join(iter_dir, "02-align-reads", "read")
        step3_prefix = os.path.join(iter_dir, "03-get-top-align", "output")
        step4_prefix = os.path.join(iter_dir, "04-check-and-iter", "output")

        # ------------------------------------------------------------------
        # Step 1 - extract reads
        # ------------------------------------------------------------------
        step1_args = argparse.Namespace(
            anchor1=current_anchor_ext,
            anchor2=args.anchor2,
            input_reads=args.input_reads,
            output_prefix=step1_prefix,
            threads=args.threads,
            skip_align=False,
            skip_extract=False,
        )
        try:
            run_extract_reads(step1_args, dry_run=dry_run)
        except PipelineError as exc:
            logger.error("Step 1 failed on iteration %d: %s", iteration, exc)
            sys.exit(1)

        # ------------------------------------------------------------------
        # Step 2 - align reads to both anchors
        # ------------------------------------------------------------------
        step1_aligned_anchor1 = f"{step1_prefix}_anchor1.other.fastq"
        step2_args = argparse.Namespace(
            anchor1=current_anchor_ext,
            anchor2=args.anchor2,
            input_reads=step1_aligned_anchor1,
            output_prefix=step2_prefix,
        )
        try:
            run_align_reads(step2_args, dry_run=dry_run)
        except PipelineError as exc:
            logger.error("Step 2 failed on iteration %d: %s", iteration, exc)
            sys.exit(1)

        # ------------------------------------------------------------------
        # Step 3 - rank alignments and build consensus
        # ------------------------------------------------------------------
        step3_input_folder = os.path.dirname(step2_prefix)   # "02-align-reads"
        step3_args = argparse.Namespace(
            input_folder=step3_input_folder,
            top_n=args.top_n,
            read_id_mapping_file=f"{step2_prefix}_read_id_mapping.csv",
            fastq_file=args.input_reads,
            output_prefix=step3_prefix,
            priority_anchor=args.priority_anchor,
        )
        try:
            run_get_top_align(step3_args, dry_run=dry_run)
        except PipelineError as exc:
            logger.error("Step 3 failed on iteration %d: %s", iteration, exc)
            sys.exit(1)

        iter_consensus = f"{step3_prefix}_consensus_seq.fasta"
        consensus_files.append(iter_consensus)

        # ------------------------------------------------------------------
        # Step 4 - check if other anchor is present; decide to loop or stop
        # ------------------------------------------------------------------
        step4_args = argparse.Namespace(
            iterated_consensus=iter_consensus,
            other_anchor_consensus=args.anchor2,
            output_prefix=step4_prefix,
            extension_direction=extension_direction,
        )
        try:
            anchor_found = run_check_and_iter(step4_args, dry_run=dry_run)
        except PipelineError as exc:
            logger.error("Step 4 failed on iteration %d: %s", iteration, exc)
            sys.exit(1)

        if anchor_found:
            logger.info("Anchor found after %d iteration(s). Moving to final assembly.", iteration)
            break

        # Prepare the new anchor for the next iteration
        new_anchor_path = f"{step4_prefix}_new_anchor.fasta"
        if not dry_run:
            try:
                require_file_nonempty(new_anchor_path, "new anchor for next iteration")
            except PipelineError as exc:
                logger.error("Cannot continue iteration: %s", exc)
                sys.exit(1)
        current_anchor_ext = new_anchor_path
        logger.info("New anchor for next iteration: %s", new_anchor_path)

    # ----------------------------------------------------------------------
    # Step 5 - concatenate consensus sequences and evaluate coverage
    # ----------------------------------------------------------------------
    check_file = f"{step4_prefix}_check_results.txt"   # last iteration’s step4
    if dry_run:
        final_pos = 0
    else:
        final_pos = _extract_query_end(check_file)

    step5_output_prefix = os.path.join(args.output_prefix, "05-concatenate-evaluate", "output")
    step5_args = argparse.Namespace(
        consensus_sequences=consensus_files,
        output_prefix=step5_output_prefix,
        final_pos=final_pos,
        fastq=args.input_reads,
        threads=args.threads,
    )
    try:
        run_concatenate_evaluate(step5_args, dry_run=dry_run)
    except PipelineError as exc:
        logger.error("Step 5 failed: %s", exc)
        sys.exit(1)

    logger.info("Full pipeline completed successfully. Output prefix: %s", args.output_prefix)


# ===========================================================================
# Argument parsers
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="pipeline",
        description="Anchor-gap assembly pipeline - unified CLI.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Global flags available to every subcommand
    parser.add_argument(
        "--log-file",
        default=None,
        metavar="PATH",
        help="Write log output to this file in addition to stdout.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be done without executing any pipeline steps.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # ------------------------------------------------------------------
    # extract-reads (step 1)
    # ------------------------------------------------------------------
    p1 = subparsers.add_parser(
        "extract-reads",
        help="Align reads to anchors and extract shared reads.",
        description="Step 1 - align reads to both anchors, extract mapped reads, find shared reads.",
    )
    _add_common_io(p1)
    p1.add_argument("--skip-align", action="store_true", help="Skip the minimap2 alignment step.")
    p1.add_argument("--skip-extract", action="store_true", help="Skip the read extraction step.")
    p1.set_defaults(func=lambda a: run_extract_reads(a, dry_run=a.dry_run))

    # ------------------------------------------------------------------
    # align-reads (step 2)
    # ------------------------------------------------------------------
    p2 = subparsers.add_parser(
        "align-reads",
        help="Align every read to both anchors and record alignment scores.",
        description="Step 2 - dual alignment of every read against both anchors.",
    )
    p2.add_argument("-1", "--anchor1", required=True, help="First anchor (left side).")
    p2.add_argument("-2", "--anchor2", required=True, help="Second anchor (right side).")
    p2.add_argument("-i", "--input-reads", required=True, help="Input FASTQ file.")
    p2.add_argument(
        "-o", "--output-prefix", default="./output",
        help="Output prefix (default: ./output).",
    )
    p2.set_defaults(func=lambda a: run_align_reads(a, dry_run=a.dry_run))

    # ------------------------------------------------------------------
    # get-top-align (step 3)
    # ------------------------------------------------------------------
    p3 = subparsers.add_parser(
        "get-top-align",
        help="Rank alignments and build a consensus sequence.",
        description="Step 3 - select top-N alignments, extract sequences, build consensus.",
    )
    p3.add_argument("-i", "--input-folder", required=True, help="Folder with alignment files.")
    p3.add_argument("-n", "--top-n", type=int, default=10, help="Number of top alignments (default: 10).")
    p3.add_argument("-m", "--read-id-mapping-file", required=True, help="CSV with read_cont→read_id mapping.")
    p3.add_argument("-f", "--fastq-file", required=True, help="Input FASTQ file.")
    p3.add_argument("-o", "--output-prefix", default="./output", help="Output prefix.")
    p3.add_argument(
        "-p", "--priority-anchor", choices=["left", "right"], default="left",
        help="Priority anchor for sequence cutting (default: left).",
    )
    p3.set_defaults(func=lambda a: run_get_top_align(a, dry_run=a.dry_run))

    # ------------------------------------------------------------------
    # check-and-iter (step 4)
    # ------------------------------------------------------------------
    p4 = subparsers.add_parser(
        "check-and-iter",
        help="Check if the other anchor is present; emit new anchor or signal completion.",
        description="Step 4 - check iterated consensus for the other anchor.",
    )
    p4.add_argument("-1", "--iterated-consensus", required=True, help="Iterated consensus FASTA.")
    p4.add_argument("-2", "--other-anchor-consensus", required=True, help="Other anchor consensus FASTA.")
    p4.add_argument("-o", "--output-prefix", default="./output", help="Output prefix.")
    p4.add_argument(
        "-d", "--extension-direction", choices=["left", "right"], default="right",
        help="Direction to extend the consensus (default: right).",
    )
    p4.set_defaults(func=lambda a: run_check_and_iter(a, dry_run=a.dry_run))

    # ------------------------------------------------------------------
    # concatenate-evaluate (step 5)
    # ------------------------------------------------------------------
    p5 = subparsers.add_parser(
        "concatenate-evaluate",
        help="Concatenate consensus sequences and evaluate read coverage.",
        description="Step 5 - assemble full consensus and measure depth.",
    )
    p5.add_argument(
        "-1", "--consensus-sequences", nargs="+", required=True,
        help="FASTA files with consensus sequences to concatenate.",
    )
    p5.add_argument("-o", "--output-prefix", default="./output", help="Output prefix.")
    p5.add_argument("--final-pos", type=int, required=True, help="Final trim position in last consensus.")
    p5.add_argument("-f", "--fastq", required=True, help="FASTQ file for coverage evaluation.")
    p5.add_argument("-t", "--threads", type=int, default=1, help="Alignment threads (default: 1).")
    p5.set_defaults(func=lambda a: run_concatenate_evaluate(a, dry_run=a.dry_run))

    # ------------------------------------------------------------------
    # run-all  (full pipeline)
    # ------------------------------------------------------------------
    pa = subparsers.add_parser(
        "run-all",
        help="Run the complete pipeline end-to-end.",
        description=(
            "Run all 5 steps in order. "
            "Steps 1-4 iterate until the other anchor is found or --max-iter is reached."
        ),
    )
    _add_common_io(pa)
    pa.add_argument(
        "--max-iter", type=int, default=10,
        help=(
            "Maximum number of extension iterations before giving up (default: 10). "
            "Pass any negative value (e.g. -1) to iterate indefinitely until the "
            "anchor is found."
        ),
    )
    pa.add_argument(
        "--top-n", type=int, default=10,
        help="Top-N alignments passed to get-top-align (default: 10).",
    )
    pa.add_argument(
        "-p", "--priority-anchor", choices=["left", "right"], default="left",
        help="Priority anchor for sequence cutting (default: left)."
             "Extension direction will be the opposite side automatically.",
    )

    pa.set_defaults(func=run_all)

    return parser


def _add_common_io(p: argparse.ArgumentParser) -> None:
    """Add the I/O arguments shared by steps 1, 2, and run-all."""
    p.add_argument("-1", "--anchor1", required=True, help="First anchor (left side).")
    p.add_argument("-2", "--anchor2", required=True, help="Second anchor (right side).")
    p.add_argument("-i", "--input-reads", required=True, help="Input FASTQ (or .fastq.gz) file.")
    p.add_argument("-o", "--output-prefix", default="./output", help="Output prefix (default: ./output).")
    p.add_argument("-t", "--threads", type=int, default=1, help="Threads for alignment (default: 1).")


# ===========================================================================
# Entry point
# ===========================================================================

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Set up logging (global flags are parsed on the root namespace)
    setup_logging(log_file=args.log_file, verbose=args.verbose)

    if args.dry_run:
        logger.info("*** DRY-RUN MODE - no commands will be executed ***")

    try:
        args.func(args)
    except PipelineError as exc:
        logger.error("Pipeline error: %s", exc)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Pipeline interrupted by user.")
        sys.exit(130)


if __name__ == "__main__":
    main()
