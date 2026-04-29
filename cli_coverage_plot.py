#!/usr/bin/env python3

import argparse
from alignment import plot_coverage
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Order the alignment scores and select the top alignments."
    )
    parser.add_argument("-i", "--coverage-file", required=True,
                        help="Input file with coverage data (output of cli_concatenate_evaluate.py).")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")
    parser.add_argument("-m", "--max-coverage", type=int, default=None,
                        help="Maximum coverage to plot (default: None, i.e. no limit).")

    args_main = parser.parse_args()

    main_pipeline(args_main)

    print('Pipeline finished!')


def main_pipeline(args: argparse.Namespace):
    """
    For calling from other scripts, e.g. cli_main.py.
    """
    print("Starting the pipeline for coverage_plot...")
    # Make sure output directory exists
    output_dir = Path(args.output_prefix.rsplit("/", 1)[0])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    ###

    plot_coverage(
        args.coverage_file, 
        f"{args.output_prefix}_coverage_plot.png", 
        max_coverage=args.max_coverage
    )

    ###


if __name__ == "__main__":
    main()

