#!/usr/bin/env python3

import argparse
from pathlib import Path
from read_anchor_check import run_dual_alignment


def main():
    parser = argparse.ArgumentParser(
        description="See the reads that align to both anchors and where they align."
    )
    parser.add_argument("-1", "--anchor1", required=True, 
                        help="First anchor (left side).")
    parser.add_argument("-2", "--anchor2", required=True, 
                        help="Second anchor (right side).")
    parser.add_argument("-i", "--input-reads", required=True,
                        help="Input file with the sequenced reads (.fastq or .fastq.gz format)")
    parser.add_argument("-r", "--read-list", required=True,
                        help="File with a list of read IDs that align to both anchors.")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for the different steps (default: 1).")

    args = parser.parse_args()

    print("Starting the pipeline for see_shared_reads...")
    
    ###

    # Make sure output directory exists
    output_dir = Path(args.output_prefix.rsplit("/", 1)[0])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define path to the list of reads that align to both anchors
    shared_reads_file:Path = Path(args.read_list)

    # Start counting reads from 1 for better readability
    read_cont = 1
    # Open the file and go through lines
    with shared_reads_file.open() as f:
        for line in f:
            read_id = line.strip()
            print(f"Read ID: {read_id}")

            read_out_pref = f"{args.output_prefix}_{read_cont}"

            # Run alignment for this read against both anchors
            _ = run_dual_alignment(
                args.anchor1,
                args.anchor2,
                args.input_reads,
                read_id,
                read_out_pref,
                label_a='Anchor 1 (left)',
                label_b='Anchor 2 (right)'
            )
            read_cont += 1
    ###

    print('Pipeline finished!')


if __name__ == "__main__":
    main()
