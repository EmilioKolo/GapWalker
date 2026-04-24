#!/usr/bin/env python3

import argparse
from pathlib import Path
from alignment import *


def main():
    parser = argparse.ArgumentParser(
        description="Extract reads that align to the anchors and extend the sequence between them."
    )
    parser.add_argument("-1", "--anchor1", required=True, 
                        help="First anchor from which the sequence will be extended.")
    parser.add_argument("-2", "--anchor2", required=True, 
                        help="Second anchor used as target for the extension.")
    parser.add_argument("-i", "--input-reads", required=True,
                        help="Input file with the sequenced reads (.fastq or .fastq.gz format)")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for the different steps (default: 1).")
    parser.add_argument("--skip-align", action="store_true",
                        help="Add to skip alignment step")
    parser.add_argument("--skip-extract", action="store_true",
                        help="Add to skip read extraction step.")

    args = parser.parse_args()

    print("Starting the pipeline...")
    
    ###

    # Make sure output directory exists
    output_dir = Path(args.output_prefix.rsplit("/", 1)[0])
    output_dir.mkdir(parents=True, exist_ok=True)

    ## 1. Align reads to anchors
    if not args.skip_align:
        print("Aligning reads to anchors...")
        align_index_minimap2(
            args.input_reads,
            args.anchor1,
            f"{args.output_prefix}_anchor1.bam",
            int(args.threads)
        )
        align_index_minimap2(
            args.input_reads,
            args.anchor2,
            f"{args.output_prefix}_anchor2.bam",
            int(args.threads)
        )
    else:
        print("Skipping alignment step...")
    
    ## 2. Extract reads that align to the anchors
    if not args.skip_extract:
        print("Extracting reads that align to the anchors...")
        extract_aligned_reads(
            f"{args.output_prefix}_anchor1.sorted.bam",
            f"{args.output_prefix}_anchor1.fastq"
        )
        extract_aligned_reads(
            f"{args.output_prefix}_anchor2.sorted.bam",
            f"{args.output_prefix}_anchor2.fastq"
        )
    else:
        print("Skipping read extraction step...")
    
    ## 3. Extract shared reads between the two fastq files
    print("Extracting shared reads between the two anchors...")
    get_shared_reads(
        f"{args.output_prefix}_anchor1.fastq",
        f"{args.output_prefix}_anchor2.fastq",
        f"{args.output_prefix}_shared.txt"
    )
    get_shared_reads(
        f"{args.output_prefix}_anchor1.other.fastq",
        f"{args.output_prefix}_anchor2.other.fastq",
        f"{args.output_prefix}_shared_other.txt"
    )

    print('Pipeline finished!')


if __name__ == "__main__":
    main()
