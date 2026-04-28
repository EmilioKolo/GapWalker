#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
from read_anchor_check import run_dual_alignment


def main():
    parser = argparse.ArgumentParser(
        description="Align reads to both anchors and see where they align."
    )
    parser.add_argument("-1", "--anchor1", required=True, 
                        help="First anchor (left side).")
    parser.add_argument("-2", "--anchor2", required=True, 
                        help="Second anchor (right side).")
    parser.add_argument("-i", "--input-reads", required=True,
                        help="Input file with the sequenced reads (only .fastq format)")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")

    args = parser.parse_args()

    print("Starting the pipeline for align_reads_to_anchors...")
    
    ###

    # Make sure output directory exists
    output_dir = Path(args.output_prefix.rsplit("/", 1)[0])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Go through each of the reads

    # Start counting reads from 1 for better readability
    START_CONT = 1
    read_cont = START_CONT
    cont = 0
    # Start a dictionary to store read_cont to read_id mapping
    read_id_dict = {}
    # Open the fastq file and go through IDs
    with Path(args.input_reads).open() as f:
        for line in f:
            if line.startswith('@') and cont % 4 == 0:
                read_id = line.strip()[1:]
                print(f"Read ID: {read_id}")

                read_out_pref = f"{args.output_prefix}_{read_cont}"
                read_id_dict[read_cont] = read_id

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
                cont += 1
                read_cont += 1
            elif cont % 4 != 0:
                cont += 1
    
    # Save the read_id_dict as a csv file for later use
    df_read_id = pd.DataFrame.from_dict(read_id_dict, orient='index', columns=['read_id'])
    df_read_id.to_csv(f"{args.output_prefix}_read_id_mapping.csv", index_label='read_cont')
    
    ###

    print('Pipeline finished!')


if __name__ == "__main__":
    main()
