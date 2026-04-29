#!/usr/bin/env python3

import argparse
from alignment import align_index_minimap2


def main():
    parser = argparse.ArgumentParser(
        description="Concatenates the different assembled sequences and evaluates alignment."
    )

    parser.add_argument("-1", "--consensus-sequences", nargs="+", required=True,
                        help="Input fasta files with the consensus sequences to concatenate and evaluate.")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output).")
    parser.add_argument("--final-pos", type=int, required=True,
                        help="Final position for the anchor on the other side of the gap.")
    parser.add_argument("-f", "--fastq", required=True,
                        help="File with reads to evaluate the coverage of the full consensus sequence.")
    parser.add_argument("-t", "--threads", type=int, default=1, 
                        help="Number of threads to use for alignment (default: 1).")

    args_main = parser.parse_args()

    main_pipeline(args_main)

    print('Pipeline finished!')


def main_pipeline(args: argparse.Namespace):
    """
    For calling from other scripts, e.g. cli_main.py.
    """
    print("Starting the pipeline for concatenate_evaluate...")
    # Open and concatenate the different consensus sequences (remove the last part of the last consensus seq)
    full_sequence:str = ''
    for idx, fasta_file in enumerate(args.consensus_sequences):
        # Start a string with the current sequence
        curr_sequence:str = ''
        # Start a boolean to check id and sequence
        id_found = False
        with open(fasta_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(">"):
                    if not id_found:
                        id_found = True
                    else:
                        print(f"WARNING: Two successive sequence IDs found in {fasta_file}. Check the file format.")
                else:
                    id_found = False
                    curr_sequence += str(line.strip())
        # Check if the current sequence is the last one
        if idx < len(args.consensus_sequences) - 1:
            full_sequence += curr_sequence
        else:
            # Remove the last part of the last consensus sequence
            full_sequence += curr_sequence[:args.final_pos]

    # Save the full consensus sequence to a fasta file
    full_consensus_fasta = f"{args.output_prefix}_full_consensus.fasta"
    with open(full_consensus_fasta, "w") as f:
        f.write(f">full_consensus\n")
        f.write(f"{full_sequence}\n")

    # Align reads to new consensus sequence and visualize depth
    align_index_minimap2(
        args.fastq, 
        full_consensus_fasta, 
        f"{args.output_prefix}_full_consensus.bam", 
        args.threads,
        coverage_plots=True
    )

    ###


if __name__ == "__main__":
    main()
