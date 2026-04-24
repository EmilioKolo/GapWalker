#!/usr/bin/env python3

import argparse
from read_anchor_check import check_sequence_presence


def main():
    parser = argparse.ArgumentParser(
        description="Checks new consensus sequence for the other anchor and iterates."
    )

    parser.add_argument("-1", "--iterated-consensus", required=True,
                        help="Input fasta file with the iterated consensus sequence (output of cli_get_top_align.py).")
    parser.add_argument("-2", "--other-anchor-consensus", required=True,
                        help="Input fasta file with the other anchor consensus sequence.")
    parser.add_argument("-f", "--fastq-file", required=True,
                        help="Input fastq file with the sequenced reads.")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")
    parser.add_argument("-d", "--extension-direction", choices=['left', 'right'], default="right",
                        help="Direction for extension of the iterated consensus sequence (default: right). Must be either 'left' or 'right'.")

    args = parser.parse_args()

    print("Starting the pipeline...")

    ###

    # 1. Check if the iterated consensus sequence contains the other anchor

    result = check_sequence_presence(
        fasta_target=args.iterated_consensus,
        fasta_query=args.other_anchor_consensus,
        min_identity=95,
        min_coverage=95
    )

    print(f"Presence check result: {result}")

    # 2.A. If it does, stop the pipeline and assemble the final sequence
    # Check for identity and/or coverage over 90%
    if result["present"]:
        print("The iterated consensus sequence contains the other anchor. Stopping the pipeline and assembling the final sequence.")
        
        ### ADD CODE TO ASSEMBLE THE FINAL SEQUENCE HERE ###
        
        return
    elif result["identity_pct"] >= 90 or result["coverage_pct"] >= 90:
        print("The iterated consensus sequence has high identity and/or coverage, but is below threshold.")
        print(result)
    else:
        print("The iterated consensus sequence does not appear to contain the other anchor. Continuing with the pipeline.")

    # 2.B. If it doesn't, extract the 2kbp sequence from the extension side of the iterated consensus sequence

    # Extract the 2kbp sequence from the extension side of the iterated consensus sequence
    new_anchor = None
    if args.extension_direction == "right":
        new_anchor = result["aligned_target"][-2000:]
    else:
        new_anchor = result["aligned_target"][:2000]
    
    # Save the new anchor sequence to a fasta file
    new_anchor_fasta = f"{args.output_prefix}_new_anchor.fasta"
    with open(new_anchor_fasta, "w") as f:
        f.write(f">new_anchor\n{new_anchor}\n")

    # 3. Repeat the pipeline with the extracted sequence as the new anchor

    ### Do this with cli_extract_reads.py ###

    ###

    print('Pipeline finished!')


if __name__ == "__main__":
    main()
