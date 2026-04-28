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
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")
    parser.add_argument("-d", "--extension-direction", choices=['left', 'right'], default="right",
                        help="Direction for extension of the iterated consensus sequence (default: right). Must be either 'left' or 'right'.")

    args_main = parser.parse_args()

    anchor_found = main_pipeline(args_main)

    if anchor_found:
        print("The iterated consensus sequence contains the other anchor. Stopping the pipeline and assembling the final sequence.")
        
        print("[ACTION NEEDED] Run cli_concatenate_evaluate.py with all the generated consensus sequences.")
        ### Do this with cli_concatenate_evaluate.py ###
    else: 
        # 3. Repeat the pipeline with the extracted sequence as the new anchor
        print("The iterated consensus sequence does not contain the other anchor.")

        print("[ACTION NEEDED] Run cli_check_and_iter.py again with the new anchor fasta file as the other anchor consensus.")
        ### Do this with cli_extract_reads.py ###

    print('Pipeline finished!')


def main_pipeline(args: argparse.Namespace):
    """
    For calling from other scripts, e.g. cli_main.py.
    """
    print("Starting the pipeline for check_and_iter...")
    # 1. Check if the iterated consensus sequence contains the other anchor
    result = check_sequence_presence(
        fasta_target=args.iterated_consensus,
        fasta_query=args.other_anchor_consensus,
        min_identity=95,
        min_coverage=95
    )

    # Save results to a file
    with open(f"{args.output_prefix}_check_results.txt", "w") as f:
        f.write(f"Present: {result['present']}\n")
        f.write(f"Score: {result['score']}\n")
        f.write(f"Identity %: {result['identity_pct']}\n")
        f.write(f"Coverage %: {result['coverage_pct']}\n")
        f.write(f"Target start: {result['target_start']}\n")
        f.write(f"Target end: {result['target_end']}\n")
        f.write(f"Query start: {result['query_start']}\n")
        f.write(f"Query end: {result['query_end']}\n")
        f.write(f"Number of matches: {result['n_matches']}\n")
        f.write(f"Number of mismatches: {result['n_mismatches']}\n")
        f.write(f"Number of gaps: {result['n_gaps']}\n")
        f.write(f"Query length: {result['query_len']}\n")
        f.write(f"Target length: {result['target_len']}\n")

    # Check for identity and coverage over 95%
    if result["present"]:
        return True
    # Check for identity and/or coverage over 90%
    elif result["identity_pct"] >= 90 or result["coverage_pct"] >= 90:
        print("The iterated consensus sequence has high identity and/or coverage, but is below threshold.")
        print(result)
    else:
        print("The iterated consensus sequence does not appear to contain the other anchor. Continuing with the pipeline.")

    # 2.B. If it doesn't, extract the 2kbp sequence from the extension side of the iterated consensus sequence

    # Get the full new sequence from the iterated consensus fasta file
    with open(args.iterated_consensus, "r") as f:
        lines = f.readlines()
        new_sequence = "".join(lines[1:]).replace("\n", "")
    # Extract the 2kbp sequence from the extension side of the iterated consensus sequence
    new_anchor = None
    if args.extension_direction == "right":
        new_anchor = new_sequence[-2000:]
    else:
        new_anchor = new_sequence[:2000]
    
    # Save the new anchor sequence to a fasta file
    new_anchor_fasta = f"{args.output_prefix}_new_anchor.fasta"
    with open(new_anchor_fasta, "w") as f:
        f.write(f">new_anchor\n{new_anchor}\n")

    return False


if __name__ == "__main__":
    main()
