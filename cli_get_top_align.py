#!/usr/bin/env python3

import argparse
from pathlib import Path
from process_alignments import *
from consensus_builder import run_consensus_builder


def main():
    parser = argparse.ArgumentParser(
        description="Order the alignment scores and select the top alignments."
    )
    parser.add_argument("-i", "--input-folder", required=True,
                        help="Input folder with the alignment files (output of cli_align_reads_to_anchors.py).")
    parser.add_argument("-n", "--top-n", type=int, default=10,
                        help="Number of top alignments to select for each read (default: 10).")
    parser.add_argument("-m", "--read-id-mapping-file", required=True,
                        help="CSV file with the mapping of read_cont to read_id (output of cli_align_reads_to_anchors.py).")
    parser.add_argument("-f", "--fastq-file", required=True,
                        help="Input fastq file with the sequenced reads.")
    parser.add_argument("-o", "--output-prefix", default="./output",
                        help="Output prefix. Assumes it's a folder if it ends with / (default: ./output)")
    parser.add_argument("-p", "--priority-anchor", choices=['left', 'right'], default="left",
                        help="Priority anchor for cutting the sequences (default: left). Must be either 'left' or 'right'.")

    args = parser.parse_args()

    print("Starting the pipeline...")

    # Make sure output directory exists
    output_dir = Path(args.output_prefix.rsplit("/", 1)[0])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    ###

    df_top_n = read_alignment_files(
        Path(args.input_folder),
        int(args.top_n),
        Path(args.read_id_mapping_file),
        priority_anchor=args.priority_anchor
    )

    # Get the sequences from the top n alignments
    df_seqs = get_seq_from_aln(
        df_top_n,
        Path(args.fastq_file),
        priority_anchor=args.priority_anchor
    )
    # Save the top n alignments with sequences to a csv file
    df_seqs.to_csv(f"{args.output_prefix}_top_{args.top_n}_alignments_with_seqs.csv",
                   index=False)

    # Define a consensus sequence for the top n alignments
    ext_direction = 'right' if args.priority_anchor == 'left' else 'left'
    # Save the sequences to a fasta file for consensus building
    with open(f"{args.output_prefix}_top_{args.top_n}_sequences.fasta", "w") as f:
        for index, row in df_seqs.iterrows():
            f.write(f">{row['read_id']}\n{row['sequence']}\n")
    # Run the consensus builder on the top n sequences
    consensus_seq, col_stats = run_consensus_builder(
        input_file=f"{args.output_prefix}_top_{args.top_n}_sequences.fasta",
        output_file=f"{args.output_prefix}_consensus_seq.fasta",
        min_base_fraction=0.4,
        gap_fraction_threshold=0.5,
        #no_mafft=True,
        save_msa=f"{args.output_prefix}_msa.fasta",
        stats=f"{args.output_prefix}_consensus_stats.tsv"
    )

    # Save the consensus sequence to a fasta file
    with open(f"{args.output_prefix}_consensus_seq.fasta", "w") as f:
        f.write(f">consensus_seq\n{consensus_seq}\n")
    
    ###

    print('Pipeline finished!')


if __name__ == "__main__":
    main()
