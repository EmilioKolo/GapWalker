#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import subprocess


def align_index_minimap2(
    fastq, 
    ref_seq, 
    bam_out:str, 
    threads:int, 
    coverage_plots:bool=True
):
    """
    Performs alignment with minimap2 and creates an index for the 
    resulting BAM file.
    """
    # Check that fasta file index exists, if not create it
    if not Path(ref_seq + ".fai").exists():
        index_fasta(ref_seq)
    # Run minimap2 to align reads to the reference sequence
    run_minimap2(fastq, bam_out, threads, ref_seq)
    # Define sorted bam
    bam_sorted = bam_out.replace(".bam", ".sorted.bam")
    # Sort the BAM file
    cmd_sort = ["samtools", "sort", bam_out, "-o", bam_sorted]
    print(f"Running command: {' '.join(cmd_sort)}")
    subprocess.run(" ".join(cmd_sort), shell=True, check=True)
    # Index the BAM file
    cmd_index = ["samtools", "index", bam_sorted]
    print(f"Running command: {' '.join(cmd_index)}")
    subprocess.run(" ".join(cmd_index), shell=True, check=True)
    if coverage_plots:
        # Generate coverage statistics
        coverage_file = bam_sorted.replace(".sorted.bam", "_coverage.txt")
        generate_coverage_stats(bam_sorted, coverage_file)
        # Plot coverage
        coverage_plot_file = bam_sorted.replace(".sorted.bam", "_coverage_plot.png")
        plot_coverage(coverage_file, coverage_plot_file)


def extract_aligned_reads(bam_file, output_fastq, aln_len:int=1500):
    """
    Extract reads that align to the anchors from the BAM file 
    and save them in FASTQ format.
    """
    # Open input BAM
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    # Create temporary BAM file to store filtered reads
    filtered_bam = bam_file.replace(".bam", ".filtered.bam")
    # Open output BAM
    outfile = pysam.AlignmentFile(filtered_bam, "wb", template=bamfile)
    for read in bamfile.fetch():
        try:
            ref_len = int(str(read.reference_length))
        except TypeError:
            ref_len = 0
        # excluding soft-clipped bases but including gaps/introns
        if ref_len > aln_len:
            outfile.write(read)
    # Close files
    bamfile.close()
    outfile.close()

    cmd = [
        "samtools", "fastq", 
        "-F", "4",     # Filter out unmapped reads
        "-F", "0x900", # Filter secondary and supplementary alignments
        "-F", "0x400", # Filter duplicate reads
        "-0", output_fastq.replace(".fastq", ".other.fastq"), 
        "-o", output_fastq, 
        filtered_bam
    ]
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(" ".join(cmd), shell=True, check=True)


def generate_coverage_stats(bam_file, output_file):
    """
    Generate coverage statistics from the BAM file using samtools depth.
    """
    cmd = [
        "samtools", "depth",
        "-a", bam_file,
        ">", output_file
    ]
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(" ".join(cmd), shell=True, check=True)


def get_shared_reads(fastq1, fastq2, output_file):
    """
    Get the shared reads between two FASTQ files and save them in a new FASTQ file.
    """
    # Read read IDs from the first FASTQ file
    read_ids1 = set()
    with open(fastq1, "r") as f:
        cont = 0
        for line in f:
            if line.startswith("@") and cont % 4 == 0:
                read_ids1.add(line.split()[0][1:])  # Remove '@' and take the first part of the header
                cont += 1
            # This handles spaces between reads
            elif cont % 4 != 0:
                cont += 1
    # Read read IDs from the second FASTQ file and find shared reads
    shared_reads = []
    with open(fastq2, "r") as f:
        cont = 0
        for line in f:
            if line.startswith("@") and cont % 4 == 0:
                read_id = line.split()[0][1:]  # Remove '@' and take the first part of the header
                if read_id in read_ids1:
                    shared_reads.append(read_id)
                cont += 1
            # This handles spaces between reads
            elif cont % 4 != 0:
                cont += 1

    # Save shared reads to output file
    with open(output_file, "w") as f:
        for read_id in shared_reads:
            f.write(f"{read_id}\n")


def index_fasta(fasta_file):
    """
    Index a FASTA file using samtools faidx.
    """
    cmd = ["samtools", "faidx", fasta_file]
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(" ".join(cmd), shell=True, check=True)


def plot_coverage(coverage_file, output_file, max_coverage=None):
    """
    Plot coverage from the coverage file and save the figure.
    """
    # Load the data
    df = pd.read_csv(coverage_file, sep='\t', header=None, 
                     names=['chr', 'pos', 'depth'])

    # Plot coverage
    plt.figure(figsize=(10, 5))

    if max_coverage is not None:
        # Define y max for plotting
        df['depth'] = df['depth'].clip(upper=max_coverage)
        # Add a horizontal line at the max coverage for visual reference
        plt.axhline(y=max_coverage, color='r', linestyle='--', 
                    label=f'Cap ({max_coverage})')
        plt.legend()
        plt.ylim(0, max_coverage)
    plt.plot(df['pos'], df['depth'], linewidth=0.5)
    plt.title('Read Depth Across Genome')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.tight_layout()
    # Save figure
    plt.savefig(output_file)
    plt.close()


def run_minimap2(fastq, bam_out, threads, ref_seq):
    """
    Run minimap2, optimised for long reads.
    """
    # Keep only mapped reads and output in BAM format
    cmd = [
        "minimap2",
        "-ax", "map-ont",
        "-t", str(threads),
        ref_seq,
        fastq,
        "|",
        "samtools", "view",
        "-F", "4",  # Filter out unmapped reads
        "-q", "30", # Filter by quality
        "-bS", "-",
        ">", bam_out
    ]
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(" ".join(cmd), shell=True, check=True)
