#!/usr/bin/env python3

import subprocess
import pysam


def align_index_minimap2(fastq, ref_seq, bam_out:str, threads:int):
    """
    Performs alignment with minimap2 and creates an index for the resulting BAM file.
    """
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


def extract_aligned_reads(bam_file, output_fastq, aln_len:int=1500):
    """
    Extract reads that align to the anchors from the BAM file 
    and save them in FASTQ format.
    """

    # Filter out reads with short alignment length
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
