# GapWalker - Long-Read Gap‑Closure Pipeline

GapWalker extends a DNA sequence from one side of a scaffolding gap to the other using
ONT long reads. Starting from two known anchor sequences, it iteratively assembles the
intervening region and stops when the far‑side anchor is reached.

## Overview

The pipeline works in five main steps (plus an optional visualisation utility):

1. **extract‑reads** - align all reads to both anchors, extract reads that map to both
2. **align‑reads** - perform a second, more precise dual alignment of each shared read
3. **get‑top‑align** - rank reads by alignment quality, extract the extension
   sequence, and build a consensus
4. **check‑and‑iter** - test if the consensus contains the other anchor; if not, use
   the last 2 kb of the new consensus as a fresh anchor for the next iteration
5. **concatenate‑evaluate** - join all consensus segments, trim at the detected anchor,
   and map original reads back to visualise coverage

A unified CLI (`python cli_main.py`) orchestrates the entire process with a single
`run‑all` command, which loops steps 1‑4 until the gap is closed (or a maximum number
of iterations is reached). Individual steps can also be executed separately for
troubleshooting or custom workflows.

## Installation

### External tools

GapWalker requires the following command‑line programs on your `PATH`:

| Tool      | Why                            | Installation                          |
|-----------|--------------------------------|---------------------------------------|
| minimap2  | long‑read alignment to anchors | `conda install -c bioconda minimap2`  |
| samtools  | BAM/FASTQ manipulation         | `conda install -c bioconda samtools`  |
| mafft     | optional, better MSA than the built‑in progressive pairwise aligner | `conda install -c bioconda mafft` |

### Python environment

**Option A - Conda (recommended)**

```bash
conda env create -f environment.yml
conda activate gapwalker
```

**Option B - pip + venv**

```bash
python3 -m venv gapwalker_env
source gapwalker_env/bin/activate
pip install -r requirements.txt
```

Note: MAFFT is optional but strongly recommended when you have more than
a few reads; the built‑in progressive aligner can handle small sets effectively but
may be slower.

## Quick start

Run the whole pipeline with default settings (max 10 iterations, left anchor as
extension source):

```bash
python cli_main.py run-all \
    -1 anchor_left.fasta \
    -2 anchor_right.fasta \
    -i reads.fastq \
    -o ./my_gap_results \
    --max-iter 10 \
    --threads 8
```