#!/usr/bin/env python3

import pandas as pd
from pathlib import Path


def get_seq_from_aln(
    df:pd.DataFrame,
    fastq_file:Path,
    priority_anchor:str='left'
) -> pd.DataFrame:
    """
    Obtains the sequences from fastq file for the alignments in the input dataframe.
    The input dataframe must be output from the read_alignment_files function.
    """
    # Check that priority_anchor is either 'left' or 'right'
    if priority_anchor.lower() not in ['left', 'right']:
        raise ValueError("priority_anchor must be either 'left' or 'right'")
    else:
        priority_anchor = priority_anchor.lower()
        non_prio_anchor = 'right' if priority_anchor == 'left' else 'left'
    # Initialize a dictionary to store the sequences for each read_id
    seq_dict = {}
    # Go through the dataframe rows
    for index, row in df.iterrows():
        read_id = row['read_id']
        # Go through the fastq file and find the sequence for the read_id
        with fastq_file.open() as f:
            for line in f:
                if line.startswith('@') and line.strip()[1:] == read_id:
                    seq = next(f).strip()
                    seq_dict[read_id] = seq
    # Cut the sequences according to the alignment positions
    for index, row in df.iterrows():
        read_id = row['read_id']
        if priority_anchor == 'left':
            anchor_end = row['anchor_left_end']
            # Get sequence from anchor end to the end of the read
            seq = seq_dict[read_id][anchor_end:]
        else:
            anchor_start = row['anchor_right_start']
            # Get sequence from the start of the read to the anchor start
            seq = seq_dict[read_id][:anchor_start]
        # Update the sequence in the dictionary
        seq_dict[read_id] = seq
    # Add the sequences to the dataframe
    df['sequence'] = df['read_id'].map(seq_dict)
    return df


def read_alignment_files(
    input_folder:Path, 
    top_n:int, 
    read_id_mapping_file:Path,
    priority_anchor:str='left'
) -> pd.DataFrame:
    """
    Reads the alignment files in the input folder and returns a dictionary 
    with the read IDs as keys and the alignment scores as values.
    """
    # Check that priority_anchor is either 'left' or 'right'
    if priority_anchor.lower() not in ['left', 'right']:
        raise ValueError("priority_anchor must be either 'left' or 'right'")
    else:
        priority_anchor = priority_anchor.lower()
        non_prio_anchor = 'right' if priority_anchor == 'left' else 'left'
    # Define anchor name in the summary files
    LEFT_ANCHOR = 'Anchor 1 (left)'
    RIGHT_ANCHOR = 'Anchor 2 (right)'
    # Get the list of alignment files and summary files in the input folder
    alignment_files = list(input_folder.glob("*alignment.txt"))
    summary_files = list(input_folder.glob("*summary.csv"))
    # Initialize a dataframe for output data
    aln_cols = ['read_id', 'anchor_left_score', 'anchor_right_score', 
                'anchor_left_identity', 'anchor_right_identity',
                'anchor_left_start', 'anchor_left_end',
                'anchor_right_start', 'anchor_right_end']
    df_all_aln = pd.DataFrame(columns=aln_cols)
    # Go through summary files
    for summary_file in summary_files:
        print(f"Processing summary file: {summary_file}")
        # Read the summary file and get the read ID and alignment scores
        df_summary = pd.read_csv(summary_file)
        # Define left anchor rows
        left_anchor_rows = df_summary[df_summary['sequence']==LEFT_ANCHOR]
        right_anchor_rows = df_summary[df_summary['sequence']==RIGHT_ANCHOR]

        # Get the information from the table (uses hardcoded column ids)
        curr_read_id = summary_file.stem.rsplit("_", maxsplit=1)[0]
        left_score = left_anchor_rows.iloc[0, 1]
        right_score = right_anchor_rows.iloc[0, 1]
        left_identity = left_anchor_rows.iloc[0, -2]
        right_identity = right_anchor_rows.iloc[0, -2]
        left_start = left_anchor_rows.iloc[0, 2]
        left_end = left_anchor_rows.iloc[0, 3]
        right_start = right_anchor_rows.iloc[0, 2]
        right_end = right_anchor_rows.iloc[0, 3]

        # Load the information into the dataframe
        df_new_row = pd.DataFrame({
            'read_id':[curr_read_id], 
            'anchor_left_score':[float(str(left_score))], 
            'anchor_right_score':[float(str(right_score))], 
            'anchor_left_identity':[float(str(left_identity))], 
            'anchor_right_identity':[float(str(right_identity))],
            'anchor_left_start':[int(str(left_start))], 
            'anchor_left_end':[int(str(left_end))],
            'anchor_right_start':[int(str(right_start))], 
            'anchor_right_end':[int(str(right_end))]
        })
        df_all_aln = pd.concat([df_all_aln, df_new_row], ignore_index=True)

    # Select the top N alignments
    score_col_name = f'anchor_{priority_anchor}_score'
    df_all_aln[score_col_name] = df_all_aln[score_col_name].astype(float)
    df_top_n = df_all_aln.nlargest(top_n, score_col_name)

    # Map read_id to read_cont and add it to the top_n dataframe
    df_read_ids = pd.read_csv(read_id_mapping_file)
    # Add "read_" to the read_cont column in df_read_ids
    df_read_ids['read_cont'] = 'read_' + df_read_ids['read_cont'].astype(str)

    # Rename the read_id column in df_top_n to read_cont for merging
    df_top_n = df_top_n.rename(columns={'read_id':'read_cont'})

    # Merge the dataframes to get the read_id in the top_n dataframe
    df_top_n = df_top_n.merge(df_read_ids, on='read_cont', how='left')

    # Reorder columns so that read_id is the first column
    cols = df_top_n.columns.tolist()
    cols = ['read_id'] + [col for col in cols if col != 'read_id']
    df_top_n = df_top_n[cols]

    return df_top_n
