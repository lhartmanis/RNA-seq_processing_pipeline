import pandas as pd
import os
import re
import numpy as np

## Functions
def parse_wall_clock_time(time_str):
    """
    Parse a wall clock time string (h:mm:ss or m:ss) and return the total time in seconds.
    """
    # Match either h:mm:ss or m:ss format
    match = re.match(r"(?:(\d+):)?(\d+):(\d+)", time_str)
    if not match:
        raise ValueError(f"Invalid time format: {time_str}")
    
    # Extract hours, minutes, and seconds
    hours = int(match.group(1)) if match.group(1) else 0
    minutes = int(match.group(2))
    seconds = int(match.group(3))
    
    # Convert to total seconds
    total_seconds = hours * 3600 + minutes * 60 + seconds
    return total_seconds

def parse_memory_logs(log_file):
    """
    Parses memory log files to extract runtime (in seconds) and memory usage (in kB).
    
    Args:
    -----
        log_file (str): file path to memory log file.

    Returns:
    --------
        tuple:
            - runtimes (list of float): List of runtimes in seconds extracted from the logs.
            - reads (list of int): List of maximum memory usage (in kB) extracted from the logs.
    
    The function reads each log file line by line and looks for:
        - "Elapsed (wall clock) time": Extracts the runtime in seconds.
        - "Maximum resident set size": Extracts the peak memory usage in kilobytes (kB).
    """
    with open(log_file, "r") as f:
        for line in f:
            if "Elapsed (wall clock) time" in line:
               time_str = line.split("):")[-1].strip()
               runtime_seconds = parse_wall_clock_time(time_str)
            
            if "Maximum resident set size" in line:
               mem_usage_mb = int(np.round(int(line.split(":")[-1].strip()) / 1024, 0))
            
            if "Percent of CPU this job got" in line:
               perc_cpu = int(line.split(":")[-1].strip()[:-1])

    return runtime_seconds, mem_usage_mb, perc_cpu

def create_stats_df(stats_path, patient):
    """
    Creates a DataFrame summarizing sequencing statistics for a given patient.

    Args:
    -----
        stats_path (str): Path to the tab-delimited file containing sequencing statistics.
        patient (str): Identifier for the patient/sample, used to rename the "Count" column.

    Returns:
    --------
        pd.DataFrame: A DataFrame with sequencing statistics, including percentages of reads mapped to exons, introns, intergenic regions, and unmapped reads. The columns are reordered for clarity, and additional columns are included if present in the input file.

    Key Features:
    -------------
    - Reads the input file and renames the "Count" column to the patient identifier.
    - Calculates percentages for exon, intron, intergenic, and unmapped reads relative to total reads.
    
    """
    if not os.path.exists(stats_path):
        print(f"Warning: File not found for patient '{patient}': {stats_path}")
        return None
    stats_df = pd.read_csv(stats_path, sep = '\t', index_col = 0).rename({"Count": patient}, axis = 1).T.rename_axis(None, axis = 1)
    stats_df.columns = [i.lower() for i in stats_df.columns.values]
    
    stats_df.loc[:, "perc_exon"]   = (stats_df.exon / stats_df.total_reads)*100
    stats_df.loc[:, "perc_intron"] = (stats_df.intron / stats_df.total_reads) * 100
    stats_df.loc[:, "perc_intergenic"] = (stats_df.intergenic / stats_df.total_reads) * 100
    stats_df.loc[:, "perc_unmapped"] = (stats_df.unassigned_unmapped / stats_df.total_reads) * 100
    
    col_order = ["total_reads", "exon", "intron", "intergenic", "unassigned_unmapped", "perc_exon", "perc_intron", "perc_intergenic", "perc_unmapped"]
    rest = sorted(stats_df.loc[:, ~stats_df.columns.isin(col_order)].columns.values)
    col_order.extend(rest)
    
    stats_df = stats_df.loc[:, col_order]
    return stats_df

def calculate_rpkm(counts, gene_length_df):
    """
    Calculates RPKM (Reads Per Kilobase of transcript per Million mapped reads) for a given set of gene counts.
    
    Args:
    -----
        counts (pd.DataFrame): A DataFrame containing raw gene counts, where rows represent genes and columns represent samples.
        gene_length_df (pd.DataFrame): A DataFrame containing gene lengths, with an "exon_length" column specifying the length of each gene's exons in base pairs.

    Returns:
    --------
        pd.DataFrame: A DataFrame containing RPKM values for each gene (rows) and sample (columns).

    Key Steps:
    ----------
    1. Converts exon lengths from base pairs to kilobases (kb).
    2. Calculates the total counts per sample in millions (reads per million).
    3. Normalizes raw counts by exon length (RPK: Reads Per Kilobase).
    4. Divides RPK values by the total counts per million to compute RPKM.
    """
    exon_kb = gene_length_df.exon_length / 1000
    sample_count_per_million = counts.sum(0) / 1e6
    rpk = counts.div(exon_kb, axis = 0).dropna(how = "all")
    rpkm = rpk.div(sample_count_per_million, axis = 1)
    return rpkm

def calculate_tpm(counts, gene_length_df):
    """
    Calculates TPM (Transcripts Per Million) for a given set of gene counts.

    Args:
    -----
        counts (pd.DataFrame): A DataFrame containing raw gene counts, where rows represent genes and columns represent samples.
        gene_length_df (pd.DataFrame): A DataFrame containing gene lengths, with an "exon_length" column specifying the length of each gene's exons in base pairs.
    
    Returns:
    --------
        pd.DataFrame: A DataFrame containing TPM values for each gene (rows) and sample (columns).
    
    Key Steps:
    ----------
    1. Converts exon lengths from base pairs to kilobases (kb).
    2. Normalizes raw counts by exon length to calculate normalized counts.
    3. Scales normalized counts by the sum of normalized counts for each sample to compute TPM.
    4. Multiplies by 1e6 to express values as "per million."
    """
    exon_kb = gene_length_df.exon_length / 1000
    normed_counts = counts.div(exon_kb, axis = 0).dropna(how = "all")
    tpm = normed_counts.div(normed_counts.sum(axis = 0), axis = 1) * 1e6
    return tpm