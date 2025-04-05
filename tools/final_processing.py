import os, re
import pandas as pd
import json
from collections import defaultdict
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from tools.pipeline_utils import parse_memory_logs, create_stats_df, calculate_rpkm, calculate_tpm
from tools.plotting_utils import plot_runtimes, plot_memory_usage, plot_cpu_usage, plot_readstats

def process_pipeline_results(output_dir, samples, gene_lengths_path):
    """
    Processes pipeline results and generates summary statistics, TPM/RPKM values, and plots.

    Args:
        output_dir (str): Base output directory for the pipeline.
        samples (list): List of sample names.
        gene_lengths_path (str): Path to the gene lengths file.

    Returns:
        None
    """
    # Define output folders
    out_folder = os.path.join(output_dir, "combined_results")
    os.makedirs(out_folder, exist_ok=True)

    exon_out = os.path.join(out_folder, "exon")
    intron_out = os.path.join(out_folder, "intron")
    combined_out = os.path.join(out_folder, "combined_exon_intron")
    stats_out = os.path.join(out_folder, "stats")
    for folder in [exon_out, intron_out, combined_out, stats_out]:
        os.makedirs(folder, exist_ok=True)

    # Load gene lengths
    gene_lengths = pd.read_csv(gene_lengths_path, sep='\t', index_col=0)

    # Initialize data structures
    res_dict = defaultdict(lambda: defaultdict(dict))
    patients, reads = [], []
    exon_counts, intron_counts, combined_counts = [], [], []
    stats_list = []

    # Process each sample
    for sample in samples:
        folder = os.path.join(output_dir, sample)
        logs = [os.path.join(folder, "results", "logs", i) for i in os.listdir(os.path.join(folder, "results", "logs")) if "memory" in i]
        fastp_results = os.path.join(folder, "results", "fastp_stats", f"{sample}_fastp.json")
        stats_path = os.path.join(folder, "results", "stats", f"{sample}_expression_stats.txt")
        expression_folders = [os.path.join(folder, "results", "expression", i) for i in os.listdir(os.path.join(folder, "results", "expression"))]

        # Collect total reads after filtering
        with open(fastp_results, "r") as fh:
            fastp_data = json.load(fh)
        total_reads = fastp_data["summary"]["after_filtering"]["total_reads"]
        patients.append(sample)
        reads.append(total_reads)

        # Collect run parameters
        for log_file in logs:
            command = os.path.basename(log_file).split("_memory_log.txt")[0]
            runtime_seconds, mem_usage_mb, perc_cpu = parse_memory_logs(log_file)
            res_dict["runtime_seconds"][command][sample] = runtime_seconds
            res_dict["mem_usage_mb"][command][sample] = mem_usage_mb
            res_dict["perc_cpu"][command][sample] = perc_cpu

        # Collect expression data
        try:
            exon_path = next(filter(re.compile("exon").search, expression_folders))
        except StopIteration:
            exon_path = None
        try:
            intron_path = next(filter(re.compile("intron").search, expression_folders))
        except StopIteration:
            intron_path = None
        try:
            combined_path = next(filter(re.compile("combined").search, expression_folders))
        except StopIteration:
            combined_path = None

        if exon_path:
            exon_data = pd.read_csv(exon_path, sep='\t', index_col=0).rename({"Count": sample}, axis=1)
            exon_counts.append(exon_data)
        if intron_path:
            intron_data = pd.read_csv(intron_path, sep='\t', index_col=0).rename({"Count": sample}, axis=1)
            intron_counts.append(intron_data)
        if combined_path:
            combined_data = pd.read_csv(combined_path, sep='\t', index_col=0).rename({"Count": sample}, axis=1)
            combined_counts.append(combined_data)

        # Collect stats
        patient_stats = create_stats_df(stats_path, sample)
        stats_list.append(patient_stats)

    func_order = ["fastp", "fastqc_raw", "fastqc_trimmed",  "star_mapping", "featurecounts", "merge_bam", "quantify_expression"]
    func_order = [i for i in func_order if i in list(res_dict["runtime_seconds"].keys())]

    # Generate DataFrames
    read_depth_df = pd.DataFrame({"patient": patients, "reads": reads}).set_index("patient")
    runtime_df    = pd.DataFrame(res_dict["runtime_seconds"]).T.loc[func_order, :].sort_index(axis = 1)
    mem_usage_df  = pd.DataFrame(res_dict["mem_usage_mb"]).T.loc[func_order, :].sort_index(axis = 1)
    cpu_df        = pd.DataFrame(res_dict["perc_cpu"]).T.loc[func_order, :].sort_index(axis = 1)

    # Combine expression data
    exon_count_df = pd.concat(exon_counts, axis=1).fillna(0).sort_index().rename_axis(None)
    intron_count_df = pd.concat(intron_counts, axis=1).fillna(0).sort_index().rename_axis(None)
    combined_count_df = pd.concat(combined_counts, axis=1).fillna(0).sort_index().rename_axis(None)

    exon_tpm = calculate_tpm(exon_count_df, gene_lengths)
    intron_tpm = calculate_tpm(intron_count_df, gene_lengths)
    combined_tpm = calculate_tpm(combined_count_df, gene_lengths)

    exon_rpkm = calculate_rpkm(exon_count_df, gene_lengths)
    intron_rpkm = calculate_rpkm(intron_count_df, gene_lengths)
    combined_rpkm = calculate_rpkm(combined_count_df, gene_lengths)

    # Save expression data
    exon_count_df.to_csv(os.path.join(exon_out, "exon_counts.txt.gz"), sep='\t')
    intron_count_df.to_csv(os.path.join(intron_out, "intron_counts.txt.gz"), sep='\t')
    combined_count_df.to_csv(os.path.join(combined_out, "combined_counts.txt.gz"), sep='\t')

    exon_tpm.to_csv(os.path.join(exon_out, "exon_tpms.txt.gz"), sep = "\t")
    intron_tpm.to_csv(os.path.join(intron_out, "intron_tpms.txt.gz"), sep = "\t")
    combined_tpm.to_csv(os.path.join(combined_out, "combined_tpms.txt.gz"), sep = "\t")

    exon_rpkm.to_csv(os.path.join(exon_out, "exon_rpkms.txt.gz"), sep = "\t")
    intron_rpkm.to_csv(os.path.join(intron_out, "intron_rpkms.txt.gz"), sep = "\t")
    combined_rpkm.to_csv(os.path.join(combined_out, "combined_rpkms.txt.gz"), sep = "\t")

    # Save stats
    combined_stats_df = pd.concat(stats_list, axis = 0)
    combined_stats_df.to_csv(os.path.join(stats_out, "expression_stats.txt.gz"), sep = '\t')

    # Generate plots
    cmap = plt.cm.viridis
    norm = mcolors.Normalize(vmin=read_depth_df.min().values[0], vmax=read_depth_df.max().values[0])
    colors = cmap(norm(read_depth_df.reads.values))
    colormap = dict(zip(read_depth_df.index.values, colors))

    pipeline_stats_plot_path = os.path.join(out_folder, "plots", "pipeline_stats")
    os.makedirs(pipeline_stats_plot_path, exist_ok=True)
    readstats_plot_path = os.path.join(out_folder, "plots", "readstats")
    os.makedirs(readstats_plot_path, exist_ok=True)

    plot_runtimes(runtime_df, pipeline_stats_plot_path, cmap, norm, colormap)
    plot_memory_usage(mem_usage_df, pipeline_stats_plot_path, cmap, norm, colormap)
    plot_cpu_usage(cpu_df, pipeline_stats_plot_path, cmap, norm, colormap)
    plot_readstats(pd.concat(stats_list, axis=0), readstats_plot_path)