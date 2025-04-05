
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import os

def plot_runtimes(runtime_df, out_path, cmap, norm, colormap):
    fig, ax = plt.subplots()
    input_df = runtime_df[::-1]
    for idx, vals in enumerate(input_df.iterrows()):
        for patient, data in vals[1].items():
            ax.scatter(x = data / 60, y = idx, color = colormap[patient], alpha = 0.7, ec = 'k', s = 60)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  # Create a ScalarMappable for the colorbar
    sm.set_array([])  # Required for ScalarMappable
    cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical')  # Add the colorbar
    cbar.set_label("Sample read depth (after trimming)")  # Label for the colorbar
    
    ax.set_xlabel("Command runtime [m]")
    ax.set_xscale("linear")
    ax.set_yticks(range(len(input_df.index.values)))
    ax.set_yticklabels(input_df.index.values)
    
    sns.despine()
    plt.savefig(os.path.join(out_path, "command_runtimes.pdf"))

def plot_memory_usage(mem_usage_df, out_path, cmap, norm, colormap):
    fig, ax = plt.subplots()
    input_df = mem_usage_df[::-1]
    for idx, vals in enumerate(input_df.iterrows()):
        for patient, data in vals[1].items():
            ax.scatter(x = data, y = idx, color = colormap[patient], alpha = 0.7, ec = 'k', s = 60)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  # Create a ScalarMappable for the colorbar
    sm.set_array([])  # Required for ScalarMappable
    cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical')  # Add the colorbar
    cbar.set_label("Sample read depth [after trimming]")  # Label for the colorbar

    ax.set_xlabel("Command RAM usage [Mb]")
    ax.set_xscale("log")

    ax.set_yticks(range(len(input_df.index.values)))
    ax.set_yticklabels(input_df.index.values)

    sns.despine()
    plt.savefig(os.path.join(out_path, "command_memusage.pdf"))

def plot_cpu_usage(cpu_df, out_path, cmap, norm, colormap):
    fig, ax = plt.subplots()
    input_df = cpu_df[::-1]

    for idx, vals in enumerate(input_df.iterrows()):
        for patient, data in vals[1].items():
            ax.scatter(x = data, y = idx, color = colormap[patient], alpha = 0.7, ec = 'k', s = 60)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  # Create a ScalarMappable for the colorbar
    sm.set_array([])  # Required for ScalarMappable
    cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical')  # Add the colorbar
    cbar.set_label("Sample read depth [after trimming]")  # Label for the colorbar

    ax.set_xlabel("Command CPU usage [%]")
    ax.set_xscale("linear")

    ax.set_yticks(range(len(input_df.index.values)))
    ax.set_yticklabels(input_df.index.values)

    sns.despine()
    plt.savefig(os.path.join(out_path, "command_cpuusage.pdf"))

def plot_readstats(combined_stats_df, plot_path):
    fig, axes = plt.subplots(2, 2, figsize = (10, 10))
    ax1, ax2, ax3, ax4 = axes.ravel()
    
    kwargs ={
    "s":50,
    "alpha":0.6,
    "ec":'k',
    "color":"#32175f"}
    
    ax1.scatter(combined_stats_df.total_reads, combined_stats_df.perc_exon, **kwargs)
    ax2.scatter(combined_stats_df.total_reads, combined_stats_df.perc_intron, **kwargs)
    ax3.scatter(combined_stats_df.total_reads, combined_stats_df.perc_intergenic, **kwargs)
    ax4.scatter(combined_stats_df.total_reads, combined_stats_df.perc_unmapped, **kwargs)
    
    ax1.set_ylabel("Percent exon mapping reads")
    ax2.set_ylabel("Percent intron mapping reads")
    ax3.set_ylabel("Percent intergenic mapping reads")
    ax4.set_ylabel("Percent unmapped reads")
    
    ax1.set_xlabel("Total reads")
    ax2.set_xlabel("Total reads")
    ax3.set_xlabel("Total reads")
    ax4.set_xlabel("Total reads")
    
    [ax.set_xscale("log") for ax in axes.ravel()]
    fig.suptitle("Combined readstats for all samples")
    
    sns.despine()
    plt.savefig(os.path.join(plot_path, "readstats.pdf"))