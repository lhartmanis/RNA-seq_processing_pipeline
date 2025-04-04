import pysam
import multiprocessing
import argparse
import logging
from collections import defaultdict

def setup_logging(log_file):
    """Set up logging to print to both console and a log file."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w")
        ]
    )

def count_reads_in_region(region, bam_path):
    """Count reads mapped to a specific region using pysam."""
    chrom, start, end = region
    bam = pysam.AlignmentFile(bam_path, "rb")
    exon_counts = defaultdict(int)
    intron_counts = defaultdict(int)
    combined_counts = defaultdict(int)
    
    stats = {
        "exon": 0,
        "intron": 0,
        "Unassigned_Unmapped": 0,
        "Unassigned_Read_Type": 0,
        "Unassigned_Singleton": 0,
        "Unassigned_MappingQuality": 0,
        "Unassigned_Chimera": 0,
        "Unassigned_FragmentLength": 0,
        "Unassigned_Duplicate": 0,
        "Unassigned_MultiMapping": 0,
        "Unassigned_Secondary": 0,
        "Unassigned_NonSplit": 0,
        "Unassigned_NoFeatures": 0,
        "Unassigned_Overlapping_Length": 0,
        "Unassigned_Ambiguity": 0,
        "Total_reads": 0
    }
    
    for read in bam.fetch(chrom, start, end):
        xf_tag = read.get_tag("XF") if read.has_tag("XF") else "Unassigned_NoFeatures"
        xt_tag = read.get_tag("XT") if read.has_tag("XT") else "Unassigned"
        
        # Convert XF tag to lowercase for consistency
        xf_tag = xf_tag.lower()
        
        if read.is_secondary or read.is_supplementary:
            stats["Unassigned_Secondary"] += 1
            continue
        
        if xf_tag == "exon":
            exon_counts[xt_tag] += 1
            combined_counts[xt_tag] += 1
            stats["exon"] += 1
        elif xf_tag == "intron":
            intron_counts[xt_tag] += 1
            combined_counts[xt_tag] += 1
            stats["intron"] += 1
        elif "unassigned" in xf_tag:
            for key in stats.keys():
                if key.lower() in xf_tag:
                    stats[key] += 1
        else:
            logging.warning(f"Unexpected XF tag: {xf_tag} for read {read.query_name}")
    
    bam.close()
    return exon_counts, intron_counts, combined_counts, stats

def count_total_reads(bam_path):
    """Count total number of reads in the BAM file."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    total_reads = sum(1 for _ in bam)
    bam.close()
    return total_reads

def merge_counts(results):
    """Merge counts from multiple threads."""
    exon_counts, intron_counts, combined_counts = defaultdict(int), defaultdict(int), defaultdict(int)
    stats = {
        "exon": 0,
        "intron": 0,
        "Unassigned_Unmapped": 0,
        "Unassigned_Read_Type": 0,
        "Unassigned_Singleton": 0,
        "Unassigned_MappingQuality": 0,
        "Unassigned_Chimera": 0,
        "Unassigned_FragmentLength": 0,
        "Unassigned_Duplicate": 0,
        "Unassigned_MultiMapping": 0,
        "Unassigned_Secondary": 0,
        "Unassigned_NonSplit": 0,
        "Unassigned_NoFeatures": 0,
        "Unassigned_Overlapping_Length": 0,
        "Unassigned_Ambiguity": 0,
        "Unassigned_Unmapped": 0, 
        "Total_reads": 0
    }
    
    for exon, intron, combined, stat in results:
        for gene, count in exon.items():
            exon_counts[gene] += count
        for gene, count in intron.items():
            intron_counts[gene] += count
        for gene, count in combined.items():
            combined_counts[gene] += count
        for key in stats:
            stats[key] += stat[key]
    
    return exon_counts, intron_counts, combined_counts, stats

def quantify_expression_parallel(bam_path, expression_folder, stats_folder, num_threads):
    """Quantify expression from a merged BAM file using multiprocessing."""
    logging.info("Starting parallel expression quantification...")
    bam = pysam.AlignmentFile(bam_path, "rb")
    regions = [(chrom, 0, bam.get_reference_length(chrom)) for chrom in bam.references]
    bam.close()
    
    num_threads = min(num_threads, len(regions))
    with multiprocessing.Pool(num_threads) as pool:
        results = pool.starmap(count_reads_in_region, [(region, bam_path) for region in regions])
    
    exon_counts, intron_counts, combined_counts, stats = merge_counts(results)
    
    # Compute total reads and assign unmapped category
    total_reads = count_total_reads(bam_path)
    assigned_unassigned_sum = sum(stats.values())
    stats["Unassigned_Unmapped"] = total_reads - assigned_unassigned_sum
    stats["Total_reads"] = total_reads
    stats["intergenic"] = stats.pop("Unassigned_NoFeatures")
    
    # Write output files
    with open(f"{expression_folder}_exon_counts.txt", "w") as f:
        f.write("GeneID\tCount\n")
        for gene, count in exon_counts.items():
            f.write(f"{gene}\t{count}\n")
    
    with open(f"{expression_folder}_intron_counts.txt", "w") as f:
        f.write("GeneID\tCount\n")
        for gene, count in intron_counts.items():
            f.write(f"{gene}\t{count}\n")
    
    with open(f"{expression_folder}_combined_counts.txt", "w") as f:
        f.write("GeneID\tCount\n")
        for gene, count in combined_counts.items():
            f.write(f"{gene}\t{count}\n")
    
    with open(f"{stats_folder}_expression_stats.txt", "w") as f:
        f.write("Category\tCount\n")
        for category, count in stats.items():
            f.write(f"{category}\t{count}\n")
    
    logging.info("Expression quantification completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quantify expression from a merged BAM file using multiprocessing.")
    parser.add_argument("-b", "--bam", required=True, help="Path to merged BAM file")
    parser.add_argument("-e", "--expression_folder", required=True, help="Path to output expression folder")
    parser.add_argument("-s", "--stats_folder", required=True, help="Path to output stats folder")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores)")
    parser.add_argument("--log", default="quantify_expression.log", help="Log file name (default: quantify_expression.log)")
    args = parser.parse_args()
    
    setup_logging(args.log)
    quantify_expression_parallel(args.bam, args.expression_folder, args.stats_folder, args.threads)

