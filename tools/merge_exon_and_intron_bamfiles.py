import pysam
import argparse
import subprocess
import os
import time
import logging
from datetime import datetime

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

def check_samtools():
    """Check if samtools is installed before starting processing."""
    try:
        subprocess.run(["samtools", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        logging.info("Samtools is installed and available.")
    except FileNotFoundError:
        logging.error("Error: samtools is not installed or not in PATH.")
        exit(1)

def check_file_exists(file_path, file_type="BAM"):
    """Check if a file exists, else log an error and exit."""
    if not os.path.exists(file_path):
        logging.error(f"Error: {file_type} file '{file_path}' not found.")
        exit(1)

def get_tag(read, tag, default="Unassigned"):
    """Retrieve a BAM tag value, return default if missing."""
    return read.get_tag(tag) if read.has_tag(tag) else default

def propagate_tags(target_read, source_read):
    """Propagate XT, XS, and XF tags from source_read to target_read."""
    target_read.set_tag("XT", get_tag(source_read, "XT"), value_type="Z")
    target_read.set_tag("XS", get_tag(source_read, "XS", "Unassigned_NoFeatures"), value_type="Z")
    target_read.set_tag("XF", get_tag(source_read, "XF", "Unassigned_NoFeatures"), value_type="Z")
    target_read.flag = source_read.flag  # Ensure correct SAM flag propagation

def sort_and_index_bam(output_bam_path, num_threads, index_output):
    """Sort and optionally index the BAM file using samtools."""
    sorted_bam_path = output_bam_path.replace(".bam", ".sorted.bam")
    logging.info(f"Sorting BAM file: {output_bam_path} -> {sorted_bam_path}")
    
    try:
        subprocess.run(["samtools", "sort", "-@", str(num_threads), "-o", sorted_bam_path, output_bam_path], check=True)
        os.replace(sorted_bam_path, output_bam_path)  # Replace original with sorted BAM
    except subprocess.CalledProcessError as e:
        logging.error(f"Error sorting BAM file: {e}")
        exit(1)
    
    if index_output:
        logging.info(f"Indexing BAM file: {output_bam_path}")
        try:
            subprocess.run(["samtools", "index", output_bam_path], check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error indexing BAM file: {e}")
            exit(1)

def merge_bam_files(exon_bam_path, intron_bam_path, output_bam_path, priority, sort_output, num_threads, index_output):
    """Merge exon and intron BAM files, ensuring correct tag propagation and SAM flag consistency."""
    
    start_time = time.time()
    logging.info("Starting BAM merging process...")

    # Check for samtools before processing
    check_samtools()

    # Check if input files exist
    check_file_exists(exon_bam_path, "Exon BAM")
    check_file_exists(intron_bam_path, "Intron BAM")

    # Open BAM files
    try:
        exon_bam = pysam.AlignmentFile(exon_bam_path, "rb")
        intron_bam = pysam.AlignmentFile(intron_bam_path, "rb")
        out_bam = pysam.AlignmentFile(output_bam_path, "wb", header=exon_bam.header)
    except Exception as e:
        logging.error(f"Failed to open BAM files: {e}")
        exit(1)

    intron_reads = {}

    logging.info("Processing intron BAM file...")
    for read in intron_bam:
        intron_reads[read.query_name] = read

    logging.info("Processing exon BAM file and writing to output...")
    for read in exon_bam:
        exon_xs = get_tag(read, "XS", "Unassigned_NoFeatures")
        intron_xs = get_tag(intron_reads.get(read.query_name, read), "XS", "Unassigned_NoFeatures")

        if exon_xs == "Assigned" and intron_xs == "Assigned":
            chosen_read = read if priority == "exon" else intron_reads[read.query_name]
            propagate_tags(read, chosen_read)
            read.set_tag("XF", priority, value_type="Z")
        elif exon_xs == "Assigned":
            propagate_tags(read, read)
            read.set_tag("XF", "exon", value_type="Z")
        elif intron_xs == "Assigned":
            propagate_tags(read, intron_reads[read.query_name])
            read.set_tag("XF", "intron", value_type="Z")
        else:
            read = read if priority == "exon" else intron_reads[read.query_name]
            read.set_tag("XT", "Unassigned", value_type="Z")
            read.set_tag("XS", read.get_tag("XS"))
            read.set_tag("XF", read.get_tag("XS"), value_type="Z")
        
        out_bam.write(read)
        if read.query_name in intron_reads:
            del intron_reads[read.query_name]  # Remove processed intron reads

    if intron_reads:
        logging.info(f"Writing {len(intron_reads)} remaining intron-only reads...")
        for read in intron_reads.values():
            propagate_tags(read, read)
            read.set_tag("XF", "intron", value_type="Z")  # Ensure intron-only reads get "intron"
            out_bam.write(read)

    exon_bam.close()
    intron_bam.close()
    out_bam.close()

    logging.info(f"Merged BAM file created: {output_bam_path}, Overlapping reads set to '{priority}'")
    
    if sort_output:
        sort_and_index_bam(output_bam_path, num_threads, index_output)
    
    total_time = time.time() - start_time
    logging.info(f"Process completed in {total_time:.2f} seconds.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge exon and intron BAM files, ensuring correct tag propagation.")
    parser.add_argument("-e", "--exon", required=True, help="Path to exon-mapped BAM file")
    parser.add_argument("-i", "--intron", required=True, help="Path to intron-mapped BAM file")
    parser.add_argument("-o", "--output", required=True, help="Path to output merged BAM file")
    parser.add_argument("-p", "--priority", choices=["exon", "intron"], default="exon",
                        help="Which feature to prioritize for overlapping reads (default: exon)")
    parser.add_argument("-s", "--sort", action="store_true",
                        help="Sort the output BAM file by coordinate (requires samtools)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Number of threads to use for sorting (default: 4)")
    parser.add_argument("--index", action="store_true",
                        help="Index the output BAM file after sorting (requires samtools)")
    parser.add_argument("--log", default="merge_bam.log",
                        help="Log file name (default: merge_bam.log)")
    args = parser.parse_args()
    setup_logging(args.log)
    merge_bam_files(args.exon, args.intron, args.output, args.priority, args.sort, args.threads, args.index)

