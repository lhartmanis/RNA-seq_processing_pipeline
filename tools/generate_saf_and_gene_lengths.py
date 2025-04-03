import pandas as pd
import pyranges as pr
import argparse
import os

def generate_safs_and_gene_lengths(gtf_file, output_dir, generate_safs = True):
    print("Reading GTF file...")
    gtf = pr.read_gtf(gtf_file)
    exons = gtf[gtf.Feature == "exon"]
    genes = gtf[gtf.Feature == "gene"]

    print("Merging overlapping regions...")
    merged_exons = {gene: pr.PyRanges(gr).merge(strand = False) for gene, gr in exons.as_df().groupby("gene_id")}  # Merge overlapping exons
    exon_lengths = {gene: (df.End - df.Start).sum() for gene, df in merged_exons.items()}

    gene_spans = genes.df.groupby("gene_id").agg(Gene_Chrom=("Chromosome", "first"), 
                                                Gene_Start=("Start", "min"), 
                                                Gene_End=("End", "max"),
                                                Strand=("Strand", "first"))

    gene_spans["gene_length"] = gene_spans["Gene_End"] - gene_spans["Gene_Start"]

    ### Collect total exon and intron lenghts for each gene
    print("Calculating total length per gene...")
    exon_lengths_df = pd.DataFrame(list(exon_lengths.items()), columns=["gene_id", "exon_length"])
    gene_lengths_df = gene_spans.merge(exon_lengths_df, on="gene_id", how="left")
    gene_lengths_df["intron_length"] = gene_lengths_df["gene_length"] - gene_lengths_df["exon_length"]

    gene_lengths_df = gene_lengths_df[["gene_id", "gene_length", "exon_length", "intron_length"]]
    
    # Save the gene lengths file
    gene_lengths_file = f"{output_dir}/gene_lengths.txt"
    print(f"Saving gene lengths to {gene_lengths_file}...")
    gene_lengths_df.to_csv(gene_lengths_file, sep="\t", index=False)

    if generate_safs:
        if not os.path.exists("shared_data"):
            os.makedirs("shared_data")
        print("Generating SAF files...")
        # Exon SAF
        print("Making exon SAF...")
        exon_saf = exons.df[["gene_id", "Chromosome", "Start", "End", "Strand"]]
        exon_saf.columns = ["GeneID", "Chr", "Start", "End", "Strand"]
        exon_saf_file = "shared_data/exon_SAF.txt"
        print(f"Saving exon SAF to {exon_saf_file}...")
        exon_saf.to_csv(exon_saf_file, sep="\t", index=False, header=False)

        # Intron SAF
        print("Making intron SAF...")
        intron_list = []
        for gene, row in gene_spans.iterrows():
            chrom, gene_start, gene_end, strand = row["Gene_Chrom"], row["Gene_Start"], row["Gene_End"], row["Strand"]
            exon_regions = merged_exons.get(gene, pd.DataFrame(columns=["Start", "End"])).as_df()  # Get exons, or empty if none
            last_end = gene_start  # Start at gene start
            
            for _, exon in exon_regions.iterrows():
                intron_start = last_end
                intron_end = exon["Start"]
                if intron_end > intron_start:  # Only keep valid introns
                    intron_list.append([gene, chrom, intron_start, intron_end, strand])
                last_end = exon["End"]  # Move to end of this exon
            
            # Final intron between last exon and gene end
            if last_end < gene_end:
                intron_list.append([gene, chrom, last_end, gene_end, strand])
                
        intron_saf = pd.DataFrame(intron_list, columns=["GeneID", "Chr", "Start", "End", "Strand"])
        intron_saf_file = f"shared_data/intron_SAF.txt"
        print(f"Saving intron SAF to {intron_saf_file}...")
        intron_saf.to_csv(intron_saf_file, sep="\t", index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SAF files and gene lengths from a GTF file.")
    parser.add_argument("-g", "--gtf_file", required=True, help="Path to GTF file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for SAF files and gene lengths")
    parser.add_argument("-s", "--generate_safs", action="store_true", help="Generate SAF files")
    args = parser.parse_args()
    generate_safs_and_gene_lengths(args.gtf_file, args.output_dir, args.generate_safs)