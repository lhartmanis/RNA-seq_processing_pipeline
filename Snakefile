import time
from datetime import datetime
import shutil
from utils import check_executable, log_message, log_pipeline_runtime, log_software_versions, ensure_directories_exist, find_bam_file, cleanup_pipeline_results

if not os.path.exists(f"{config['output_dir']}"):
    os.makedirs(f"{config['output_dir']}")

if not os.path.exists(f"{config['temp_dir']}"):
    os.makedirs(f"{config['temp_dir']}")

shared_data_dir = f"{config['output_dir']}/shared_data"

# Determine the list of samples
if "samples_file" in config and config["samples_file"]:
    # Read sample names from the file
    with open(config["samples_file"], "r") as f:
        config["samples"] = [line.strip() for line in f if line.strip()]
elif "samples" in config and config["samples"]:
    # Use the samples list from the YAML file
    config["samples"] = config["samples"]
else:
    raise ValueError("No samples specified. Please provide either 'samples_file' or 'samples' in the configuration yaml file.")

# Ensure config["samples"] is always a list
if isinstance(config["samples"], str):
    config["samples"] = [config["samples"]]

# Ensure the required directories exist
ensure_directories_exist(
    base_dirs=["results/exon", "results/intron", "results/fastqc", "snakemake_checkpoints", "results/trimmed_fastq", "results/alignment", "results/expression", "results/stats"],
    output_dir=config["output_dir"],
    samples=config["samples"]
)

rule all:
    input:
        # Ensure SAF creation is completed
        f"{shared_data_dir}/saf_creation_done.txt" if config["generate_saf"] else f"{shared_data_dir}/saf_creation_skipped.txt",
        # Ensure renaming is completed
        expand("{output_dir}/{sample}/snakemake_checkpoints/renaming_done.txt", sample=config["samples"], output_dir=config["output_dir"]),
        # Include trimmed reads only if fastp is run
        expand("{output_dir}/{sample}/results/trimmed_fastq/{sample}_trimmed_R1.fastq.gz", sample=config["samples"], output_dir=config["output_dir"]) if config.get("run_fastp", False) else [],
        expand("{output_dir}/{sample}/results/trimmed_fastq/{sample}_trimmed_R2.fastq.gz", sample=config["samples"], output_dir=config["output_dir"]) if config.get("run_fastp", False) and config["read_type"] == "paired" else [],
        # FastQC outputs (before trimming)
        expand("{output_dir}/{sample}/results/fastqc/{sample}_R1_fastqc.html", sample=config["samples"], output_dir=config["output_dir"]) if config["fastqc_before_trimming"] else [],
        expand("{output_dir}/{sample}/results/fastqc/{sample}_R2_fastqc.html", sample=config["samples"], output_dir=config["output_dir"]) if config["fastqc_before_trimming"] and config["read_type"] == "paired" else [],
        # FastQC outputs (after trimming)
        expand("{output_dir}/{sample}/results/fastqc/{sample}_trimmed_R1_fastqc.html", sample=config["samples"], output_dir=config["output_dir"]) if config["fastqc_after_trimming"] else [],
        expand("{output_dir}/{sample}/results/fastqc/{sample}_trimmed_R2_fastqc.html", sample=config["samples"], output_dir=config["output_dir"]) if config["fastqc_after_trimming"] and config["read_type"] == "paired" else [],
        # Quantification outputs
        expand("{output_dir}/{sample}/results/expression/{sample}_exon_counts.txt", sample=config["samples"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/results/expression/{sample}_intron_counts.txt", sample=config["samples"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/results/expression/{sample}_combined_counts.txt", sample=config["samples"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/results/stats/{sample}_expression_stats.txt", sample=config["samples"], output_dir=config["output_dir"]),
        # Aligned BAM files
        expand("{output_dir}/{sample}/results/alignment/{sample}_intron_exon_aligned.bam", sample=config["samples"], output_dir=config["output_dir"]),
        # Cleanup and reverting outputs
        expand("{output_dir}/{sample}/snakemake_checkpoints/cleanup_done.txt", sample=config["samples"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/snakemake_checkpoints/reverting_done.txt", sample=config["samples"], output_dir=config["output_dir"])

rule generate_saf_and_gene_lengths:
    output:
        saf_creation_started = f"{shared_data_dir}/saf_creation_started.txt",
        saf_creation_status = f"{shared_data_dir}/saf_creation_done.txt" if config["generate_saf"] else f"{shared_data_dir}/saf_creation_skipped.txt",
        gene_lengths = f"{shared_data_dir}/gene_lengths.txt",
        exon_saf = f"{shared_data_dir}/exon_SAF.txt" if config["generate_saf"] else temp("dummy_exon_saf"),
        intron_saf = f"{shared_data_dir}/intron_SAF.txt" if config["generate_saf"] else temp("dummy_intron_saf")
    params:
        gtf_file = config["annotation_gtf"],
        output_dir = shared_data_dir,
        generate_saf = config["generate_saf"],
        script = "tools/generate_saf_and_gene_lengths.py",
        python_exec = config["executables"]["python"],
        saf_flag = "-s" if config["generate_saf"] else ""
    run:        
        if params.generate_saf:
            # Touch the "saf_creation_started.txt" file
            with open(output.saf_creation_started, "w") as f:
                f.write("SAF creation started.\n")

            # Run the SAF and gene lengths generation script
            shell(
                """
                {params.python_exec} {params.script} \
                    -g {params.gtf_file} \
                    -o {params.output_dir} \
                    {params.saf_flag}
                """
            )

            # Touch the "saf_creation_done.txt" file
            with open(output.saf_creation_status, "w") as f:
                f.write("SAF creation completed.\n")
        else:
            # Touch the "saf_creation_skipped.txt" file
            with open(output.saf_creation_status, "w") as f:
                f.write("SAF creation skipped.\n")

rule rename_files:
    input:
        # Ensure the input directory exists
        directory = config["input_dir"]
    output:
        renaming_done ="{output_dir}/{sample}/snakemake_checkpoints/renaming_done.txt",
        mapping_file = "{output_dir}/{sample}/renaming_map.json"
    run:
        from utils import rename_files
        # Perform renaming for the sample
        rename_files(
            input_dir = input.directory,
            config_file = "config.yaml",
            mapping_file = output.mapping_file
        )
        # Touch the renaming_done.txt file to indicate completion
        with open(output.renaming_done, "w") as f:
            f.write("Done renaming files!")

if config.get("run_fastp", False):
    rule fastp:
        input:
            renaming_done = "{output_dir}/{sample}/snakemake_checkpoints/renaming_done.txt"
        output:
            trimmed_read1 = "{output_dir}/{sample}/results/trimmed_fastq/{sample}_trimmed_R1.fastq.gz",
            trimmed_read2 = "{output_dir}/{sample}/results/trimmed_fastq/{sample}_trimmed_R2.fastq.gz" if config["read_type"] == "paired" else None,
            fastp_done = "{output_dir}/{sample}/snakemake_checkpoints/fastp_done.txt"
        params:
            fastp_exec = config["executables"]["fastp"],
            read1 = lambda wildcards: f"{config['input_dir']}/{wildcards.sample}_R1.fastq.gz",
            read2_param = lambda wildcards: f"-I {config['input_dir']}/{wildcards.sample}_R2.fastq.gz" if config["read_type"] == "paired" else "",
            trimmed_read2_param = lambda wildcards: f"-O {wildcards.output_dir}/{wildcards.sample}/results/trimmed_fastq/{wildcards.sample}_trimmed_R2.fastq.gz if config['read_type'] == 'paired' else ''",
            detect_adapter_for_pe = lambda wildcards: "--detect_adapter_for_pe" if config["read_type"] == "paired" else "",
            threads = config["fastp_threads"],
            extra = config.get("fastp_params", "")
        resources:
            mem_mb = lambda wildcards: int(os.path.getsize(f"{config['input_dir']}/{wildcards.sample}_R1.fastq.gz") / 1e6 * 2)  # 2x input file size in MB
        shell:
            """
            {params.fastp_exec} -i {params.read1} \
            {params.read2_param} \
            -o {output.trimmed_read1} \
            {params.trimmed_read2_param} \
            -w {params.threads} {params.detect_adapter_for_pe} {params.extra} && \
            touch {output.fastp_done}
            """
else:
    rule fastp_skipped:
        input:
            renaming_done = "{output_dir}/{sample}/snakemake_checkpoints/renaming_done.txt"
        output:
            fastp_skipped = "{output_dir}/{sample}/snakemake_checkpoints/fastp_skipped.txt"
        run:
            with open(output.fastp_skipped, "w") as f:
                f.write("Fastp was skipped for this sample.")

rule fastqc_raw:
    input:
        renaming_done = "{output_dir}/{sample}/snakemake_checkpoints/renaming_done.txt"
    output:
        raw_fastqc = [
            "{output_dir}/{sample}/results/fastqc/{sample}_R1_fastqc.html",
            "{output_dir}/{sample}/results/fastqc/{sample}_R2_fastqc.html"
        ] if config["fastqc_before_trimming"] and config["read_type"] == "paired" else (
            ["{output_dir}/{sample}/results/fastqc/{sample}_R1_fastqc.html"] if config["fastqc_before_trimming"] else []
        )
    params:
        raw = lambda wildcards: [
            f"{config['input_dir']}/{wildcards.sample}_R1.fastq.gz",
            f"{config['input_dir']}/{wildcards.sample}_R2.fastq.gz"
        ] if config["fastqc_before_trimming"] and config["read_type"] == "paired" else (
            [f"{config['input_dir']}/{wildcards.sample}_R1.fastq.gz"] if config["fastqc_before_trimming"] else []
        ),
        fastqc_exec = config["executables"]["fastqc"],
        out_folder = lambda wildcards: f"{wildcards.output_dir}/{wildcards.sample}/results/fastqc/"
    resources:
        mem_mb = config.get("fastqc_mem_mb", 1000)
    shell:
        """
        {params.fastqc_exec} {params.raw} -t 2 -o {params.out_folder}
        """

rule fastqc_trimmed:
    input:
        fastp_status = lambda wildcards: f"{wildcards.output_dir}/{wildcards.sample}/snakemake_checkpoints/fastp_done.txt" if config.get("run_fastp", False) else f"{wildcards.output_dir}/{wildcards.sample}/snakemake_checkpoints/fastp_skipped.txt"
    output:
        trimmed_fastqc = [
            "{output_dir}/{sample}/results/fastqc/{sample}_trimmed_R1_fastqc.html",
            "{output_dir}/{sample}/results/fastqc/{sample}_trimmed_R2_fastqc.html"
        ] if config["fastqc_after_trimming"] and config["read_type"] == "paired" else (
            ["{output_dir}/{sample}/results/fastqc/{sample}_trimmed_R1_fastqc.html"] if config["fastqc_after_trimming"] else []
        )
    params:
        trimmed = lambda wildcards: [
            f"{wildcards.output_dir}/{wildcards.sample}/results/trimmed_fastq/{wildcards.sample}_trimmed_R1.fastq.gz",
            f"{wildcards.output_dir}/{wildcards.sample}/results/trimmed_fastq/{wildcards.sample}_trimmed_R2.fastq.gz"
        ] if config["fastqc_after_trimming"] and config["read_type"] == "paired" else (
            [f"{wildcards.output_dir}/{wildcards.sample}/results/trimmed_fastq/{wildcards.sample}_trimmed_R1.fastq.gz"] if config["fastqc_after_trimming"] else []
        ),
        fastqc_exec = config["executables"]["fastqc"],
        out_folder = lambda wildcards: f"{wildcards.output_dir}/{wildcards.sample}/results/fastqc/"
    resources:
        mem_mb = 1000
    shell:
        """
        {params.fastqc_exec} {params.trimmed} -t 2 -o {params.out_folder}
        """

rule star_mapping:
    input:
        read1 = lambda wildcards: f"{wildcards.output_dir}/{wildcards.sample}/results/trimmed_fastq/{wildcards.sample}_trimmed_R1.fastq.gz" if config["run_fastp"] else f"{config['input_dir']}/{wildcards.sample}_R1.fastq.gz",
        read2 = lambda wildcards: f"{wildcards.output_dir}/{wildcards.sample}/results/trimmed_fastq/{wildcards.sample}_trimmed_R2.fastq.gz" if config["run_fastp"] and config["read_type"] == "paired" else f"{config['input_dir']}/{wildcards.sample}_R2.fastq.gz" if config["read_type"] == "paired" else None
    output:
        "{output_dir}/{sample}/results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    threads: config["star_threads"]
    params:
        tmp_dir = lambda wildcards: f"{config.get('temp_dir')}/{wildcards.sample}",
        genome_dir = config.get("star_genome_index", "/path/to/genomeDir"),
        gtf_file = config.get("annotation_gtf", "/path/to/annotation.gtf"),
        overhang = config.get("sjdb_overhang", 99),
        multimap_max = config.get("outFilterMultimapNmax", 50),
        sorting_threads = config.get("outBAMsortingThreadN", 10),
        adapter_seq = config.get("clip3AdapterSeq", "CTGTCTCTTATACACATCT"),
        quant_mode = config.get("quant_mode", "GeneCounts"),
        star_params = config.get("star_params", ""),
        star_exec = config["executables"]["star"]
    
    resources:
        mem_mb = config.get("star_mem_mb", 50000)  # Adjust memory as needed
    run:
        import shutil
        import os
        # Remove the temporary directory if it already exists, STAR crashes otherwise
        if os.path.exists(params.tmp_dir):
            shutil.rmtree(params.tmp_dir)
            print(f"Removed existing temporary directory: {params.tmp_dir}")

        # Run the STAR command
        shell(
            """
            echo "Running STAR for sample: {wildcards.sample}" && \
            echo " {params.star_exec} --runThreadN {threads} \
                 --genomeDir {params.genome_dir} \
                 --readFilesIn {input.read1} {input.read2} \
                 --readFilesCommand zcat \
                 --outFilterMultimapNmax {params.multimap_max} \
                 --outFileNamePrefix {wildcards.output_dir}/{wildcards.sample}/results/alignment/{wildcards.sample}_ \
                 --outSAMtype BAM SortedByCoordinate \
                 --readFilesType Fastx \
                 --clip3pAdapterSeq {params.adapter_seq} \
                 --outTmpDir {params.tmp_dir} \
                 --outSAMunmapped Within \
                 --outSAMmultNmax 1 \
                 --outBAMsortingThreadN {params.sorting_threads} \
                 --sjdbGTFfile {params.gtf_file} \
                 --sjdbOverhang {params.overhang} \
                 --quantMode {params.quant_mode} \
                 {params.star_params}" && \
            {params.star_exec} --runThreadN {threads} \
                 --genomeDir {params.genome_dir} \
                 --readFilesIn {input.read1} {input.read2} \
                 --readFilesCommand zcat \
                 --outFilterMultimapNmax {params.multimap_max} \
                 --outFileNamePrefix {wildcards.output_dir}/{wildcards.sample}/results/alignment/{wildcards.sample}_ \
                 --outSAMtype BAM SortedByCoordinate \
                 --readFilesType Fastx \
                 --clip3pAdapterSeq {params.adapter_seq} \
                 --outTmpDir {params.tmp_dir} \
                 --outSAMunmapped Within \
                 --outSAMmultNmax 1 \
                 --outBAMsortingThreadN {params.sorting_threads} \
                 --sjdbGTFfile {params.gtf_file} \
                 --sjdbOverhang {params.overhang} \
                 --quantMode {params.quant_mode} \
                 {params.star_params}
            """
        )
        
        # Move STAR log files to a STAR_output folder
        results_dir = f"{wildcards.output_dir}/{wildcards.sample}/results/alignment"
        star_output_dir = os.path.join(results_dir, "STAR_output")
        os.makedirs(star_output_dir, exist_ok=True)

        star_files = [
            f"{wildcards.sample}_Log.final.out",
            f"{wildcards.sample}_Log.out",
            f"{wildcards.sample}_Log.progress.out",
            f"{wildcards.sample}_ReadsPerGene.out.tab",
            f"{wildcards.sample}_SJ.out.tab",
            f"{wildcards.sample}__STARgenome"
        ]
        for file in star_files:
            src = os.path.join(results_dir, file)
            dst = os.path.join(star_output_dir, file)
            if os.path.exists(src):
                shutil.move(src, dst)

rule featurecounts:
    input:
        bam = lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/results/alignment/{wildcards.sample}_Aligned.sortedByCoord.out.bam",
        saf_creation_status = f"{shared_data_dir}/saf_creation_done.txt" if config["generate_saf"] else f"{shared_data_dir}/saf_creation_skipped.txt",
        exon_saf = f"{shared_data_dir}/exon_SAF.txt" if config["generate_saf"] else config["exon_saf"],
        intron_saf = f"{shared_data_dir}/intron_SAF.txt" if config["generate_saf"] else config["intron_saf"]
    output:
        exon = "{output_dir}/{sample}/results/exon/exon_counts.txt",
        intron = "{output_dir}/{sample}/results/intron/intron_counts.txt"
    threads: config.get("featurecounts_threads", 5)
    params:
        exon_dir = "{output_dir}/{sample}/results/exon",
        intron_dir = "{output_dir}/{sample}/results/intron",
        featurecounts_params = config.get("featurecounts_params", "-p --largestOverlap --primary"),
        featurecounts_exec = config["executables"]["featurecounts"]
    shell:
        """
        mkdir -p {params.exon_dir} {params.intron_dir} && \
        chmod u+w {params.exon_dir} {params.intron_dir} && \
        {params.featurecounts_exec} -F SAF -a {input.exon_saf} -o {output.exon} -T {threads} -R BAM {params.featurecounts_params} {input.bam} && \
        {params.featurecounts_exec} -F SAF -a {input.intron_saf} -o {output.intron} -T {threads} -R BAM {params.featurecounts_params} {input.bam}
        """

rule rename_featurecounts_bam:
    input:
        exon_counts = "{output_dir}/{sample}/results/exon/exon_counts.txt",
        intron_counts = "{output_dir}/{sample}/results/intron/intron_counts.txt",
    output:
        exon_bam = "{output_dir}/{sample}/results/exon/{sample}_exon.bam",
        intron_bam = "{output_dir}/{sample}/results/intron/{sample}_intron.bam"
    params:
        exon_dir = "{output_dir}/{sample}/results/exon",
        intron_dir = "{output_dir}/{sample}/results/intron"

    run:
        from utils import rename_featurecounts_bam
        exon_bam, intron_bam = rename_featurecounts_bam(
            wildcards.sample,
            params.exon_dir,
            params.intron_dir
        )
        assert exon_bam == output.exon_bam, f"Exon BAM file mismatch: {exon_bam} != {output.exon_bam}"
        assert intron_bam == output.intron_bam, f"Intron BAM file mismatch: {intron_bam} != {output.intron_bam}"


rule merge_bam:
    input:
        exon = "{output_dir}/{sample}/results/exon/{sample}_exon.bam",
        intron = "{output_dir}/{sample}/results/intron/{sample}_intron.bam"
    output:
        merged = "{output_dir}/{sample}/results/alignment/{sample}_intron_exon_aligned.bam",
        index = "{output_dir}/{sample}/results/alignment/{sample}_intron_exon_aligned.bam.bai"
    params:
        python_exec = config["executables"]["python"],
        priority = config.get("merge_priority", "exon"),  # Default priority is "exon"
        threads = config.get("merge_threads", 10),  # Default threads for sorting
        index_flag = lambda wildcards: "--index" if config.get("merge_index", True) else "",  # Dynamically set the --index flag
        log=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/results/{wildcards.sample}_merge_bam.log"
    shell:
        """
        {params.python_exec} tools/merge_exon_and_intron_bamfiles.py \
            -e {input.exon} \
            -i {input.intron} \
            -o {output.merged} \
            -p {params.priority} \
            -s \
            -t {params.threads} \
            {params.index_flag} \
            --log {params.log}
        """

rule quantify_expression:
    input:
        bam = lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/results/alignment/{wildcards.sample}_intron_exon_aligned.bam"
    output:
        exon_counts = "{output_dir}/{sample}/results/expression/{sample}_exon_counts.txt",
        intron_counts = "{output_dir}/{sample}/results/expression/{sample}_intron_counts.txt",
        combined_counts = "{output_dir}/{sample}/results/expression/{sample}_combined_counts.txt",
        expression_stats = "{output_dir}/{sample}/results/stats/{sample}_expression_stats.txt",
        log="{output_dir}/{sample}/results/{sample}_quant.log"
    params:
        python_exec = config["executables"]["python"],
        expression_folder = "{output_dir}/{sample}/results/expression/{sample}",
        stats_folder = "{output_dir}/{sample}/results/stats/{sample}",
        threads = config.get("quant_threads", 10)  
    shell:
        """
        {params.python_exec} tools/quantify_expression_generate_stats.py \
            -b {input.bam} \
            -e {params.expression_folder} \
            -s {params.stats_folder} \
            -t {params.threads} \
            --log {output.log}
        """
# This might need to change from aligned bam to count files.
rule cleanup:
    input:
        aligned_bam = "{output_dir}/{sample}/results/alignment/{sample}_intron_exon_aligned.bam"
    output:
        cleanup_done = "{output_dir}/{sample}/snakemake_checkpoints/cleanup_done.txt"
    params:
        output_dir = config["output_dir"],
        sample = lambda wildcards: wildcards.sample
    run:
        from utils import cleanup_pipeline_results
        cleanup_pipeline_results(params.output_dir, params.sample)
        # Touch the cleanup_done.txt file to indicate completion
        with open(output.cleanup_done, "w") as f:
            f.write("Cleanup completed for sample.")

rule revert_renaming:
    input:
        "{output_dir}/{sample}/snakemake_checkpoints/cleanup_done.txt"  
    output:
        reverting_done = "{output_dir}/{sample}/snakemake_checkpoints/reverting_done.txt"
    params:
        mapping_file = lambda wildcards: f"{wildcards.output_dir}/{wildcards.sample}/renaming_map.json"
    run:
        from utils import revert_renaming
        revert_renaming(params.mapping_file)
        with open(output.reverting_done, "w") as f:
            f.write("Reverting renaming completed.")
