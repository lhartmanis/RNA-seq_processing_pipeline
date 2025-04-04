# RNA-seq Processing Pipeline

This repository contains a Snakemake-based pipeline for processing RNA-seq data. The pipeline automates tasks such as quality control, read trimming, alignment and feature quantification, making it easy to process RNA-seq datasets reproducibly and efficiently.

---

## **Features**
- **Quality Control:** Runs FastQC before and after trimming.
- **Read Trimming:** Uses `fastp` for adapter trimming and quality filtering.
- **Alignment:** Aligns reads to a reference genome using STAR.
- **Feature Quantification:** Quantifies exon and intron counts using featureCounts.
- **BAM Merging:** Merges exon and intron BAM files into a single file.
- **Expression Quantification:** Generates combined expression statistics.
- **Cleanup:** Removes intermediate files to save disk space.
- **Reverting Renaming:** Reverts renamed files to their original names.

---

## **Requirements**
### **Software**
The following software must be installed and available in your `PATH`:
- [Snakemake](https://snakemake.readthedocs.io/)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [STAR](https://github.com/alexdobin/STAR)
- [featureCounts](http://subread.sourceforge.net/)
- Python 3.x:
    - pyranges
    - pandas
    - numpy

### **Input Files**
- Paired-end or single-end FASTQ files are placed in the `input dir` folder and named according to: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz` (for paired-end reads). 
- Reference genome index for STAR.

## **Installation**
1. Clone this repository:
   ```bash
   git clone https://github.com/lhartmanis/RNA-seq_processing_pipeline.git
   cd RNA_seq_processing_pipeline
2. Install the required software and dependencies.
3. Edit the config.yaml file to specify your input files, output directory, and pipeline parameters.

## **Usage**
### **1. Configure the Pipeline**
Edit the `config.yaml` file to include:
- Input directory (`input_dir`) containing FASTQ files.
- Output directory (`output_dir`) for pipeline results.
- List of sample names (`samples`).
- Paths to reference genome and annotation files.

### **2. Run the Pipeline**
Run the pipeline using Snakemake:
```bash
snakemake --cores <number_of_cores>
```

### **3. Capture Logs**
To save logs to a file:
```bash
snakemake --cores <number_of_cores> > pipeline.log 2>&1
```

## **Pipeline Workflow**
1. **Renaming Files:** Renames input FASTQ files to a standardized format.
2. **Quality Control:** Runs FastQC on raw reads.
3. **Read Trimming:** Trims adapters and filters low-quality reads using `fastp`.
4. **Alignment:** Aligns trimmed reads to the reference genome using STAR.
5. **Feature Quantification:** Quantifies exon and intron counts using featureCounts.
6. **BAM Merging:** Merges exon and intron BAM files into a single BAM file.
7. **Expression Quantification:** Generates combined expression statistics to account for exonic and intronic mapping reads.
8. **Cleanup:** Removes intermediate files to save disk space.
9. **Reverting Renaming:** Reverts renamed files to their original names.

## **Configuration**
The pipeline is configured using the `config.yaml` file. Below is an example configuration:

```yaml
input_dir: "/path/to/input"
output_dir: "/path/to/output"
samples:
  - sample1
  - sample2
run_fastp: true
fastqc_before_trimming: true
fastqc_after_trimming: true
read_type: "paired"
star_genome_index: "/path/to/star/genomeDir"
annotation_gtf: "/path/to/annotation.gtf"
exon_saf: "/path/to/exon_SAF"
intron_saf: "/path/to/intron_SAF"
executables:
  fastp: "fastp"
  fastqc: "fastqc"
  star: "STAR"
  featurecounts: "featureCounts"
  python: "python3"
```

## **Output**
The pipeline generates the following outputs for each sample:

### **Quality Control Reports**
- `{output_dir}/{sample}/results/fastqc/{sample}_R1_fastqc.html`
- `{output_dir}/{sample}/results/fastqc/{sample}_R2_fastqc.html` (if paired-end)
- `{output_dir}/{sample}/results/fastqc/{sample}_trimmed_R1_fastqc.html`
- `{output_dir}/{sample}/results/fsatqc/{sample}_trimmed_R2_fastqc.html` (if paired-end)

### **Trimmed Reads**
- `{output_dir}/{sample}/results/trimmed_fastq/{sample}_trimmed_R1.fastq.gz`
- `{output_dir}/{sample}/results/trimmed_fastq/{sample}_trimmed_R2.fastq.gz` (if paired-end)

### **Aligned BAM Files**
- `{output_dir}/{sample}/results/alignment/{sample}_Aligned.sortedByCoord.out.bam`

### **Quantification Results**
- `{output_dir}/{sample}/results/expression/{sample}_exon_counts.txt`
- `{output_dir}/{sample}/results/expression/{sample}_intron_counts.txt`
- `{output_dir}/{sample}/results/expression/{sample}_combined_counts.txt`
- `{output_dir}/{sample}/results/stats/{sample}_expression_stats.txt`

### **Cleanup and Logs**
- `{output_dir}/{sample}/cleanup_done.txt`
- `{output_dir}/{sample}/reverting_done.txt`
- Log files for each step are stored in the `results` directory.

## **Contributing**
Contributions are welcome! If you would like to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bug fix:
   ```bash
   git checkout -b feature-name
   ````

3. Make your changes and commit them:
    ```bash
    git commit -m "Description of changes"

4. Push your changes to your fork;
    ```bash
    git push origin feature-name
    ````
5. Submit a pull request to the main repository

## **License**
This project is licensed under the MIT License. See the `LICENSE` file for details.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## **Acknowledgments**
This pipeline was made possible thanks to the following tools and resources:
- [Snakemake](https://snakemake.readthedocs.io/) for workflow management.
- [fastp](https://github.com/OpenGene/fastp) for read trimming and quality filtering.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control.
- [STAR](https://github.com/alexdobin/STAR) for RNA-seq alignment.
- [featureCounts](http://subread.sourceforge.net/) for feature quantification.
- The open-source community for providing these tools and resources.

Special thanks to the contributors and maintainers of these tools for their invaluable work.
