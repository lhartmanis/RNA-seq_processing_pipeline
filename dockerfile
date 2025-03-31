# Use a specific Ubuntu base image
FROM ubuntu:20.04

# Set non-interactive mode to avoid prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Define tool versions
ENV FASTP_VERSION=0.23.2 \
    FASTQC_VERSION=0.11.9 \
    STAR_VERSION=2.7.10a \
    SUBREAD_VERSION=2.0.3 \
    SAMTOOLS_VERSION=1.17 \
    PYTHON_VERSION=3.9 \
    SNAKEMAKE_VERSION=7.32.4 \
    PYSAM_VERSION=0.21.0

# Update system and install core dependencies
RUN apt-get update && apt-get install -y \
    wget \
    python3=${PYTHON_VERSION} \
    python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools with specific versions
RUN wget -qO- https://github.com/OpenGene/fastp/releases/download/v${FASTP_VERSION}/fastp-${FASTP_VERSION}-linux64 | tee /usr/local/bin/fastp && chmod +x /usr/local/bin/fastp
RUN wget -qO /usr/local/bin/fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && chmod +x /usr/local/bin/fastqc
RUN wget -qO- https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz | tar xz -C /usr/local/bin && chmod +x /usr/local/bin/STAR
RUN wget -qO- https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz | tar xz -C /usr/local/bin && chmod +x /usr/local/bin/featureCounts
RUN wget -qO- https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xj -C /usr/local/bin && chmod +x /usr/local/bin/samtools

# Install Python dependencies
RUN pip3 install --no-cache-dir snakemake==${SNAKEMAKE_VERSION} pysam==${PYSAM_VERSION}

# Set the working directory
WORKDIR /pipeline

# Copy the pipeline scripts, Snakefile, and config files into the container
COPY . /pipeline

# Define an entrypoint so the container runs Snakemake by default
ENTRYPOINT ["snakemake"]

# Default command to use all available cores
CMD ["--cores", "ALL"]

#run:
# docker build -t rnaseq_pipeline
# docker run --rm -v $(pwd):/pipeline rnaseq_pipeline --cores 8
