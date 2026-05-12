FROM ubuntu:22.04

LABEL maintainer="abdelrahmanabdelmohsensaber@example.com"
LABEL description="Gene prediction benchmark environment"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    python3 \
    python3-pip \
    ncbi-blast+ \
    prodigal \
    augustus \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install pandas numpy biopython pysam scikit-learn xgboost

# Install MetaGeneAnnotator
RUN wget http://metagene.cb.k.u-tokyo.ac.jp/metagene_annotator.zip && \
    unzip metagene_annotator.zip && \
    cd metagene_annotator && make && cp metagene_annotator /usr/local/bin/

# Install BRAKER2
RUN git clone https://github.com/Gaius-Augustus/BRAKER.git && \
    cd BRAKER && make && cp braker.pl /usr/local/bin/

# Copy scripts
COPY scripts/ /opt/scripts/
RUN chmod +x /opt/scripts/*.py

# Set working directory
WORKDIR /data

# Default command
CMD ["python3", "/opt/scripts/calculate_gpcs.py"]