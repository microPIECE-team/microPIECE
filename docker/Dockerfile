FROM ubuntu

LABEL maintainer="frank.foerster@ime.fraunhofer.de" \
      description="Container for the MiCPi package" \
      version="0.1"

RUN apt update && \
    apt install \
       --yes \
       --no-install-recommends \
       build-essential \
       wget \
       bwa \
       bowtie \
       bowtie2 \
       ncbi-blast+ \
       python-cutadapt \
       emboss \
       openjdk-8-jre-headless \
       && \
       rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN apt update && apt install --yes --no-install-recommends \
       zlib1g-dev \
       libncurses5-dev \
       libbz2-dev \
       liblzma-dev \
       && \
    wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 && \
    tar xjf samtools-1.4.1.tar.bz2 && \
    cd samtools-1.4.1 && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd .. && \
    rm -rf samtools-1.4.1 samtools-1.4.1.tar.bz2 /var/lib/apt/lists/*



VOLUME /data
WORKDIR /data