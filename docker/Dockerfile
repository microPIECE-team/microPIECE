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

RUN wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
    tar xzf bedtools.tar.gz && \
    cd bedtools2 && make && \
    mv bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedtools2 bedtools.tar.gz

RUN wget http://www.bioinf.uni-leipzig.de/Software/proteinortho/proteinortho_v5.16b.tar.gz && \
    tar xzf proteinortho_v5.16b.tar.gz && \
    rm -rf proteinortho_v5.16b.tar.gz
ENV PATH=/opt/proteinortho_v5.16b/:"$PATH"

RUN wget -O /usr/local/bin/xa2multi.pl https://raw.githubusercontent.com/lh3/bwa/master/xa2multi.pl && \
    chmod +x /usr/local/bin/xa2multi.pl

RUN wget https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.4.3/ViennaRNA-2.4.3.tar.gz && \
    tar -zxvf ViennaRNA-2.4.3.tar.gz && \
    cd ViennaRNA-2.4.3 && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf ViennaRNA-2.4.3.tar.gz ViennaRNA-2.4.3

VOLUME /data
WORKDIR /data