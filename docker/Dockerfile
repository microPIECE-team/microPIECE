FROM ubuntu

LABEL maintainer="frank.foerster@ime.fraunhofer.de" \
      description="Container for the microPIECE package" \
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

RUN apt update && apt install --yes --no-install-recommends bamtools libgsl-dev && \
    wget http://smithlabresearch.org/downloads/piranha-1.2.1.tar.gz && \
    tar xzf piranha-1.2.1.tar.gz && \
    cd piranha-1.2.1 && \
    ./configure && \
    make all && \
    make test && \
    make install && \
    cd /opt && rm -rf piranha-1.2.1.tar.gz /var/lib/apt/lists/*

ENV PATH=/opt/piranha-1.2.1/bin/:"$PATH"

ENV PATH="$PATH":/opt/mirdeep2/bin
ENV PERL_MB_OPT="--install_base /root/perl5"
ENV PERL_MM_OPT="INSTALL_BASE=/root/perl5"
ENV PERL5LIB=/opt/mirdeep2/lib/perl5

RUN apt update && apt install --yes --no-install-recommends git && \
    git clone https://github.com/rajewsky-lab/mirdeep2.git && \
    cd mirdeep2 && \
    perl install.pl && \
    perl install.pl && \
    cd /opt && rm -rf /var/lib/apt/lists/*

RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-02-12.tar.gz && \
    tar xzf gmap-gsnap-2018-02-12.tar.gz && \
    cd gmap-2018-02-12 && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd /opt && rm -rf gmap-gsnap-2018-02-12.tar.gz gmap-2018-02-12

RUN apt update && apt install --yes --no-install-recommends git && \
    git clone https://github.com/lpantano/seqbuster.git && \
    cd /opt && rm -rf /var/lib/apt/lists/*

ENV MIRALIGNERDIR=/opt/seqbuster/modules/miraligner/
ENV MIRALIGNER="$MIRALIGNERDIR"/miraligner.jar

RUN wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz && \
    tar xzf miRanda-aug2010.tar.gz && \
    cd miRanda-3.3a && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd /opt && rm -rf miRanda-aug2010.tar.gz miRanda-3.3a

RUN apt update && apt install --yes --no-install-recommends git && \
    mkdir gffread && cd gffread && \
    git clone https://github.com/gpertea/gclib && \
    cd gclib && git checkout b790ac157971c5e6da77d0c76ee9ebf26ec4a5ef && cd .. && \
    git clone --branch v0.9.12 https://github.com/gpertea/gffread && \
    cd gffread && \
    make && \
    cd /opt && rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/gffread/gffread/:"$PATH"

RUN apt update && apt install --yes --no-install-recommends \
    liblog-log4perl-perl

VOLUME /data
WORKDIR /data
