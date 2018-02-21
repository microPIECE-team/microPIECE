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

VOLUME /data
WORKDIR /data