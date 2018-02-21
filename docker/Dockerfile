FROM ubuntu

LABEL maintainer="frank.foerster@ime.fraunhofer.de" \
      description="Container for the MiCPi package" \
      version="0.1"

RUN apt update && \
    apt install \
       --yes \
       --no-install-recommends \
       build-essential \
       && \
       rm -rf /var/lib/apt/lists/*