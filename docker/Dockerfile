ARG baseimage_version=latest
FROM micropiece/micropiece-dockerbase:${baseimage_version}

ARG micropiece_version=master

LABEL maintainer="frank.foerster@ime.fraunhofer.de" \
      description="Container for the microPIECE package" \
      version="1.0"

RUN apt update && apt install --yes --no-install-recommends git && \
    git clone --branch ${micropiece_version} https://github.com/microPIECE-team/microPIECE.git /opt/microPIECE && \
    rm -rf /opt/microPIECE/.git && \
    rm -rf /var/lib/apt/lists/*
ENV PATH=/opt/microPIECE/:"$PATH"
