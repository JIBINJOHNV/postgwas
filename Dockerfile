# ----------------------------------------------
# Base image: micromamba with Ubuntu base
# ----------------------------------------------
FROM mambaorg/micromamba:1.4.2
USER root
# Install initial build tools + git (required for vgraph)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc g++ make git && \
    rm -rf /var/lib/apt/lists/*
# Copy environment file
COPY environment.yml /tmp/environment.yml
# ----------------------------------------------
# Create the conda environment
# ----------------------------------------------
ARG MAMBA_ENV_NAME=postgwas
RUN micromamba create -y -n ${MAMBA_ENV_NAME} -f /tmp/environment.yml && \
    micromamba clean --all --yes
# All further RUN commands executed inside env
SHELL ["micromamba", "run", "-n", "postgwas", "/bin/bash", "-c"]

# ----------------------------------------------
# Install BCFtools and Freeseek Plugins
# ----------------------------------------------
WORKDIR /opt/tools
# Install prerequisites: SuiteSparse, wget, bzip2, libcurl, Zlib, Bzip2, and Lzma development files
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget bzip2 libcurl4-openssl-dev libsuitesparse-dev zlib1g-dev libbz2-dev liblzma-dev && \
    rm -rf /var/lib/apt/lists/*
# Download bcftools + htslib
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xvjf bcftools-1.22.tar.bz2 && \
    tar -xvjf htslib-1.22.tar.bz2
# Download Freeseek plugins into bcftools plugins directory
WORKDIR /opt/tools/bcftools-1.22/plugins
RUN wget http://raw.githubusercontent.com/freeseek/score/master/score.c && \
    wget http://raw.githubusercontent.com/freeseek/score/master/score.h && \
    wget http://raw.githubusercontent.com/freeseek/score/master/munge.c && \
    wget http://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    wget http://raw.githubusercontent.com/freeseek/score/master/metal.c && \
    wget http://raw.githubusercontent.com/freeseek/score/master/blup.c && \
    wget http://raw.githubusercontent.com/freeseek/score/master/pgs.c && \
    wget http://raw.githubusercontent.com/freeseek/score/master/pgs.mk
# ---
# Build HTSlib + bcftools
# ---
WORKDIR /opt/tools/htslib-1.22
RUN ./configure --enable-libcurl --prefix=/usr/local && \
    make && \
    make install
WORKDIR /opt/tools/bcftools-1.22
# PATCH: Comment out CHOLMOD v4+ incompatible code in pgs.c (lines 2254 and 2255)
RUN sed -i '2254s/^/\/\//' plugins/pgs.c && \
    sed -i '2255s/^/\/\//' plugins/pgs.c
RUN export CFLAGS="-I/usr/include/suitesparse" && \
    ./configure --prefix=/usr/local --with-htslib=/opt/tools/htslib-1.22 && \
    make && \
    make install

# ----------------------------------------------
# Install the postgwas package
# ----------------------------------------------
WORKDIR /opt/postgwas
COPY . /opt/postgwas
# Install MAGMA (Linux binary) — your folder tree shows magma/magma exists
RUN chmod +x magma/magma && \
    mv magma/magma /usr/local/bin/magma
# Install postgwas into conda env
RUN pip install --no-cache-dir .
# ----------------------------------------------
# Create runtime user (non-root)
# ----------------------------------------------
USER root
RUN useradd -m pguser
USER pguser
# ----------------------------------------------
# ENTRYPOINT — always run inside conda env
# This allows:
#   docker run image postgwas
#   docker run image postgwas finemap --help
#   docker run image bash
# ----------------------------------------------
ENTRYPOINT ["micromamba", "run", "-n", "postgwas"]
CMD ["postgwas", "--help"]