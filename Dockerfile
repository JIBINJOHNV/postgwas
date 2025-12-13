# =====================================================================
# Base Image: micromamba (Debian bullseye)
# =====================================================================
FROM mambaorg/micromamba:1.4.2

USER root

# =====================================================================
# System dependencies
# =====================================================================
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc g++ make wget git curl bzip2 \
        libcurl4-openssl-dev libsuitesparse-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        ca-certificates pkg-config \
        libssl-dev libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

# =====================================================================
# Install R + required packages (for SuSiE)
# =====================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('data.table','argparse','susieR','Matrix','dplyr'), repos='https://cloud.r-project.org/')"

# =====================================================================
# Install Docker CLI (STATIC BINARY — avoids libc mismatch)
# =====================================================================
RUN curl -L https://download.docker.com/linux/static/stable/x86_64/docker-24.0.6.tgz -o /tmp/docker.tgz && \
    tar -xzvf /tmp/docker.tgz -C /tmp && \
    mv /tmp/docker/docker /usr/local/bin/docker && \
    chmod +x /usr/local/bin/docker && \
    rm -rf /tmp/docker /tmp/docker.tgz

# =====================================================================
# Create Conda environment
# =====================================================================
COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n postgwas -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set default environment for interactive shells
RUN sed -i 's/micromamba activate base/micromamba activate postgwas/g' /etc/skel/.bashrc

# All following RUN commands executed inside postgwas environment
SHELL ["micromamba", "run", "-n", "postgwas", "/bin/bash", "-c"]

# =====================================================================
# Build HTSlib + BCFtools + Freeseek Plugins
# =====================================================================
WORKDIR /opt/tools

# Download sources
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xvjf bcftools-1.22.tar.bz2 && \
    tar -xvjf htslib-1.22.tar.bz2

# Build HTSlib
WORKDIR /opt/tools/htslib-1.22
RUN ./configure --enable-libcurl --prefix=/usr/local && \
    make -j && make install

WORKDIR /opt/tools/bcftools-1.22/plugins
RUN wget https://raw.githubusercontent.com/freeseek/score/master/score.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/score.h && \
    wget https://raw.githubusercontent.com/freeseek/score/master/munge.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/metal.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/blup.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/pgs.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/pgs.mk

# Patch pgs.c for SuiteSparse compatibility
WORKDIR /opt/tools/bcftools-1.22
RUN sed -i '2254s/^/\/\//' plugins/pgs.c && \
    sed -i '2255s/^/\/\//' plugins/pgs.c

# Build BCFtools + plugins
RUN ./configure \
        --prefix=/usr/local \
        --with-htslib=/opt/tools/htslib-1.22 \
        CPPFLAGS="-I/usr/include/suitesparse" \
        CFLAGS="-I/usr/include/suitesparse" && \
    make -j && make install

# =====================================================================
# Install PostGWAS
# =====================================================================
WORKDIR /opt/postgwas
COPY . /opt/postgwas

# Install MAGMA
RUN chmod +x magma/magma && mv magma/magma /usr/local/bin/magma

# Install PostGWAS itself
RUN pip install --upgrade pip && \
    pip install --no-deps --no-cache-dir -e .

# =====================================================================
# Switch back to host shell
# =====================================================================
SHELL ["/bin/bash", "-c"]

# Create runtime user
RUN useradd -m pguser
USER pguser

# =====================================================================
# Entrypoint
# =====================================================================
ENTRYPOINT ["micromamba", "run", "-n", "postgwas"]
CMD ["postgwas", "--help"]
