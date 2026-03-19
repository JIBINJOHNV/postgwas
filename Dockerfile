# =====================================================================
# Base Image
# =====================================================================
FROM mambaorg/micromamba:1.4.2

USER root
SHELL ["/bin/bash", "-c"]

# =====================================================================
# System dependencies (minimal)
# =====================================================================
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget curl git bzip2 unzip zip \
    build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libbz2-dev liblzma-dev zlib1g-dev \
    libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libfontconfig1-dev libeigen3-dev \
    binutils ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# =====================================================================
# 1. MAIN ENV (NO MKL)
# =====================================================================
COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n postgwas -f /tmp/environment.yml && \
    micromamba install -y -n postgwas -c conda-forge \
    gcc_linux-64 gxx_linux-64 sysroot_linux-64 \
    r-base r-xml2 r-matrix r-data.table r-survey && \
    micromamba run -n postgwas Rscript -e "\
    install.packages('remotes', repos='https://cloud.r-project.org'); \
    remotes::install_github('jianyanglab/gsmr2', dependencies=TRUE)" && \
    \
    # REMOVE BUILD DEPS
    micromamba remove -n postgwas -y \
    gcc_linux-64 gxx_linux-64 sysroot_linux-64 && \
    \
    # CLEAN ENV (CRITICAL)
    find /opt/conda/envs/postgwas -name "*.a" -delete && \
    find /opt/conda/envs/postgwas -type d -name "tests" -exec rm -rf {} + && \
    find /opt/conda/envs/postgwas -name "*.pyc" -delete && \
    \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/postgwas/bin:$PATH"

# =====================================================================
# 2. VEP ENV (CLEANED)
# =====================================================================
RUN micromamba create -y -n vep -c conda-forge -c bioconda \
    ensembl-vep=113 && \
    \
    find /opt/conda/envs/vep -name "*.a" -delete && \
    find /opt/conda/envs/vep -type d -name "tests" -exec rm -rf {} + && \
    \
    micromamba clean --all --yes

# =====================================================================
# 3. LDSC ENV (LEAN)
# =====================================================================
RUN micromamba create -y -n ldsc -c conda-forge -c bioconda \
    python=2.7 numpy=1.16 scipy=1.2.1 pandas=0.24.2 \
    bitarray && \
    \
    find /opt/conda/envs/ldsc -name "*.a" -delete && \
    micromamba clean --all --yes

WORKDIR /opt
RUN git clone --depth 1 https://github.com/bulik/ldsc.git

# =====================================================================
# 4. ENRICHER ENV (LEANED)
# =====================================================================
RUN micromamba create -y -n enricher -c conda-forge -c bioconda \
    omnipath gseapy zeep rpy2 pandas requests \
    networkx matplotlib scipy "numpy<2.0" \
    gprofiler-official && \
    \
    find /opt/conda/envs/enricher -name "*.a" -delete && \
    find /opt/conda/envs/enricher -type d -name "tests" -exec rm -rf {} + && \
    \
    micromamba clean --all --yes

# =====================================================================
# Build HTSlib + BCFtools (STRIPPED)
# =====================================================================
WORKDIR /opt/tools

RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xjf bcftools-1.22.tar.bz2 && \
    tar -xjf htslib-1.22.tar.bz2 && \
    cd htslib-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j && make install && \
    cd ../bcftools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j && make install && \
    \
    strip /usr/local/bin/* || true && \
    \
    cd / && rm -rf /opt/tools

# =====================================================================
# External tools (STRIPPED)
# =====================================================================
WORKDIR /tmp/tools

RUN wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz && \
    tar -xzf ldstore_v2.0_x86_64.tgz && \
    mv ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 /usr/local/bin/ldstore && \
    chmod +x /usr/local/bin/ldstore && \
    strip /usr/local/bin/ldstore || true && \
    \
    wget http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz && \
    tar -xzf finemap_v1.4.2_x86_64.tgz && \
    mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /usr/local/bin/finemap && \
    chmod +x /usr/local/bin/finemap && \
    strip /usr/local/bin/finemap || true && \
    \
    rm -rf /tmp/tools

# =====================================================================
# Install PostGWAS
# =====================================================================
WORKDIR /opt/postgwas
COPY . /opt/postgwas

RUN micromamba run -n postgwas pip install --no-cache-dir -e . && \
    micromamba clean --all --yes

# =====================================================================
# Final
# =====================================================================
ENV HOME=/tmp
RUN chown -R mambauser:mambauser /opt/postgwas

USER mambauser
ENTRYPOINT ["micromamba", "run", "-n", "postgwas"]