# ===============================
# Base Image
# ===============================
FROM mambaorg/micromamba:1.5.1

# Switch to root to install system packages
USER root

# Install compiler tools required for vgraph
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc g++ make && \
    rm -rf /var/lib/apt/lists/*

# Switch back to micromamba's user (correct user = mamba)
USER mamba

# ===============================
# Conda/Mamba Environment
# ===============================

ARG MAMBA_ENV_NAME=postgwas

# Copy environment.yml with proper permissions
COPY environment.yml /tmp/environment.yml
RUN chown mamba:mamba /tmp/environment.yml

# Create conda environment
RUN micromamba create -y -n ${MAMBA_ENV_NAME} -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate environment for all following RUN commands
SHELL ["micromamba", "run", "-n", "postgwas", "/bin/bash", "-c"]

# ===============================
# Install PostGWAS + MAGMA
# ===============================

WORKDIR /opt/postgwas
COPY . /opt/postgwas

# Install MAGMA (Linux static binary included inside repo)
RUN chmod +x magma/magma && \
    mv magma/magma /usr/local/bin/magma && \
    mv magma/aux /usr/local/bin/aux

# Install postgwas package
RUN pip install -e .

# Default command
CMD ["micromamba", "run", "-n", "postgwas", "postgwas", "--help"]
