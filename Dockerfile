# ----------------------------------------------
# Base image: micromamba with Ubuntu base
# ----------------------------------------------
FROM mambaorg/micromamba:1.4.2

USER root

# Install build tools + git (required for vgraph)
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
