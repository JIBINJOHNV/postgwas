FROM mambaorg/micromamba:1.4.2

# ----------------------------------------
# Create a user manually (fixes your error)
# ----------------------------------------
USER root
RUN useradd -m pguser

# System-level build tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc g++ make git && \
    rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml /tmp/environment.yml

# ----------------------------------------
# Create environment as root (safe)
# ----------------------------------------
RUN micromamba create -y -n postgwas -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Use environment for all future RUN commands
SHELL ["micromamba", "run", "-n", "postgwas", "/bin/bash", "-c"]

# ----------------------------------------
# Install your package
# ----------------------------------------
WORKDIR /opt/postgwas
COPY . /opt/postgwas

# Install MAGMA
RUN chmod +x magma/magma && \
    mv magma/magma /usr/local/bin/magma

# Install postgwas
RUN pip install -e .

# ----------------------------------------
# Switch to clean non-root user for runtime
# ----------------------------------------
USER pguser

CMD ["micromamba", "run", "-n", "postgwas", "postgwas", "--help"]
