FROM mambaorg/micromamba:1.4.2

USER root
RUN useradd -m pguser

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc g++ make git && \
    rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n postgwas -f /tmp/environment.yml && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "postgwas", "/bin/bash", "-c"]

WORKDIR /opt/postgwas
COPY . .

# Install MAGMA
RUN chmod +x magma/magma && \
    mv magma/magma /usr/local/bin/magma

# Install postgwas (editable mode)
RUN pip install --no-cache-dir -e .

USER pguser

CMD ["micromamba", "run", "-n", "postgwas", "postgwas", "--help"]
