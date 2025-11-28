FROM mambaorg/micromamba:1.5.8 AS base

COPY docker_environment.yml /tmp/docker_environment.yml

ARG username=mopitaner

ENV DEBIAN_FRONTEND=noninteractive

RUN micromamba install -y -n base -f /tmp/docker_environment.yml && \
    micromamba clean --all --yes

# Install SeuratDisk inside the "base" env
RUN micromamba run -n base R -q -e "\
  if (!requireNamespace('remotes', quietly = TRUE)) \
    install.packages('remotes', repos = 'https://cloud.r-project.org'); \
  remotes::install_github('mojaveazure/seurat-disk')"

USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_PREFIX="/opt/conda" \
    CONDA_DEFAULT_ENV="mopitas-mapra" \
    CONDA_EXE="/opt/conda/bin/micromamba" \
    CONDA_PYTHON_EXE="/opt/conda/bin/python"

RUN apt-get update && apt-get install -y --no-install-recommends git curl \
    && curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && npm install -g @angular/cli \
    && apt-get clean


# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash $username && mkdir -p /workspaces/mopitas-mapra \
    && chown -R $username:$username /workspaces

USER $username
ENV USER=$username
WORKDIR /workspaces/mopitas-mapra
