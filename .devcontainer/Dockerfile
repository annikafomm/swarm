FROM mambaorg/micromamba:1.5.8 AS base

COPY docker_environment.yml /tmp/docker_environment.yml
COPY tangram.yml /tmp/tangram.yml

# COPY multiome.yml /tmp/multiome.yml

ARG username=mopitaner

ENV DEBIAN_FRONTEND=noninteractive

RUN micromamba install -y -n base -f /tmp/docker_environment.yml && \
    micromamba clean --all --yes
RUN micromamba create -y -n tangram -f /tmp/tangram.yml && \
    micromamba clean --all --yes

RUN micromamba run -n base R -q -e "\
  if (!requireNamespace('remotes', quietly = TRUE)) \
    install.packages('remotes', repos = 'https://cloud.r-project.org'); \
  if (!requireNamespace('scCustomize', quietly = TRUE)) \
    install.packages('scCustomize', repos = 'https://cloud.r-project.org'); \
  remotes::install_github('mojaveazure/seurat-disk')"

RUN rm -rf /opt/conda/lib/R/library/Seurat \
           /opt/conda/lib/R/library/SeuratObject && \
    micromamba install -n base -y -c conda-forge --force-reinstall \
      r-seurat=5.3.0 \
      r-seuratobject=5.2.0

# RUN micromamba create -y -n multiome -f /tmp/multiome.yml && \
#     micromamba clean --all --yes


# # RUN micromamba run -n multiome R -q -e "\
# #   if (!requireNamespace('remotes', quietly = TRUE)) \
# #     install.packages('remotes', repos = 'https://cloud.r-project.org'); \
# #   remotes::install_github('mojaveazure/seurat-disk')"

#  sadly this does not work because of permission issues to install in multiome env
# RUN micromamba run -n multiome R -q -e "\
#   options(repos = c(CRAN='https://cloud.r-project.org')); \
#   lib <- file.path(Sys.getenv('CONDA_PREFIX'), 'lib/R/library'); \
#   dir.create(lib, recursive=TRUE, showWarnings=FALSE); \
#   if (isFALSE(requireNamespace('remotes', quietly=TRUE, lib.loc=lib))) install.packages('remotes', lib=lib); \
#   remotes::install_github('mojaveazure/seurat-disk', lib=lib); \
#   cat('SeuratDisk in multiome? ', requireNamespace('SeuratDisk', quietly=TRUE, lib.loc=lib), '\n')"


# Install SeuratDisk inside the "base" env
# RUN micromamba run -n base R -q -e "\
#   if (!requireNamespace('remotes', quietly = TRUE)) \
#     install.packages('remotes', repos = 'https://cloud.r-project.org'); \
#   remotes::install_github('mojaveazure/seurat-disk')"

RUN micromamba run -n base R -q -e "\
  if (!requireNamespace('remotes', quietly = TRUE)) \
    install.packages('remotes', repos = 'https://cloud.r-project.org'); \
  if (!requireNamespace('scCustomize', quietly = TRUE)) \
    install.packages('scCustomize', repos = 'https://cloud.r-project.org'); \
  remotes::install_github('mojaveazure/seurat-disk')"

#  aparently the wrong versions get installed, so here we install them correctly..
RUN rm -rf /opt/conda/lib/R/library/Seurat \
           /opt/conda/lib/R/library/SeuratObject && \
    micromamba install -n base -y -c conda-forge --force-reinstall \
      r-seurat=5.3.0 \
      r-seuratobject=5.2.0

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

# Install SSH to push with git (when cloning with ssh)
RUN apt-get install -y --no-install-recommends openssh-client

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash $username && mkdir -p /workspaces/mopitas-mapra \
    && chown -R $username:$username /workspaces

USER $username
ENV USER=$username
WORKDIR /workspaces/mopitas-mapra
