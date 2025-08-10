FROM rocker/r-ver:4.3.0

ARG username=mopitaner

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y git curl \
    && curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y nodejs \
    && apt-get clean \
    && npm install -g @angular/cli \
    && useradd -m $username \
    && apt-get install -y python3 python3-pip python3-venv \
    && apt-get install -y r-base r-base-dev \
    && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scanpy \
    anndata

RUN R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'readr', 'devtools', 'BiocManager'), repos='https://cran.rstudio.com/', dependencies=TRUE)" \
    && R -e "BiocManager::install(version='3.18')" \
    && R -e "BiocManager::install(c('AUCell', 'viper', 'SPONGE'), dependencies=TRUE, ask=FALSE, update=FALSE)"

USER $username

WORKDIR /workspaces/mopitas-mapra

