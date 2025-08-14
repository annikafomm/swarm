FROM rocker/r-ver:4.5.0

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
    && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scanpy \
    anndata
RUN R -e "install.packages('BiocManager', repos='https://cran.rstudio.com/', dependencies=TRUE)" \
    && R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'readr', 'devtools'), repos='https://cran.rstudio.com/', dependencies=TRUE)" \
    && R -e "BiocManager::install(version='3.21')" \
    && R -e "BiocManager::install(c('AUCell', 'viper'), dependencies=TRUE, ask=FALSE, update=FALSE)" \
    && R -e "BiocManager::install('SPONGE', dependencies=TRUE, ask=FALSE, update=FALSE)"

    
USER $username


WORKDIR /workspaces/mopitas-mapra

