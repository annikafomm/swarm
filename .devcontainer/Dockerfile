FROM ubuntu:22.04

ARG username=mopitaner

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y git curl \
    && curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y nodejs \
    && apt-get clean \
    && npm install -g @angular/cli \
    && useradd -m $username

USER $username

WORKDIR /workspaces/mopitas-mapra

COPY --chown=$username ../frontend/package*.json ./frontend/

RUN cd frontend && npm install
