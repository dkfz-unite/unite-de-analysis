FROM ubuntu:latest

RUN apt-get update
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get install -y r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev

COPY . /usr/local/src/r_scripts
WORKDIR /usr/local/src/r_scripts
RUN Rscript 000_packages_installation.R
