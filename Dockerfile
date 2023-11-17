FROM r-base

RUN apt-get update
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev

COPY ./*.R /usr/local/src/r_scripts/
WORKDIR /usr/local/src/r_scripts
RUN Rscript 000_packages_installation.R
