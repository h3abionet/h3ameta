FROM quay.io/aptible/alpine
MAINTAINER Scott Hazelhurst "scott.hazelhurst@wits.ac.za"
RUN  apk update && apk add perl && apk add git && apk add && apk add bash && \
     apk add make && apk add g++  && \
     git clone https://github.com/DerrickWood/kraken2.git && \
     cd kraken2 && \
     ./install_kraken2.sh  /usr/local/bin
LABEL tool=kraken2 version=latest
