FROM ubuntu:18.04

MAINTAINER Scott Hazelhurst "scott.hazelhurst@wits.ac.za"
RUN  apt-get -y update  
RUN  apt-get -y install ncbi-blast+ git make g++ && \
     git clone https://github.com/DerrickWood/kraken2.git && \
        cd kraken2 && \
        ./install_kraken2.sh  /usr/local/bin 
RUN apt-get -y install rsync

LABEL tool=kraken2 version=latest
