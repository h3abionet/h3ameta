FROM ubuntu:18.04
MAINTAINER Gerrit Botha "gerrit.botha@uct.ac.za"

# Got setup from here: https://hub.docker.com/r/pbelmann/krona/dockerfile

RUN apt-get -qq update && apt-get -qq -y install curl bzip2 wget make bash \
    && wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar \
    && tar -xvf KronaTools-2.7.tar \
    && rm KronaTools-2.7.tar \
    && cd KronaTools-2.7 \
    && ./install.pl

RUN KronaTools-2.7/updateTaxonomy.sh
