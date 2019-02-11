# Dockerfile for taxonomic classification model workflow
# Eli Moss
# elimoss@stanford.edu
# February 2019

FROM neurodebian:jessie
MAINTAINER Eli Moss "elimoss@stanford.edu"

# install debian packages
RUN apt-get update
RUN apt-get install -y eatmydata
RUN eatmydata apt-get install -y wget bzip2 \
   ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 \
   git git-annex-standalone
RUN apt-get clean
#
RUN apt-get install -y g++ build-essential
#
#
# # install anaconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	-O ~/anaconda.sh && \
	bash ~/anaconda.sh -b -p /usr/local/anaconda && \
	rm ~/anaconda.sh

# # set anaconda path
ENV PATH="/usr/local/anaconda/bin:${PATH}"
#
# # requirements for taxonomic classification
RUN conda install -y -c conda-forge -c bioconda -c r \
	kraken2 krona kraken ncurses datrie r-ggplot2 r-doby r-rcolorbrewer r-scales r-plyr r-stringi
#
RUN mkdir /usr/local/anaconda/bin/taxonomy
RUN ktUpdateTaxonomy.sh
#
# #install bracken from source
RUN git clone https://github.com/jenniferlu717/Bracken.git
RUN bash Bracken/install_bracken.sh
RUN cp -r Bracken/* /usr/local/bin/
