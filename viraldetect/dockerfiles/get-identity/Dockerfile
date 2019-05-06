FROM ubuntu:18.04
MAINTAINER Gerrit Botha "gerrit.botha@uct.ac.za"

# Followed instructions here: https://github.com/conda/conda-docker/tree/master/miniconda3/debian

RUN apt-get -qq update && apt-get -qq -y install curl bzip2 wget \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /usr/local/envs/get-identity/bin/:$PATH

RUN wget https://raw.githubusercontent.com/grbot/h3ameta/gerrit-dev/viraldetect/dockerfiles/get-identity/get-identity.py
COPY get-identity.py /usr/local/bin/
RUN chmod +x /usr/local/bin/get-identity.py
