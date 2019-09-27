FROM ubuntu:18.04
MAINTAINER Gerrit Botha "gerrit.botha@uct.ac.za"

# Followed instructions here: https://github.com/conda/conda-docker/blob/master/miniconda2/debian/Dockerfile
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 wget \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=2 \
    && conda update conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /usr/local/envs/final-report-long-reads/bin/:$PATH

RUN wget https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/long-reads/bin/final-report-long-reads.py
COPY final-report-long-reads.py /usr/local/bin/
RUN chmod +x /usr/local/bin/final-report-long-reads.py

RUN pip install pandas
