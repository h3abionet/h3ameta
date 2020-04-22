#############################################################
# Dockerfile to build a container for pairfq
#############################################################

FROM ubuntu:18.04
MAINTAINER Gerrit Botha "gerrit.botha@uct.ac.za"

RUN apt-get -m update && apt-get install -y wget unzip zip curl perl-base cpanminus \
libc-dev build-essential

RUN mkdir pairfq
WORKDIR /pairfq
RUN curl -sL git.io/pairfq_lite > pairfq_lite

RUN chmod +x pairfq_lite

# Switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

WORKDIR /
ENV PATH /pairfq:$PATH
