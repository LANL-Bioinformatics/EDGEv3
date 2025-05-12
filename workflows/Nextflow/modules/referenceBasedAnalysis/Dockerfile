# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:main AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda

RUN conda init bash && \
    . ~/.bashrc && \
    conda create --name ref2g && \
    conda activate ref2g

# install dependencies for ref2g
RUN conda install -n ref2g -c conda-forge python=3.9 r-base ghostscript
RUN conda install -n ref2g -c bioconda bowtie2=2.5.4 samtools=1.21 minimap2=2.1.1 bcftools=1.21 samclip=0.4.0 bwa=0.7.19 perl-parallel-forkmanager=2.03 mummer4=4.0.1
#use latest miccr (v0.0.3) with Chienchi's fix for tokenizing alignment file
RUN git clone --depth 1 https://github.com/chienchi/miccr.git
RUN conda install -c conda-forge conda-pack

RUN conda-pack -n ref2g -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv
COPY --from=build /miccr /miccr
RUN apt-get update && apt-get install procps -y && apt-get clean
ENV PERL5LIB=/venv/lib/perl5/core_perl
ENV PATH=/venv/bin:/miccr:/miccr/utils:$PATH

SHELL ["/bin/bash", "-c"]
CMD /bin/bash
