# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:main AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name binning \
    && conda activate binning 

# install dependencies for phageFinder
RUN conda install -n binning -c bioconda maxbin2=2.2.7
RUN conda install -c conda-forge conda-pack

RUN conda-pack -n binning -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash