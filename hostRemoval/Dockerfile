# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

# install dependencies for hostRemoval
RUN conda install -c bioconda bwa=0.7.17
RUN conda install -c conda-forge r-base=4.4.0

ADD bin/host_reads_removal_by_mapping.pl /opt/conda/bin
ADD bin/removal_stats.pl /opt/conda/bin

CMD ["/bin/bash"]