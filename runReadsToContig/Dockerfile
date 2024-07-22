# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

# install dependencies
RUN conda install -c conda-forge r-base
RUN conda install -c conda-forge python=3.11
RUN conda install -c bioconda perl-json
RUN conda install -c bioconda samclip=0.4.0
RUN conda install -c bioconda bwa
RUN conda install -c bioconda bowtie2
RUN conda install -c bioconda samtools=1.6

ADD bin/bam_to_fastq.pl /opt/conda/bin
ADD bin/contig_stats.pl /opt/conda/bin
ADD bin/ContigCoverageFold_plots_from_samPileup.pl /opt/conda/bin
ADD bin/fastq_utility.pm /opt/conda/bin
ADD bin/runReadsToContig.pl /opt/conda/bin
ADD bin/tab2Json_for_dataTable.pl /opt/conda/bin

CMD ["/bin/bash"]