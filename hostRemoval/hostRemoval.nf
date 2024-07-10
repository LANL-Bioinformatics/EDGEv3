#!/usr/bin/env nextflow
//to run: nextflow run hostRemoval.nf -params-file [JSON parameter file]


process hostRemoval {
    publishDir "$params.outDir", mode: 'copy'

    input:
    path reads
    path ref

    output:
    //-fasta option does work, but files still have .fastq extensions 
    path "${params.prefix}.{1,2,unpaired}.fastq"
    path "${params.prefix}.mapping?E.log"
    path "${params.prefix}.stats.txt", emit: cleanstats
    path "host.fast{a,q}", optional: true

    script:
    
    def files = (params.inputFiles != null && params.inputFiles.size() > 1) ? "-p ${reads[0]} ${reads[1]} " : "-u $reads "

    def refFile = ref.name != "NO_FILE" ? "-ref $ref " : ""
    def prefix = params.prefix != "host_clean" ? "-prefix $params.prefix " : ""
    def fasta = params.fasta != null ? "-fasta $params.fasta " : ""
    def similarity = params.similarity != null ? "-similarity $params.similarity " : ""
    def minScore = params.bwaMemOptions != null ? "$params.bwaMemOptions " : "-T 50 "
    def ontFlag = params.fastqSource.equalsIgnoreCase("nanopore") ? "-x ont2d " : ""
    ontFlag = params.fastqSource.equalsIgnoreCase("pacbio") ? "-x pacbio " : ontFlag
    minScore = ontFlag != "" ? "-T $params.minLen " : minScore
    def bwaMemOptions = "-bwaMemOptions \"${ontFlag} ${minScore}\" "
    def cpu = params.cpus != null ? "-cpu $params.cpus " : ""
    
    """
    host_reads_removal_by_mapping.pl\
    $refFile\
    $prefix\
    $fasta\
    $cpu\
    -host \
    $bwaMemOptions\
    $files\
    -o .
    """
}

process hostRemovalStats {
    publishDir "$params.outDir", mode: 'copy'

    input:
    path stats
    path hosts

    output:
    path "hostclean.stats.txt"
    path "HostRemovalStats.pdf"

    script:
    def hosts = "-host $hosts "
    def stats = "-s $stats"

    """
    removal_stats.pl\
    $hosts\
    $stats
    """
}
workflow {
    if (params.h != null) {
        "perl host_reads_removal_by_mapping.pl -help".execute().text
    }
    else {
        "mkdir nf_assets".execute().text
        "touch nf_assets/NO_FILE".execute().text
        providedRef = channel.fromPath(params.host, checkIfExists:true) //custom error on non-existence?
        hostRemoval(channel.fromPath(params.inputFiles).collect(), providedRef.collect())
        hostRemovalStats(hostRemoval.out.cleanstats, providedRef.collect())
    }
}