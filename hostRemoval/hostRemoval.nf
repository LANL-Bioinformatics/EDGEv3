#!/usr/bin/env nextflow
//to run: nextflow run hostRemoval.nf -params-file [JSON parameter file]


process hostRemoval {
    publishDir(
        path: "$params.outDir/HostRemoval",
        mode: 'copy',
    )

    input:
    path reads
    each path(ref)

    output:
    //-fasta option does work, but files still have .fastq extensions 
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.{1,2,unpaired}.fastq"
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.mapping?E.log"
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.stats.txt", emit: cleanstats
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/host.fast{a,q}", optional: true

    script:
    
    def files = (params.inputFiles != null && params.inputFiles.size() > 1) ? "-p $reads " : "-u $reads "

    def refFile = ref.name != "NO_FILE" ? "-ref $ref " : ""
    def prefix = "-prefix ${ref.name.take(ref.name.lastIndexOf('.'))}.clean "
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

    mkdir ${ref.name.take(ref.name.lastIndexOf('.'))}
    mv host.fastq ./${ref.name.take(ref.name.lastIndexOf('.'))}
    mv ${ref.name.take(ref.name.lastIndexOf('.'))}.* ./${ref.name.take(ref.name.lastIndexOf('.'))}
    """
}

process hostRemovalStats {
    //TODO: check output for multiple host removals
    publishDir "$params.outDir/HostRemoval", mode: 'copy'

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
        providedRef = channel.fromPath(params.host, checkIfExists:true)
        hostRemoval(channel.fromPath(params.inputFiles).collect(), providedRef)
        hostRemovalStats(hostRemoval.out.cleanstats, providedRef.collect())
    }
}