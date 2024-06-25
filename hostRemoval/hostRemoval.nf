#!/usr/bin/env nextflow
//to run: nextflow run hostRemoval.nf -params-file [JSON parameter file]

//overwriteable parameters
params.inputFiles = null
params.ref = null
params.outDir = null
params.prefix = "host_clean" //perl host_removal script default
params.fasta = null
params.host = null
params.bwaMemOptions = null
params.cpu = null
params.h = null

process hostRemoval {
    publishDir ".", mode:copy

    input:
    path "reads"

    output:
    //-fasta option does work, but files still have .fastq extensions 
    path "$params.outdir/${params.prefix}.{1,2,unpaired}.fast{a,q}"
    path "$params.outdir/${params.prefix}.mappingPE.log"
    path "$params.outdir/${params.prefix}.stats.txt"
    path "$params.outdir/host.fast{a,q}", optional: true

    def outDir = params.outDir != null ? "-o ./$params.outDir " : ""
    def files = (params.inputFiles != null && params.inputFiles.size() > 1) ? "-p reads1 reads2 " : "-u reads "

    def ref = params.ref != null ? "-ref ./$params.ref " : ""
    def prefix = params.prefix != "host_clean" ? "-prefix ./$params.prefix " : ""
    def fasta = params.fasta != null ? "-fasta ./$params.fasta " : ""
    def host = params.host != null ? "-host ./$params.host " : ""
    def bwaMemOptions = params.bwaMemOptions != null ? "-d ./$params.bwaMemOptions " : ""
    def cpu = params.cpu != null ? "-d ./$params.cpu " : ""
    
    script:
    """
    perl host_reads_removal_by_mapping.pl\
$ref\
$prefix\
$fasta\
$host\
$bwaMemOptions\
$files\
$outDir
    """

}

workflow {
    if params.h != null {
        "perl host_reads_removal_by_mapping.pl -help".execute().text
    }
    else {
        hostRemoval(channel.fromPath(params.inputFiles, relative:true).collect())
    }
}