#!/usr/bin/env nextflow
//to run: nextflow [OPT: -log /path/to/log file] run sra2fastq.nf -params-file [JSON parameter file]
//not supporting filesize or run count restrictions


accessions_ch = Channel.of(params.accessions)

process SRA2FASTQ {
    //debug true

    tag "$accession"
    publishDir "$params.outdir", mode: 'copy'
    errorStrategy "finish" //complete any processes that didn't fail

    input: 

    val accession //single accession string

    output:
    path "$accession/*.fastq.gz"
    path "$accession/*_metadata.txt"
    path "$accession/sra2fastq_temp/*", optional: true

    script: 
    //conditionally create command-line options based on non-empty parameters, for use in the command below
    def clean = params.clean != "" ? "--clean True" : "" 
    def platform_restrict = params.platform_restrict != "" ? "--platform_restrict $params.platform_restrict" : ""

    //invoke sra2fastq.py with those options
    """
    sra2fastq.py $accession \
    $clean \
    $platform_restrict \
    """
}


workflow {

    fastq_ch = SRA2FASTQ(accessions_ch.flatten().unique())

}