#!/usr/bin/env nextflow
//to run: nextflow [OPT: -log /path/to/log file] run sra2fastq.nf -params-file [JSON parameter file]
//not supporting filesize or run count restrictions



process SRA2FASTQ {

    tag "$accession"
    publishDir "$params.outDir/SRA_Download", mode: 'copy'
    errorStrategy "finish" //complete any processes that didn't fail

    input: 

    val accession //single accession string

    output:
    path "$accession/*.fastq.gz"
    path "$accession/*_metadata.txt"
    path "$accession/sra2fastq_temp/*", optional: true

    script: 
    //conditionally create command-line options based on non-empty parameters, for use in the command below
    def clean = params.clean != null ? "--clean True" : "" 
    def platform_restrict = params.platform_restrict != null ? "--platform_restrict $params.platform_restrict" : ""

    //invoke sra2fastq.py with those options
    """
    sra2fastq.py $accession \
    $clean \
    $platform_restrict \
    """
}


workflow {
    accessions_ch = channel.of(params.accessions)
    fastq_ch = SRA2FASTQ(accessions_ch.flatten().unique())

}