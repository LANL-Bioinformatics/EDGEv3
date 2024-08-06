#!/usr/bin/env nextflow
//to run: nextflow run hostRemoval.nf -params-file [JSON parameter file]


//clean input FASTQ files of reads that map to provided host references, one process per given host
process hostRemoval {
    publishDir(
        path: "$params.outDir/HostRemoval",
        mode: 'copy'
    )
    
    tag "${ref.name.take(ref.name.lastIndexOf('.'))}"

    input:
    path reads
    each path(ref)

    output:
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.1.fastq", emit: cleaned1
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.2.fastq", emit: cleaned2
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.unpaired.fastq", emit: cleanedSingleton
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.mapping?E.log"
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/${ref.name.take(ref.name.lastIndexOf('.'))}.clean.stats.txt", emit: cleanstats
    path "${ref.name.take(ref.name.lastIndexOf('.'))}/host.fastq"

    script:
    
    def files = (params.inputFiles != null && params.inputFiles.size() > 1) ? "-p $reads " : "-u $reads "

    def refFile = ref.name != "NO_FILE" ? "-ref $ref " : ""
    def prefix = "-prefix ${ref.name.take(ref.name.lastIndexOf('.'))}.clean "
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

//merge cleaned FASTQ files into files cleaned of ALL reads mapping to ANY provided host reference
process mergeCleaned {
    publishDir(
        path: "$params.outDir/HostRemoval",
        mode: 'copy'
    )

    tag "${cleanedFiles[0].name.tokenize('.')[-2]}"

    input:
    path cleanedFiles

    output:
    path "hostclean.{1,2}.fastq"

    
    script:
    """
    seqkit common $cleanedFiles -n > hostclean.${cleanedFiles[0].name.tokenize('.')[-2]}.fastq
    """
}

//Concatenate leftover unpaired reads that didn't map to a host reference, and remove any name duplicates (i.e., leftovers appearing in multiple cleanings)
process mergeCleanUnpaired {
    publishDir(
        path: "$params.outDir/HostRemoval",
        mode: 'copy'
    )
    
    input:
    path remainingUnpairedReads

    output:
    path "hostclean.unpaired.fastq"

    script:
    """
    cat $remainingUnpairedReads > hostclean.unpaired.fastq
    seqkit rmdup -n hostclean.unpaired.fastq
    """
}

process hostRemovalStats {
    publishDir "$params.outDir/HostRemoval", mode: 'copy'

    input:
    path stats

    output:
    path "hostclean.stats.txt" //issue with counting here
    path "HostRemovalStats.pdf"

    script:
    def stats = "-stats $stats "

    """
    removal_stats.pl\
    $stats\
    """
}


workflow {
    if (params.h != null) {
        //help option for host removal script <- necessary?
        "perl host_reads_removal_by_mapping.pl -help".execute().text
    }
    else {

        //setup
        "mkdir nf_assets".execute().text
        "touch nf_assets/NO_FILE".execute().text
        providedRef = channel.fromPath(params.host, checkIfExists:true)

        //remove host reads in parallel
        hostRemoval(channel.fromPath(params.inputFiles).collect(), providedRef.collect())

        //merge clean paired-end reads (intersection)
        mergeCleaned(hostRemoval.out.cleaned1.collect()
            .concat(
                hostRemoval.out.cleaned2.collect(), 
            )
        )

        //merge clean unpaired reads (removing any duplicates by read name)
        mergeCleanUnpaired(hostRemoval.out.cleanedSingleton.collect())

        //calculate overall stats and create PDF
        hostRemovalStats(hostRemoval.out.cleanstats.collect())
    }
}