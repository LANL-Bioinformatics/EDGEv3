#!/usr/bin/env nextflow

//double-checks that any provided adapter file is in FASTA format
process adapterFileCheck {
    label "qc"
    label "small"

    input:
    path adapterFile

    output:
    stdout

    script:
    """
    isFasta.pl $adapterFile
    """
}

//main QC process. puts parameters together and runs FaQCs.
process qc {
    label "qc"
    label "medium"

    publishDir(
        path: "${settings["outDir"]}/QcReads",
        mode: 'copy'
    )

    input:
    val settings
    path fastq
    val validAdapter
    path adapter
    val avgLen

    output:
    path "QC.{1,2}.trimmed.fastq", optional:true, emit: pairedQC
    path "QC.unpaired.trimmed.fastq", optional:true, emit: unpairedQC
    path "QC_qc_report.pdf", optional: true
    path "QC.stats.txt", optional: true
    path "QC.log", emit: log

    script:
    //adjust minLength
    def min = settings["minLen"]
    if(settings["minLen"] < 1) {
        min = Math.abs(settings["minLen"] * avgLen.toInteger())
    }

    def qcSoftware = "FaQCs"

    
    def inputArg = settings["pairedFile"] ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-u $fastq"

    def adapterArg = ""
    if(adapter.name != "NO_FILE3" && validAdapter == "Yes"){
        adapterArg = "--adapter --artifactFile $adapter"
    } 

    def polyA = settings["polyA"] ? "--polyA" : ""
    def phiX = settings["filtPhiX"] ? "--phiX" : ""

    """
    $qcSoftware $inputArg \
    -q ${settings["trimQual"]} --min_L $min --avg_q ${settings["avgQual"]} \
    -n ${settings["numN"]} --lc ${settings["filtLC"]} --5end ${settings["trim5end"]} --3end ${settings["trim3end"]} \
    --split_size 1000000 -d . -t ${task.cpus} \
    $polyA \
    $adapterArg \
    $phiX
    1>QC.log 2>&1
    """
}

workflow FAQCS {
    take:
    settings
    fastq
    avgLen


    main:
 
    //adapter setup
    adapter_ch = channel.fromPath(settings["artifactFile"], checkIfExists:true)
    //checks to see if the provided adapter file is a valid FASTA
    adapterFileCheck(adapter_ch)

    //main QC process
    qc(settings, fastq, adapterFileCheck.out, adapter_ch, avgLen)
    

    trimmed = channel.empty()
    if(settings["pairedFile"]) {
        trimmed = qc.out.pairedQC
    }
    else {
        trimmed = qc.out.unpairedQC
    }
    paired = qc.out.pairedQC
    unpaired = qc.out.unpairedQC
    
    emit:
    trimmed

}