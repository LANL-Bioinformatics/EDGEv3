#!/usr/bin/env nextflow


//double-checks that any provided adapter file is in FASTA format
process adapterFileCheck {
    label "qc"
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
//EDGE currently uses a custom script (illumina_fastq_QC.pl) to handle QC for long reads,
//but it was unable to create report files when I attempted using it. For now, all input reads go through FaQCs.
process qc {
    label "qc"
    publishDir(
        path: "${settings["qcOutDir"]}",
        mode: 'copy'
    )

    input:
    val settings
    path paired
    path unpaired
    val validAdapter
    path adapter
    val avgLen

    output:
    path "QC.{1,2}.trimmed.fastq", optional:true, emit: pairedQC
    path "QC.unpaired.trimmed.fastq", optional:true, emit: unpairedQC
    path "QC_qc_report.pdf", optional: true, emit: qcReport
    path "QC.stats.txt", optional: true, emit: qcStats
    path "QC.log", emit: log

    script:
    //adjust minLength
    min = settings["minLen"]
    if(settings["minLen"] < 1) {
        min = Math.abs(settings["minLen"] * avgLen.toInteger())
    }

    def qcSoftware = "FaQCs"
    // if(params.ontFlag || params.pacbioFlag) {
    //     qcSoftware = "illumina_fastq_QC.pl"
    // }
    def pairedArg = paired.name != "NO_FILE" ? "-1 ${paired[0]} -2 ${paired[1]}" : ""
    // if(pairedArg != "" && (params.ontFlag || params.pacbioFlag)) {
    //     pairedArg = "-p $paired"
    // }
    def unpairedArg = unpaired.name != "NO_FILE2" ? "-u $unpaired" : ""

    def adapterArg = ""
    if(adapter.name != "NO_FILE3" && validAdapter == "Yes"){
        adapterArg = "--adapter --artifactFile $adapter"
    } 

    polyA = settings["trimPolyA"] ? "--polyA" : ""
    phiX = settings["filtPhiX"] ? "--phiX" : ""

    def trim = ""
    // if(params.ontFlag || params.pacbioFlag) {
    //     trim = "--trim_only"
    // }

    """
    $qcSoftware $pairedArg $unpairedArg \
    -q ${settings["trimQual"]} --min_L $min --avg_q ${settings["avgQual"]} \
    -n ${settings["numN"]} --lc ${settings["filtLC"]} --5end ${settings["trim5end"]} --3end ${settings["trim3end"]} \
    --split_size 1000000  -d . -t ${settings["cpus"]} \
    $polyA \
    $trim \
    $adapterArg \
    $phiX \
    1>QC.log 2>&1
    """
}

// process nanoplot {
//     label "qc"
//     publishDir(
//         path: "${settings["outDir"]}/QcReads",
//         mode: 'copy'
//     )
//     input:
//     val settings
//     path unpaired

//     output:
//     path "*" //lots of output plots

//     script:
//     """
//     NanoPlot --fastq $unpaired --N50 --loglength -t ${settings["cpus"]} -f pdf --outdir . 2>/dev/null
//     """

// }

// //Porechop for removing adapters from ONT or PacBio reads
// process porechop {
//     label "qc"
//     publishDir(
//         path: "${settings["outDir"]}/QcReads",
//         mode: 'copy'
//     )


//     input:
//     val settings
//     path trimmed
//     path log
//     output:
//     path "*.porechop.fastq", emit: porechopped
    
//     script:
//     """
//     porechop -i $trimmed -o ./QC.unpaired.porechop.fastq -t ${settings["cpus"]} > $log
//     """
// }

workflow FAQCS {
    take:
    settings
    paired
    unpaired
    avgLen


    main:
 
    //adapter setup
    adapter_ch = channel.fromPath(settings["artifactFile"], checkIfExists:true)
    //checks to see if the provided adapter file is a valid FASTA
    adapterFileCheck(adapter_ch)

    //main QC process
    qc(settings, paired, unpaired, adapterFileCheck.out, adapter_ch, avgLen)
    

    paired = qc.out.pairedQC
    unpaired = qc.out.unpairedQC
    qcReport = qc.out.qcReport
    qcStats = qc.out.qcStats
    
    emit:
    paired
    unpaired
    qcReport
    qcStats

}