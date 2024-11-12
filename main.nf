#!/usr/bin/env nextflow

include {SRA2FASTQ} from './modules/sra2fastq/sra2fastq.nf'
include {COUNTFASTQ} from './modules/countFastq/countFastq.nf'
include {FAQCS} from './modules/runFaQCs/runFaQCs.nf'

workflow {

    //input specification

    fastqFiles = channel.fromPath(params.shared.inputFastq, checkIfExists:true)

    if(params.modules.sra2fastq) {
        SRA2FASTQ(params.sra2fastq.plus(params.shared))
        fastqFiles = fastqFiles.concat(SRA2FASTQ.out.fastq).flatten()
    }
    
    COUNTFASTQ(params.shared, fastqFiles.collect())

    avgLen = COUNTFASTQ.out.avgReadLen
    fastqFiles = COUNTFASTQ.out.fastqFiles


    if(params.modules.faqcs) {
        FAQCS(params.faqcs.plus(params.shared), fastqFiles,avgLen)
    }

}