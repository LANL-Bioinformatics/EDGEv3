#!/usr/bin/env nextflow

include {SRA2FASTQ} from './modules/sra2fastq/sra2fastq.nf'
include {COUNTFASTQ} from './modules/countFastq/countFastq.nf'
include {FAQCS} from './modules/runFaQCs/runFaQCs.nf'
include {HOSTREMOVAL} from './modules/hostRemoval/hostRemoval.nf'
include {ASSEMBLY} from './modules/runAssembly/runAssembly.nf'
include {READSTOCONTIGS} from './modules/runReadsToContig/runReadsToContig.nf'
include {READSTAXONOMYASSIGNMENT} from './modules/readsTaxonomyAssignment/readsTaxonomyAssignment.nf'
include {CONTIGSTAXONOMYASSIGNMENT} from './modules/contigsTaxonomyAssignment/contigsTaxonomyAssignment.nf'

workflow {

    //input specification    
    fastqFiles = channel.fromPath(params.shared.inputFastq, checkIfExists:true)
    contigs = channel.empty()
    if(params.r2c.useAssembledContigs) {
        contigs = channel.fromPath(params.shared.inputContigs, checkIfExists:true)
    }


    if(params.modules.sra2fastq) {
        SRA2FASTQ(params.sra2fastq.plus(params.shared))
        fastqFiles = fastqFiles.concat(SRA2FASTQ.out.fastq).flatten()
    }
    
    COUNTFASTQ(params.shared, fastqFiles.collect())

    avgLen = COUNTFASTQ.out.avgReadLen
    paired = COUNTFASTQ.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE","${projectDir}/nf_assets/NO_FILE2"])
    unpaired = COUNTFASTQ.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE")


    if(params.modules.faqcs) {
        FAQCS(params.faqcs.plus(params.shared), paired, unpaired,avgLen)

        paired = FAQCS.out.paired
        unpaired = FAQCS.out.unpaired
    }

    if(params.modules.hostRemoval) {
        HOSTREMOVAL(params.hostRemoval.plus(params.shared),paired,unpaired)
        paired = HOSTREMOVAL.out.paired.ifEmpty(params.pairedFiles)
        unpaired = HOSTREMOVAL.out.unpaired.ifEmpty(params.unpairedFiles)
    }

    if(params.modules.runAssembly && !params.r2c.useAssembledContigs) {
        ASSEMBLY(params.assembly.plus(params.shared), paired, unpaired, avgLen)
        contigs = ASSEMBLY.out.outContigs
        READSTOCONTIGS(params.r2c.plus(params.shared), paired, unpaired, contigs)
    }

    //should always run if contigs were provided or generated
    READSTOCONTIGS(params.r2c.plus(params.shared), paired, unpaired, contigs)

    if(params.modules.readsTaxonomyAssignment) {
        READSTAXONOMYASSIGNMENT(params.readsTaxonomy.plus(params.shared).plus(params.faqcs), paired, unpaired, avgLen)
    }

    if(params.modules.contigsTaxonomyAssignment) {
        CONTIGSTAXONOMYASSIGNMENT(params.contigsTaxonomy.plus(params.shared), contigs, READSTOCONTIGS.out.covTable)
    }

}