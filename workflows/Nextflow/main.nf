#!/usr/bin/env nextflow

include {SRA2FASTQ} from './modules/sra2fastq/sra2fastq.nf'
include {COUNTFASTQ} from './modules/countFastq/countFastq.nf'
include {PROCESSCONTIGS} from './modules/processProvidedContigs/processProvidedContigs.nf'
include {FAQCS} from './modules/runFaQCs/runFaQCs.nf'
include {HOSTREMOVAL} from './modules/hostRemoval/hostRemoval.nf'
include {ASSEMBLY} from './modules/runAssembly/runAssembly.nf'
include {READSTOCONTIGS} from './modules/runReadsToContig/runReadsToContig.nf'
include {READSTAXONOMYASSIGNMENT} from './modules/readsTaxonomyAssignment/readsTaxonomyAssignment.nf'
include {CONTIGSTAXONOMYASSIGNMENT} from './modules/contigsTaxonomyAssignment/contigsTaxonomyAssignment.nf'
include {ANNOTATION} from './modules/runAnnotation/runAnnotation.nf'

workflow {

    //input specification
    fastqFiles = channel.empty()
    if(params.shared.inputFastq != null) {
        fastqFiles = channel.fromPath(params.shared.inputFastq, checkIfExists:true)
    }
    
    contigs = channel.empty()
    annContigs = channel.empty()
    if(params.shared.inputContigs != "${projectDir}/nf_assets/NO_FILE3" || params.shared.assembledContigs != "${projectDir}/nf_assets/NO_FILE3") {
        if(params.shared.inputContigs != "${projectDir}/nf_assets/NO_FILE3") {
            contigs = channel.fromPath(params.shared.inputContigs, checkIfExists:true)
        }
        else if(params.shared.assembledContigs != "${projectDir}/nf_assets/NO_FILE3") {
            contigs = channel.fromPath(params.shared.assembledContigs, checkIfExists:true)
        }
        PROCESSCONTIGS(params.shared.plus(params.assembly).plus(params.annotation).plus(params.modules), contigs)
        annContigs = PROCESSCONTIGS.out.annotationContigs
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

    coverageTable = channel.empty()
    if(params.modules.runAssembly) {
        //assemble if not already using assembled or provided contigs
        if (params.inputContigs == "${projectDir}/nf_assets/NO_FILE3" && params.assembledContigs == "${projectDir}/nf_assets/NO_FILE3") {
            ASSEMBLY(params.assembly.plus(params.shared).plus(params.annotation).plus(params.modules), paired, unpaired, avgLen)
            contigs = ASSEMBLY.out.outContigs
            annContigs = ASSEMBLY.out.annotationContigs
        }
        //run validation alignment if reads were provided
        READSTOCONTIGS(params.r2c.plus(params.shared), paired, unpaired, contigs)
        coverageTable = READSTOCONTIGS.out.covTable
    }


    if(params.modules.readsTaxonomyAssignment) {
        READSTAXONOMYASSIGNMENT(params.readsTaxonomy.plus(params.shared).plus(params.faqcs), paired, unpaired, avgLen)
    }

    if(params.modules.contigsTaxonomyAssignment) {
        CONTIGSTAXONOMYASSIGNMENT(params.contigsTaxonomy.plus(params.shared), contigs, coverageTable.ifEmpty{"DNE"})
    }

    if(params.modules.annotation) {
        ANNOTATION(params.annotation.plus(params.shared), annContigs)
    }
}