#!/usr/bin/env nextflow

include {SRA2FASTQ} from './modules/sra2fastq/sra2fastq.nf'
include {COUNTFASTQ} from './modules/countFastq/countFastq.nf'
include {COUNTFASTQ_SRA} from './modules/countFastq/countFastq.nf'
include {PROCESSCONTIGS} from './modules/processProvidedContigs/processProvidedContigs.nf'
include {FAQCS} from './modules/runFaQCs/runFaQCs.nf'
include {HOSTREMOVAL} from './modules/hostRemoval/hostRemoval.nf'
include {ASSEMBLY} from './modules/runAssembly/runAssembly.nf'
include {READSTOCONTIGS} from './modules/runReadsToContig/runReadsToContig.nf'
include {READSTAXONOMYASSIGNMENT} from './modules/readsTaxonomyAssignment/readsTaxonomyAssignment.nf'
include {CONTIGSTAXONOMYASSIGNMENT} from './modules/contigsTaxonomyAssignment/contigsTaxonomyAssignment.nf'
include {ANNOTATION} from './modules/runAnnotation/runAnnotation.nf'
include {PHAGEFINDER} from './modules/phageFinder/phageFinder.nf'
include {ANTISMASH} from './modules/runAntiSmash/runAntiSmash.nf'
include {BINNING} from './modules/readsBinning/readsBinning.nf'
include {REPORT} from './modules/report/report.nf'

workflow {

    //input specification
    fastqFiles = channel.empty()

    if(params.shared.inputFastq.size() != 0) {
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

    COUNTFASTQ(params.shared, fastqFiles.collect())

    avgLen = COUNTFASTQ.out.avgReadLen
    counts = COUNTFASTQ.out.counts
    paired = COUNTFASTQ.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
    unpaired = COUNTFASTQ.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
    if(params.modules.sra2fastq) {
        SRA2FASTQ(params.sra2fastq.plus(params.shared))
        COUNTFASTQ_SRA(SRA2FASTQ.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"]), SRA2FASTQ.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2"))
        avgLen = COUNTFASTQ_SRA.out.avgReadLen
        paired = COUNTFASTQ_SRA.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
        unpaired = COUNTFASTQ_SRA.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
    }

    qcStats = channel.empty()
    qcReport = channel.empty()
    if(params.modules.faqcs) {
        FAQCS(params.faqcs.plus(params.shared), paired, unpaired,avgLen)

        paired = FAQCS.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
        unpaired = FAQCS.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
        qcStats = FAQCS.out.qcStats
        qcReport = FAQCS.out.qcReport
    }

    hostRemovalReport = channel.empty()
    if(params.modules.hostRemoval) {

        HOSTREMOVAL(params.hostRemoval.plus(params.shared),paired,unpaired)
        paired = HOSTREMOVAL.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
        unpaired = HOSTREMOVAL.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
        hostRemovalReport = HOSTREMOVAL.out.hostRemovalReport
    }

    coverageTable = channel.empty()
    abundances = channel.empty()
    contigStatsReport = channel.empty()
    contigPlots = channel.empty()
    alnStats = channel.empty()
    if(!params.binning.binningAbundFile.endsWith("NO_FILE3")) { //user provided abundance file
        abundances = channel.fromPath(params.binning.binningAbundFile, checkIfExists:true)
    }
    if(params.modules.runAssembly) {
        //assemble if not already using assembled or provided contigs
        if (params.shared.inputContigs == "${projectDir}/nf_assets/NO_FILE3" && params.shared.assembledContigs == "${projectDir}/nf_assets/NO_FILE3") {
            ASSEMBLY(params.assembly.plus(params.shared).plus(params.annotation).plus(params.modules), paired, unpaired, avgLen)
            contigs = ASSEMBLY.out.outContigs
            annContigs = ASSEMBLY.out.annotationContigs
        }
        //run validation alignment if reads were provided
        if(params.shared.inputFastq.size() != 0 && params.sra2fastq.accessions.size() == 0) {
            READSTOCONTIGS(params.r2c.plus(params.shared), paired, unpaired, contigs)
            alnStats= READSTOCONTIGS.out.alnStats
            coverageTable = READSTOCONTIGS.out.covTable
            contigStatsReport = READSTOCONTIGS.out.contigStatsReport
            contigPlots = READSTOCONTIGS.out.contigPlots
            if(params.binning.binningAbundFile.endsWith("NO_FILE3")) { //user did not provide abundance file and assembly was run
                abundances = READSTOCONTIGS.out.magnitudes
            }
        }
    }

    rtaReports = channel.empty()
    if(params.modules.readsTaxonomyAssignment) {
        READSTAXONOMYASSIGNMENT(params.readsTaxonomy.plus(params.shared).plus(params.faqcs), paired, unpaired, avgLen)
        rtaReports = rtaReports.concat(READSTAXONOMYASSIGNMENT.out.trees, READSTAXONOMYASSIGNMENT.out.heatmaps).collect()
    }

    ctaReport = channel.empty()
    if(params.modules.contigsTaxonomyAssignment) {
        CONTIGSTAXONOMYASSIGNMENT(params.contigsTaxonomy.plus(params.shared), contigs, coverageTable.ifEmpty{"DNE"})
        ctaReport = CONTIGSTAXONOMYASSIGNMENT.out.ctaReport
    }

    antismashInput = contigs
    annStats = channel.empty()
    if(params.modules.annotation) {
        ANNOTATION(params.annotation.plus(params.shared), annContigs)
        annStats = ANNOTATION.out.annStats

        if(params.modules.phageFinder && (params.annotation.taxKingdom == null || !(params.annotation.taxKingdom.equalsIgnoreCase("viruses")))) {
            PHAGEFINDER(params.shared, ANNOTATION.out.gff, ANNOTATION.out.faa, ANNOTATION.out.fna)
        }

        antismashInput = ANNOTATION.out.gbk
    }

    if(params.modules.secondaryMetaboliteAnalysis) {
        ANTISMASH(params.shared.plus(params.SMA), antismashInput)
    }

    if(params.modules.readsBinning) {
        BINNING(params.shared.plus(params.binning), contigs, abundances)
    }
    //TODO: channel.empty() parameters here indicate files from upstream processes not yet implemented into report generation
    REPORT(
        params.shared.plus(params.modules), 
        counts.ifEmpty{file("DNE")},
        qcStats.ifEmpty{file("DNE1")},
        qcReport.ifEmpty{file("DNE2")},
        hostRemovalReport.ifEmpty{file("DNE3")},
        rtaReports.ifEmpty{file("DNE4")}, 
        ctaReport.ifEmpty{file("DNE5")},
        contigStatsReport.ifEmpty{file("DNE6")},
        contigPlots.ifEmpty{file("DNE7")},
        annStats.ifEmpty{file("DNE8")},
        alnStats.ifEmpty{file("DNE9")}
    )
}