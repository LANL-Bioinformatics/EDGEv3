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
include {PHAGEFINDER} from './modules/phageFinder/phageFinder.nf'
include {ANTISMASH} from './modules/runAntiSmash/runAntiSmash.nf'
include {BINNING} from './modules/contigBinning/contigBinning.nf'
include {REPORT} from './modules/report/report.nf'

workflow {
    //input specification

    //reads as input
    //allows multiple unpaired read files, or multiple paired read files, but not both
    paired = channel.empty()
    unpaired = channel.empty()
    if(params.shared.pairedFile) {
        paired = channel.fromPath([params.shared.inputFastq, params.shared.inputFastq2].transpose().flatten(), checkIfExists:true).collect()
    }
    else {
        unpaired = channel.fromPath(params.shared.inputFastq, checkIfExists:true).collect()
    }
    
    //contigs as input or pre-assembled from reads
    contigs = channel.empty()
    annContigs = channel.empty()
    if(params.shared.inputContigs != "${projectDir}/nf_assets/NO_FILE3" || params.shared.assembledContigs != "${projectDir}/nf_assets/NO_FILE3") {
        if(params.shared.inputContigs != "${projectDir}/nf_assets/NO_FILE3") {
            contigs = channel.fromPath(params.shared.inputContigs, checkIfExists:true)
        }
        else if(params.shared.assembledContigs != "${projectDir}/nf_assets/NO_FILE3") {
            contigs = channel.fromPath(params.shared.assembledContigs, checkIfExists:true)
        }
        PROCESSCONTIGS(params.shared.plus(params.assembly).plus(params.annotation).plus(params.modules).plus(params.outputLocations), contigs)
        annContigs = PROCESSCONTIGS.out.annotationContigs
    }

    //reads processing
    COUNTFASTQ(paired.ifEmpty("${projectDir}/nf_assets/NO_FILE"), unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2"))
    avgLen = COUNTFASTQ.out.avgReadLen
    counts = COUNTFASTQ.out.counts
    paired = COUNTFASTQ.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
    unpaired = COUNTFASTQ.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")

    //SRA download and processing
    if(params.modules.sra2fastq) {
        SRA2FASTQ(params.sra2fastq.plus(params.shared).plus(params.outputLocations))
        COUNTFASTQ(SRA2FASTQ.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"]), SRA2FASTQ.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2"))
        avgLen = COUNTFASTQ_SRA.out.avgReadLen
        paired = COUNTFASTQ_SRA.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
        unpaired = COUNTFASTQ_SRA.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
    }

    //QC
    qcStats = channel.empty()
    qcReport = channel.empty()
    if(params.modules.faqcs) {
        FAQCS(params.faqcs.plus(params.shared).plus(params.outputLocations), paired, unpaired,avgLen)

        paired = FAQCS.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
        unpaired = FAQCS.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
        qcStats = FAQCS.out.qcStats
        qcReport = FAQCS.out.qcReport
    }

    //Host reads removal
    hostRemovalReport = channel.empty()
    if(params.modules.hostRemoval) {

        HOSTREMOVAL(params.hostRemoval.plus(params.shared).plus(params.outputLocations),paired,unpaired)
        paired = HOSTREMOVAL.out.paired.ifEmpty(["${projectDir}/nf_assets/NO_FILE"])
        unpaired = HOSTREMOVAL.out.unpaired.ifEmpty("${projectDir}/nf_assets/NO_FILE2")
        hostRemovalReport = HOSTREMOVAL.out.hostRemovalReport
    }

    //Assembly and validation alignment
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
            ASSEMBLY(params.assembly.plus(params.shared).plus(params.annotation).plus(params.modules).plus(params.outputLocations), paired, unpaired, avgLen)
            contigs = ASSEMBLY.out.outContigs
            annContigs = ASSEMBLY.out.annotationContigs
        }
        //run validation alignment if reads were provided
        if(params.shared.inputFastq.size() != 0 && params.sra2fastq.accessions.size() == 0) {
            READSTOCONTIGS(params.r2c.plus(params.shared).plus(params.outputLocations), paired, unpaired, contigs)
            alnStats= READSTOCONTIGS.out.alnStats
            coverageTable = READSTOCONTIGS.out.covTable
            contigStatsReport = READSTOCONTIGS.out.contigStatsReport
            contigPlots = READSTOCONTIGS.out.contigPlots
            if(params.binning.binningAbundFile.endsWith("NO_FILE3")) { //user did not provide abundance file and assembly was run
                abundances = READSTOCONTIGS.out.magnitudes
            }
        }
    }

    //Reads-based taxonomic classification
    rtaReports = channel.empty()
    if(params.modules.readsTaxonomyAssignment) {
        READSTAXONOMYASSIGNMENT(params.readsTaxonomy.plus(params.shared).plus(params.faqcs).plus(params.outputLocations), paired, unpaired, avgLen)
        rtaReports = rtaReports.concat(READSTAXONOMYASSIGNMENT.out.trees, READSTAXONOMYASSIGNMENT.out.heatmaps).collect()

    }

    //Contig-based taxonomic classification
    ctaReport = channel.empty()
    if(params.modules.contigsTaxonomyAssignment) {
        CONTIGSTAXONOMYASSIGNMENT(params.contigsTaxonomy.plus(params.shared).plus(params.outputLocations), contigs, coverageTable.ifEmpty{"DNE"})
        ctaReport = CONTIGSTAXONOMYASSIGNMENT.out.ctaReport
    }

    //Annotation and PhageFinder
    antismashInput = contigs
    annStats = channel.empty()
    if(params.modules.annotation) {
        ANNOTATION(params.annotation.plus(params.shared).plus(params.outputLocations), annContigs)
        annStats = ANNOTATION.out.annStats

        if(params.modules.phageFinder && (params.annotation.taxKingdom == null || !(params.annotation.taxKingdom.equalsIgnoreCase("viruses")))) {
            PHAGEFINDER(params.shared.plus(params.outputLocations), ANNOTATION.out.gff, ANNOTATION.out.faa, ANNOTATION.out.fna)
        }

        antismashInput = ANNOTATION.out.gbk
    }

    //secondary metabolite analysis
    if(params.modules.secondaryMetaboliteAnalysis) {
        ANTISMASH(params.shared.plus(params.SMA).plus(params.outputLocations), antismashInput)
    }

    //binning
    if(params.modules.binning) {
        BINNING(params.shared.plus(params.binning).plus(params.outputLocations), contigs, abundances)
    }

    //report generation
    REPORT(
        params.shared.plus(params.modules).plus(params.outputLocations), 
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