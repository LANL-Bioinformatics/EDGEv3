#!/usr/bin/env nextflow

process checkAndAlignToReference {
    label "r2g"
    label "medium"
    publishDir(
        path: "${settings["refBasedOutDir"]}",
        mode: 'copy'
    )

    input:
    val settings
    path reference
    val platform
    path paired
    path unpaired

    output:
    path "*"
    path "readsToRef.gaps", emit: gaps
    path "readsToRef.vcf", optional:true, emit: vcf
    path "*.sort.bam", emit: bam

    script:
    def taxKingdom = settings["taxKingdom"] != null ? "-kingdom ${settings["taxKingdom"]}" : ""
    def pairedFiles = paired.name != "NO_FILE" ? "-p \"$paired\"" : ""
    def unpairedFiles = unpaired.name != "NO_FILE2" ? "-u $unpaired" : ""
    def platformArg = platform != null ? "-plat $platform" : ""
    def alnOptions = settings["r2gAlignerOptions"] != "" ? "-alnOpt ${settings["r2gAlignerOptions"]}" : ""
    def minMapQual = settings["r2gMinMapQual"] != null ? "-minmap ${settings["r2gMinMapQual"]}" : ""
    def maxClip = settings["r2gMaxClip"] != null ? "-maxclip ${settings["r2gMaxClip"]}" : ""
    def extractMapped = settings["r2gExtractMapped"] != false ? "-x-mapped ${settings["r2gExtractMapped"]}" : "" 
    
    //variant call configurations
    def variantCall = settings["r2gVariantCall"] != false ? "-vc 1" : ""
    def vcQual = settings["r2gVariantCallMinQual"] != null ? "-vc-qual ${settings["r2gVariantCallMinQual"]}" : ""
    def vcPloidy = settings["r2gVariantCallPloidy"] != null ? "-vc-ploidy ${settings["r2gVariantCallPloidy"]}" : ""
    

    //consensus sequence configurations
    def doConsensus = settings["r2gGetConsensus"] != false ? "-consensus 1" : ""
    def consensusMapQual = settings["r2gConsensusMinMapQual"] != null ? "-c-mapq ${settings["r2gConsensusMinMapQual"]}" : ""
    def consensusMinCov = settings["r2gConsensusMinCov"] != null ? "-c-mincov ${settings["r2gConsensusMinCov"]}" : ""
    def consensusMaxCov = settings["r2gConsensusMaxCov"] != null ? "-c-maxcov ${settings["r2gConsensusMaxCov"]}" : ""
    def consensusAltProp = settings["r2gConsensusAltProp"] != null ? "-c-altprop ${settings["r2gConsensusAltProp"]}" : ""
    def consensusAltIndelProp = settings["r2gConsensusAltIndelProp"] != null ? "-c-indelprop ${settings["r2gConsensusAltProp"]}" : ""
    def consensusMinBaseQual = settings["r2gConsensusMinBaseQ"] != null ? "-c-baseq ${settings["r2gConsensusMinBaseQ"]}" : ""
    def consensusDisableBAQ = settings["r2gConsensusDisableBAQ"] != false ? "-c-baq 0" : ""
    def consensusPCRdedup = settings["r2gConsensusPCRdedup"] != false ? "-c-dedup 1" : "" 
    def consensusHomopolymerFilt = settings["r2gConsensusHomopolymerFilt"] != false ? "-c-polymer 1" : ""
    def consensusStrandBiasFilt = settings["r2gConsensusStrandBiasFilt"] != false ? "-c-sb 1" : ""
    def consensusVarlogOpt = settings["r2gConsensusVarlogOpt"] != false ? "-c-varlog 1" : ""
    def consensusCompOpt = settings["r2gConsensusCompOpt"] != false ? "-c-compopt 1" : ""

    """
    ref_pipeline.pl -ref $reference \
    $pairedFiles \
    $unpairedFiles \
    -t ${task.cpus} \
    -proj ${settings["projName"]} \
    -out \$PWD \
    $taxKingdom \
    $platformArg \
    -aln ${settings["r2gAligner"]} \
    $alnOptions \
    $minMapQual \
    $maxClip \
    $variantCall \
    $vcQual \
    $vcPloidy \
    $extractMapped \
    $doConsensus \
    $consensusMapQual \
    $consensusMinCov \
    $consensusMaxCov \
    $consensusAltProp \
    $consensusAltIndelProp \
    $consensusMinBaseQual \
    $consensusDisableBAQ \
    $consensusPCRdedup \
    $consensusHomopolymerFilt \
    $consensusStrandBiasFilt \
    $consensusVarlogOpt \
    $consensusCompOpt
    """

    
}

process retrieveUnmappedReads {
    label "r2g"
    label "small"

    input:
    val settings
    path reference
    path paired
    path bams //staged into workdir for retrieve_unmapped.pl

    output:
    stdout emit: count
    path "singleEnd.fastq", emit: singleUnmapped, optional: true
    path "pairedEnd.*.fastq", emit: pairedUnmapped, optional: true
    

    script:
    def pairedFiles =  paired.name != "NO_FILE" ? "-paired" : ""

    """
    mkdir UnmappedReads
    retrieve_unmapped.pl \
    -ref $reference \
    $pairedFiles
    """

}

process mapUnmapped {
    label "r2g"
    label "medium"
    input:
    val settings
    path unmappedPaired
    path unmappedUnpaired
    val count
    val platform

    when:
    !count.contains("Total Unmapped:0")

    output:
    
    script:
    def ontFlag = (platform != null && platform.contains("NANOPORE")) ?  "-x ont2d -T ${settings["minLen"] != null ? settings["minLen"] : 50} " : ""
    def pbFlag =  (platform != null && platform.contains("PACBIO")) ? "-x pacbio -T ${settings["minLen"] != null ? settings["minLen"] : 50} " : ""
    def pairedFiles = unmappedPaired.name != "NO_FILE" ? "-p \"$unmappedPaired\"" : ""
    def unpairedFiles = unmappedUnpaired.name != "NO_FILE2" ? "-u $unmappedUnpaired" : ""
    def maxClip = settings["r2gMaxClip"] != null ? "-max_clip ${settings["r2gMaxClip"]}" : ""

    """
    runReadsToContig.pl \
    -c 0 \
    -cpu ${task.cpus} \
    $pairedFiles \
    $unpairedFiles \
    $maxClip \
    -bwa_options \'$ontFlag $pbFlag\' \
    -d . -pre UnmappedReads -ref ${settings["refSeqDB"]} \
    &>mapping.log
    """
}


workflow REFERENCEBASEDANALYSIS {

    take:
    settings
    platform
    paired
    unpaired

    main:

    reference = channel.fromPath(settings["referenceGenome"], checkIfExists: true)
    checkAndAlignToReference(settings, reference, platform, paired, unpaired)
    if((settings["r2gMapUnmapped"].toBoolean())|| (settings["r2gExtractUnmapped"].toBoolean())) {
        retrieveUnmappedReads(settings, reference, paired, checkAndAlignToReference.out.bam)
        if(settings["r2gMapUnmapped"]) {
            mapUnmapped(settings, 
                retrieveUnmappedReads.out.pairedUnmapped.ifEmpty(["${projectDir}/nf_assets/NO_FILE"]),  
                retrieveUnmappedReads.out.singleUnmapped.ifEmpty("${projectDir}/nf_assets/NO_FILE2"), 
                retrieveUnmappedReads.out.count,
                platform)
        }
    }

    //emit:

}
