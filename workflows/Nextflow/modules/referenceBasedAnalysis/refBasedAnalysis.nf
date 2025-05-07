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

    output:
    

    script:
    def pairedFiles =  paired.name != "NO_FILE" ? "-paired" : ""

    """
    mkdir UnmappedReads
    mkdir readsMappingToRef
    retrieve_unmapped.pl -ref $reference \
    $pairedFiles
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
    if((settings["rg2MapUnmapped"]) || (settings["r2gExtractUnmapped"])) {
        retrieveUnmappedReads(settings, reference, paired)
    }

    //emit:

}
