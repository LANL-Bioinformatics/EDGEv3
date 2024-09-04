#!/usr/bin/env nextflow

process unsetKingdom {
    input:
    path contigs

    output:
    stdout

    shell:
    '''
    a=$(grep -v ">" !{contigs} | wc -m)
    b=$(grep -v ">" !{contigs} | wc -l)
    c=$((a-b))
    if [ c -gt 580000 ]; then echo "Bacteria"; else echo "Viruses"; fi
    '''
}

process prokkaAnnotate {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )

    input:
    path contigs
    path protein
    path hmm
    val kingdom

    output:
    path "*"

    script:
    def kingdom = kingdom.trim()
    def protein = protein.name == "NO_FILE" ? "" : "--protein $params.customProtein"
    def hmmPrep = hmm.name == "NO_FILE2" ? "" : "hmmpress $hmm" //TODO: test functionality
    def hmm = hmm.name == "NO_FILE2" ? "" : "--hmms $hmm"
    def evalue = params.evalue == null ? "" : "--evalue $params.evalue"
    def gcode = params.gcode == null ? "" : "--gcode $params.gcode"
    def locustag = params.projName == null ? "" : "--locustag $params.projName"
    def prefix = params.projName == null ? "" : "--prefix $params.projName"
    def cpu = params.cpus == null ? "" : "--cpus $params.cpus"
    def taxKingdom = kingdom.equalsIgnoreCase("metagenome") ? "--kingdom Bacteria --metagenome" : "--kingdom $kingdom"

    """
    $hmmPrep

    prokka --quiet --force \
    $protein \
    $hmm \
    $evalue \
    $gcode \
    $locustag \
    $prefix \
    $cpu \
    --outdir . \
    $taxKingdom \
    $contigs 2>>Annotation.log 
    """
}

process rattAnnotate {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )

    input:
    path contigs
    path gbk

    output:
    path "*"

    shell:
    //happens in work directory
    //RATT version needs to be custom EDGE version
    '''
    mkdir -p ./RATT/source
    cp !{gbk} ./RATT/source/source.gbk
    genbank2embl.pl ./RATT/source/source.gbk
    runRATT.sh $PWD/RATT/source !{contigs} !{params.projName} Species 1>>Annotation.log
    '''
}


workflow {

    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text

    contig_ch = channel.fromPath(params.inputContigs, checkIfExists:true)
    kingdom_ch = channel.of(params.taxKingdom)
    hmm_ch = channel.fromPath(params.customHMM, checkIfExists:true)
    prot_ch = channel.fromPath(params.customProtein, checkIfExists:true)
    if (params.taxKingdom == null) {
        kingdom_ch = unsetKingdom(contig_ch)
    }
    //TODO: add safety for bad source reference GenBank file
    if (params.annotateProgram =~ /(?i)prokka/) {

        prokkaAnnotate(contig_ch, prot_ch, hmm_ch, kingdom_ch)
    }
    else if (params.annotateProgram =~ /(?i)ratt/) {
        gb_ch = channel.fromPath(params.sourceGBK, checkIfExists:true)
        rattAnnotate(contig_ch, gb_ch)
    }
}