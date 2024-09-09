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
    path "Annotation.log"
    path "${params.projName}.gff", emit: gff
    path "${params.projName}.err"
    path "${params.projName}.faa"
    path "${params.projName}.ffn"
    path "${params.projName}.fna"
    path "${params.projName}.fsa"
    path "${params.projName}.gbk"
    path "${params.projName}.log"
    path "${params.projName}.sqn"
    path "${params.projName}.tbl"
    path "${params.projName}.tsv"
    path "${params.projName}.txt"
    path "*.hmm.*", optional:true

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

    cat ${params.projName}.log >> Annotation.log
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
    path "Annotation.log"
    path "${params.projName}.gff", emit: gff
    path "${params.projName}.faa"
    path "${params.projName}.fna"
    path "${params.projName}.gbk"

    shell:
    //happens in work directory
    //RATT version needs to be custom EDGE version
    '''
    mkdir -p ./RATT/source
    cp !{gbk} ./RATT/source/source.gbk
    genbank2embl.pl ./RATT/source/source.gbk
    cd RATT
    runRATT.sh $PWD/source ../!{contigs} !{params.projName} Species 1>>../Annotation.log 2>&1
    cd ..
    cat ./RATT/*final.embl | fix_RATT_embl_feature.pl - > RATT/all.embl && embl2genbank.pl RATT/all.embl !{params.projName}.gbk 
    genbank2fasta.pl -translation !{params.projName}.gbk > !{params.projName}.faa
    genbank2fasta.pl -genome !{params.projName}.gbk > !{params.projName}.fna
    genbank2gff3.pl -e 3 --outdir stdout --DEBUG --typesource contig $PWD/!{params.projName}.gbk >!{params.projName}.gff
    '''
}


process annPlot {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )
    
    input:
    path gff

    output:
    path "plot_gff3.log"
    path "annotation_stats_plots.pdf"
    path "${params.projName}.txt", optional:true

    script:
    def rattReport = params.annotateProgram.equalsIgnoreCase("ratt") ? "awk \'\$1 ~ /CDS|RNA/ {print \$1\": \"\$2}' plot_gff3.log > ${params.projName}.txt" : ""
    """
    plot_gff3_stats.pl --input $gff --title $params.projName --prefix ./annotation_stats --outfmt PDF 1>plot_gff3.log 2>&1
    $rattReport
    """
}

process keggPlots {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Annotation",
        mode: 'copy'
    )

    input:
    path gff
    
    output:
    path "kegg_map/*"
    path "kegg_map.log"
    
    script:
    //TODO: will eventually need UI integration.
    """
    check_server_up.pl --url "http://rest.kegg.jp" && \
    opaver_anno.pl -g $gff -o ./kegg_map -p $params.projName > kegg_map.log 2>&1
    """
    }


workflow {

    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text

    contig_ch = channel.fromPath(params.inputContigs, checkIfExists:true)
    kingdom_ch = channel.of(params.taxKingdom)
    hmm_ch = channel.fromPath(params.customHMM, checkIfExists:true)
    prot_ch = channel.fromPath(params.customProtein, checkIfExists:true)
    if (params.annotateProgram =~ /(?i)prokka/) {
        if (params.taxKingdom == null) {
        kingdom_ch = unsetKingdom(contig_ch)
        }
        prokkaAnnotate(contig_ch, prot_ch, hmm_ch, kingdom_ch)
        annPlot(prokkaAnnotate.out.gff)
        if(params.keggView == true) {
            keggPlots(prokkaAnnotate.out.gff)
        }
    }
    else if (params.annotateProgram =~ /(?i)ratt/) {
        gb_ch = channel.fromPath(params.sourceGBK, checkIfExists:true)
        rattAnnotate(contig_ch, gb_ch)
        annPlot(rattAnnotate.out.gff)
        if(params.keggView == true) {
            keggPlots(rattAnnotate.out.gff)
        }
    }


}