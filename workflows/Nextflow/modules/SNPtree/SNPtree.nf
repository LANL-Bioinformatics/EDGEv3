#!/usr/bin/env nextflow

process prepareSNPphylogeny {
    label 'phyl'
    containerOptions "--bind=${settings["snpDBbase"]}:/venv/bin/database"

    publishDir(
        path: "${settings["phylogenyOutDir"]}",
        mode: 'copy',
        saveAs: {
            filename ->
            if(filename.startsWith("results/tables/${settings["projName"]}")) {
                filename.drop(15)
            }
            else if(filename.startsWith("results/alignments/${settings["projName"]}")) {
                filename.drop(19)
            }
            else if(filename.endsWith(".fasttree")) {
                null
            }
            else {
                null
            }
        }
    )

    input:
    val settings
    path paired
    path unpaired
    path contigs
    output:
    path "results/trees/{*all.fasttree,RAxML_bestTree.*_all}", emit: allTree
    path "results/trees/{*cds.fasttree,RAxML_bestTree.*_cds}", emit:cdsTree
    path "results/alignments/${settings["projName"]}_*", optional: true
    path "results/tables/${settings["projName"]}_*", optional:true
    path "results/RaXML_*", optional: true
    path "results/trees/RAxML_bipartitions.*_all_best", emit: bootstrapTree, optional: true
    path "annotation.txt", emit: phyloAnn

    script:

    def kingdom = settings["taxKingdom"] != null ? "-kingdom ${settings["taxKingdom"]}" : ""
    def bootstrap = settings["phameBootstrap"] != false ? "-bootstrap -bootstrap_n ${settings["phameBootstrapNum"]}" : ""
    def db = settings["snpDBname"] != null ? "-db ${settings["snpDBname"]}" : ""
    def genomeNames = settings["snpGenomes"].size() != 0 ? "-genomesList ${settings["snpGenomes"].join(",")}" : ""
    def genomeFiles = settings["snpGenomesFiles"].size() != 0 ? "-genomesFiles ${settings["snpGenomesFiles"].join(",")}" : ""
    def reference = settings["snpRefGenome"] != null ? "-reference ${settings["snpRefGenome"]}" : ""
    def pair = paired.name != "NO_FILE" ? "-p $paired" : ""
    def single = unpaired.name != "NO_FILE2" ? "-s $unpaired" : ""
    def contig = contigs.name != "NO_FILE3" ? "-c $contigs" : ""
    def lr = (settings["fastqSource"] != null && (settings["fastqSource"].equalsIgnoreCase("nanopore") || settings["fastqSource"].equalsIgnoreCase("pacbio"))) ? "-nanopore" : ""

    """
    prepare_SNP_phylogeny.pl \
    -o . \
    -n ${settings["projName"]} \
    -tree ${settings["treeMaker"]} \
    -cpu ${settings["cpus"]} \
    -kingdom ${settings["taxKingdom"]} \
    -map "/venv/bin/database/SNPdb/reference.txt" \
    -bwa_id_map "/venv/bin/database/bwa_index/id_mapping.txt"\
    -bwa_genome_index "/venv/bin/database/bwa_index/NCBI-Bacteria-Virus.fna" \
    $bootstrap \
    $db \
    -db_path "/venv/bin/database" \
    $genomeNames \
    $genomeFiles \
    $reference \
    $pair \
    $contig \
    $single \
    $lr >> log.txt

    phame phame.ctrl 1>> log.txt 2>&1
    """

}

process prepareXMLphylogeny {
    label 'phyl'
    
    publishDir(
        path: "${settings["phylogenyOutDir"]}",
        mode: 'copy',
        saveAs: {
            filename ->
            if(filename.endsWith("_cds.fasttree.nwk") || filename.endsWith("cds.nwk")) {
                "SNPphyloTree.cds.nwk"
            }
            else if(filename.endsWith("_cds.fasttree.xml") || filename.endsWith("cds.xml")) {
                "SNPphyloTree.cds.xml"
            }
            else if(filename.endsWith("fasttree.nwk")|| filename.endsWith("all.nwk")) {
                "SNPphyloTree.all.nwk"
            }
            else if(filename.endsWith("fasttree.xml")|| filename.endsWith("all.xml")) {
                "SNPphyloTree.all.xml"
            }
        }
    )

    input:
    val settings
    each path(tree)
    path ann

    output:
    path "*.nwk"
    path "*.xml"

    script:
    """
    newickToPhyloXML.pl -m -i  \
    $tree -o . \
    $ann
    """
}




workflow PHYLOGENETICANALYSIS {

    take:
    settings
    paired
    unpaired
    contigs

    main:
    prepareSNPphylogeny(settings,paired,unpaired,contigs)
    resultTree = prepareSNPphylogeny.out.bootstrapTree.concat(prepareSNPphylogeny.out.allTree).first()
    if(settings["snpRefGenome"] != null || (settings["snpDBname"] != null && !(settings["snpDBname"]).equalsIgnoreCase("hantavirus"))) {
        resultTree = resultTree.concat(prepareSNPphylogeny.out.cdsTree)
    }
    prepareXMLphylogeny(settings,
                        resultTree,
                        prepareSNPphylogeny.out.phyloAnn)


}