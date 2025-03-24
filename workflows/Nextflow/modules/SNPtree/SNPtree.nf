#!/usr/bin/env nextflow

process prepareSNPphylogeny {

    publishDir(
        path: "${settings["phylogenyOutDir"]}",
        mode: 'copy'
    )

    input:
    val settings
    path paired
    path unpaired
    path contigs
    output:
    path "results/trees/{*all.fasttree,RAxML_bestTree.*_all}", emit: allTree
    path "results/${settings["projName"]}_*"
    path "results/RaXML_*", optional: true
    path "results/RAxML_bipartitions.*_all_best", emit: bootstrapTree, optional: true
    path "annotation.txt", emit: phyloAnn

    script:

    def kingdom = settings["taxKingdom"] != null ? "-kingdom ${settings["taxKingdom"]}" : ""
    def bootstrap = settings["phameBootstrap"] != false ? "-bootstrap -bootstrap_n ${settings["phameBootstrapNum"]}" : ""
    def db = settings["snpDBname"] != null ? "-db ${settings["snpDBname"]}" : ""
    def genomeNames = settings["snpGenomes"].size() != 0 ? "-genomesList ${settings["snpGenomes"]}" : ""
    def genomeFiles = settings["snpGenomesFiles"].size() != 0 ? "-genomesFiles ${settings["snpGenomesFiles"]}" : ""
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
    -map "${settings["snpDBbase"]}/SNPdb/reference.txt" \
    -bwa_id_map "${settings["snpDBbase"]}/bwa_index/id_mapping.txt"\
    -bwa_genome_index "${settings["snpDBbase"]}/bwa_index/NCBI-Bacteria-Virus.fna" \
    $bootstrap \
    $db \
    -db_path "${settings["snpDBbase"]}" \
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
    publishDir(
        path: "${settings["phylogenyOutDir"]}",
        mode: 'copy'
    )

    input:
    val settings
    path tree
    path ann

    output:
    path "*"

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
    resultTree = prepareSNPphylogeny.out.bootstrapTree.concat(prepareSNPphylogeny.out.allTree)
    prepareXMLphylogeny(settings,
                        resultTree.first(),
                        prepareSNPphylogeny.out.phyloAnn)


}