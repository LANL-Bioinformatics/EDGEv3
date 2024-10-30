//André Watson
//Aug 2024
//apwat @ lanl.gov

//base process. Takes a FASTA file containing contigs and performs taxonomic analysis with MICCR (https://github.com/chienchi/miccr).
process contigTaxonomy {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Taxonomy",
        mode: 'copy'
    )
    input:
    path contigs

    output:
    path "${params.projName}.log"
    path "log.txt"
    path "${params.projName}.lca_ctg.tsv", emit: taxLcaResult
    path "${params.projName}.ctg.tsv", emit: taxResult
    path "${params.projName}.unclassified.fasta"
    path "${params.projName}.paf"
    
    script:
    """
    miccr.py -x asm10 -d $params.dbPath -t $params.cpus -p $params.projName -i $contigs 1>log.txt 2>&1 
    get_unclassified_fasta.pl -in $contigs -classified ${params.projName}.lca_ctg.tsv -output ${params.projName}.unclassified.fasta -log log.txt
    """
}

//adds multi-level taxonomic classification to results file. Takes in a .ctg.tsv file produced by MICCR.
process addLineage {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Taxonomy",
        mode: 'copy'
    )
    
    input:
    path taxResult

    output:
    path "*.lineage", emit: lineage

    script:
    def dbFolder = params.dbPath.take(params.dbPath.lastIndexOf("/")) //get folder path containing DB
    """
    add_lineage.py $dbFolder $taxResult > ${taxResult.name}.lineage
    """
}

//creates taxonomy classification graphs. Takes lineage file, .lca_ctg.tsv file produced by MICCR, 
//and a coverage table (see https://github.com/chienchi/miccr/blob/master/utils/README.md), or from workflow runReadsToContig 
process plotAndTable {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/Taxonomy",
	mode: 'copy'
    )
    input:
    path lineage
    path covTable
    path lcaResult
    
    output:
    path "${params.projName}.ctg_class.LCA.json"
    path "summary_by_*.txt"
    path "*.pdf"
    
    script:
    """
    classification_plot.R $lineage $params.projName $covTable
    tab2Json_for_dataTable.pl -project_dir $params.outDir -mode contig -limit $params.rowLimit $lcaResult > ${params.projName}.ctg_class.LCA.json
    """
}

workflow {
    contigs = channel.fromPath(params.contigFile, checkIfExists:true)
    coverageTable = channel.fromPath(params.coverageTable, checkIfExists:true)
    contigTaxonomy(contigs)
    addLineage(contigTaxonomy.out.taxResult)
    plotAndTable(addLineage.out.lineage, coverageTable, contigTaxonomy.out.taxLcaResult)
}