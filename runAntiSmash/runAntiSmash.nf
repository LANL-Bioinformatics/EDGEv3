process antiSmash {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/AntiSmash",
        mode: 'copy'
    )
    input:
    path input
    output:
    path "*"
    script:
    def clusterblast = params.clusterblast == true ? "--cb-general" : ""
    def subclusterblast = params.subclusterblast == true ? "--cb-subclusters" : ""
    def knownclusterblast = params.knownclusterblast == true ? "--cb-knownclusters" : ""
    def mibig = params.mibig == true ? "--cc-mibig" : ""
    def smcogs = params.smcogs == true ? "--smcog-trees" : ""
    def asf = params.asf == true ? "--asf" : ""
    def rre = params.rre == true ? "--rre" : ""
    def fullhmmer = params.fullhmm == true ? "--fullhmmer" : ""
    def tigrfam = params.tigrfam == true ? "--tigrfam" : ""
    def pfam2go = params.pfam2go == true ? "--pfam2go" : ""
    def cassis = (params.cassis == true && params.taxon.equalsIgnoreCase("fungi")) ? "--cassis" : ""
    //antismash no longer supports fungal genefinding with prodigal-m (and will drop glimmerHMM in future releases)
    def genefinding = params.taxon.equalsIgnoreCase("fungi") ? "--genefinding-tool none" : "--genefinding-tool prodigal-m"

    """
    antismash -c $params.numCPU --taxon $params.taxon \
    --logfile antismashLog.txt --output-dir ./output \
    --html-title $params.projName --database $params.database \
    $clusterblast \
    $subclusterblast \
    $knownclusterblast \
    $mibig \
    $smcogs \
    $asf \
    $rre \
    $fullhmmer \
    $tigrfam \
    $pfam2go \
    $cassis \
    $genefinding \
    $input
    """

}


workflow {

    antiSmash(channel.fromPath(params.input,checkIfExists:true))

}