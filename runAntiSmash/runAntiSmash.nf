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
    def clusterblast = params.antismash-clusterblast == true ? "--cb-general" : ""
    def subclusterblast = params.antismash-subclusterblast == true ? "--cb-subclusters" : ""
    def knownclusterblast = params.antismash-knownclusterblast == true ? "--cb-knownclusters" : ""
    def mibig = params.antismash-mibig == true ? "--cc-mibig" : ""
    def smcogs = params.antismash-smcogs == true ? "--smcogs-trees" : ""
    def asf = params.antismash-asf == true ? "--asf" : ""
    def rre = params.antismash-rre == true ? "--rre" : ""
    def fullhmmer = params.antismash-fullhmm == true ? "--fullhmmer" : ""
    def tigrfam = params.antismash-tigrfam == true ? "--tigrfam" : ""
    def pfam2go = params.antismash-pfam2go == true ? "--pfam2go" : ""
    def cassis = (params.antismash-cassis == true && params.antismash-taxon.equalsIgnoreCase("fungi")) ? "--cassis" : ""

    """
    antismash -c $params.numCPU --taxon ${params.antismash-taxon} \
    --logfile antismashLog.txt --output-dir $params.outDir \
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
    --genefinding-tool prodigal-m $input
    """

}


workflow {

    antiSmash(channel.fromPath(params.input,checkIfExists:true))

}