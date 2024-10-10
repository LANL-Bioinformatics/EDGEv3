process antiSmashFungi {
    publishDir(
        path: "$params.outDir/AssemblyBasedAnalysis/AntiSmash",
        mode: 'copy'
    )

    input:
    path input
    path gff

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
    def cassis = params.cassis == true ? "--cassis" : ""
    """
    antismash -c $params.numCPU --taxon fungi \
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
    --genefinding-tool glimmerhmm --genefinding-gff3 $gff\
    $input
    """
}
process antiSmashBacteria {
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

    """
    antismash -c $params.numCPU --taxon bacteria \
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
    --genefinding-tool prodigal-m \
    $input
    """

}


workflow {

    if (params.taxon.equalsIgnoreCase("bacteria")) {
        antiSmashBacteria(channel.fromPath(params.input,checkIfExists:true))
    }
    else if (params.taxon.equalsIgnoreCase("fungi")) {
        antiSmashFungi(channel.fromPath(params.input,checkIfExists:true), channel.fromPath(params.gffInput, checkIfExists:true))
    }
}