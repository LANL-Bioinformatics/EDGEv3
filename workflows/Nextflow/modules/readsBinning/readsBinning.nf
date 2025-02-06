process runBinning {

    publishDir(
        path: "${settings["outDir"]}/AssemblyBasedAnalysis/Binning",
        mode: 'copy'
    )

    input:
    val settings
    path contigs

    output:

    script:
    abundanceFile = channel.fromPath(settings["abundFile"], checkIfExists:true)

    """
    run_MaxBin.pl \
    -contig $contigs \
    -out ${project_name}_bin \
    -abund $abundanceFile -thread ${settings['cpus']} \
    -plotmarker -min_contig_length ${settings["binningMinLength"]} \
    -max_iteration ${settings["binningMaxItr"]} \
    -prob_threshold ${settings["binningProb"]} \
    -markerset ${settings["binningMarkerSet"]}
    """
}


workflow BINNING {
    take:
    settings
    contigs


    main:
    
    runBinning(settings, contigs)

}