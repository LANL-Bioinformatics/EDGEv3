process runBinning {
    label 'binning'
    publishDir(
        path: "${settings["readBinningOutDir"]}",
        mode: 'copy'
    )

    input:
    val settings
    path contigs
    path abundances

    output:
    path "*" //many outputs

    script:
    """
    run_MaxBin.pl \
    -contig $contigs \
    -out ${settings['projName']}_bin \
    -abund $abundances -thread ${settings['cpus']} \
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
    abundances


    main:
    
    runBinning(settings, contigs, abundances)

}