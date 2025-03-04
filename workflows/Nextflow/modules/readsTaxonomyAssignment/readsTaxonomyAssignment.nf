#!/usr/bin/env nextflow

//main RTA process
process readsTaxonomy {
    label 'rta'

    containerOptions "--compat --cleanenv \
                        --bind=${settings["baseDB"]}:/venv/bin/../../../database \
                        --bind=${settings["baseDB"]}:/venv/opt/krona/taxonomy"

    publishDir(
        path: "${settings["readsTaxonomyOutDir"]}",
	mode: 'copy',
	pattern: "{*.log,log**,report**}"
    )
    
    input:
    val settings
    path paired
    path unpaired
    path taxonomyConfig
    path errorlog

    output:
    //outputs are many and variable
    path "error.log", emit: logfile
    path "taxonomyProfiling.log"
    path "report/**/**/*.svg", emit: svgs
    path "log**"
    path "report**"

    script:
    def numCPU = settings["cpus"] != null ? settings["cpus"] : 8
    """
    cat $paired $unpaired > allReads.fastq
    microbial_profiling.pl -o . \
    -s $taxonomyConfig \
    -c $numCPU \
    allReads.fastq 2>>$errorlog

    svg2pdf.sh  report/*/*/*.svg 2>>$errorlog
    """
}

//creates RTA config file based on input settings
process readsTaxonomyConfig {
    label 'rta'

    containerOptions "--compat --cleanenv \
                        --bind=${settings["baseDB"]}:/venv/bin/../../../database \
                        --bind=${settings["baseDB"]}:/venv/opt/krona/taxonomy"

    input:
    val settings
    val avgLen

    output:
    path "error.log", emit: errorlog
    path "microbial_profiling.settings.ini", emit: config

    script:
    def bwaScoreCut = 30
    if (settings["fastqSource"] != null && (settings["fastqSource"].equalsIgnoreCase("nanopore") || settings["fastqSource"].equalsIgnoreCase("pacbio"))) {
        if (settings["minLen"] > 1000) {
            bwaScoreCut = settings["minLen"]
        } 
        else {
            bwaScoreCut=1000
        }
    }
    else{
        bwaScoreCut = (avgLen as Integer)*0.8
    }
    bwaScoreCut = bwaScoreCut as Integer
    tools = settings["enabledTools"] != null ? "-tools \'${settings["enabledTools"]}\' " : ""
    splitTrimMinQ = settings["splitTrimMinQ"] != null ? "-splitrim-minq ${settings["splitTrimMinQ"]} " : ""

    base = settings["baseDB"] != null ? "-base-db ${settings["baseDB"]}" : ""

    bwa = settings["custom_bwa_db"] != null ? "-bwa-db ${settings["custom_bwa_db"]} " : ""
    metaphlan = settings["custom_metaphlan_db"] != null ? "-metaphlan-db ${settings["custom_metaphlan_db"]} " : ""
    kraken = settings["custom_kraken_db"] != null ? "-kraken-db ${settings["custom_kraken_db"]} " : ""
    centrifuge = settings["custom_centrifuge_db"] != null ? "-centrifuge-db ${settings["custom_centrifuge_db"]} " : ""
    pangia = settings["custom_pangia_db"] != null ? "-pangia-db ${settings["custom_pangia_db"]} " : ""
    diamond = settings["custom_diamond_db"] != null ? "-diamond-db ${settings["custom_diamond_db"]} " : ""

    gottcha_speDB_v = settings["custom_gottcha_speDB_v"] != null ? "-gottcha-v-speDB ${settings["custom_gottcha_speDB_v"]} " : ""
    gottcha_speDB_b = settings["custom_gottcha_speDB_b"] != null ? "-gottcha-b-speDB ${settings["custom_gottcha_speDB_b"]} " : ""
    gottcha_strDB_v = settings["custom_gottcha_strDB_v"] != null ? "-gottcha-v-strDB ${settings["custom_gottcha_strDB_v"]} " : ""
    gottcha_strDB_b = settings["custom_gottcha_strDB_b"] != null ? "-gottcha-b-strDB ${settings["custom_gottcha_strDB_b"]} " : ""
    gottcha_genDB_v = settings["custom_gottcha_genDB_v"] != null ? "-gottcha-v-genDB ${settings["custom_gottcha_genDB_v"]} " : ""
    gottcha_genDB_b = settings["custom_gottcha_genDB_b"] != null ? "-gottcha-b-genDB ${settings["custom_gottcha_genDB_b"]} " : ""

    gottcha2_genDB_v = settings["custom_gottcha2_genDB_v"] != null ? "-gottcha2-v-genDB ${settings["custom_gottcha2_genDB_v"]} " : ""
    gottcha2_speDB_v = settings["custom_gottcha2_speDB_v"] != null ? "-gottcha2-v-speDB ${settings["custom_gottcha2_speDB_v"]} " : ""
    gottcha2_speDB_b = settings["custom_gottcha2_speDB_b"] != null ? "-gottcha2-b-speDB ${settings["custom_gottcha2_speDB_b"]} " : ""

    np = (settings["fastqSource"] != null && settings["fastqSource"].equalsIgnoreCase("nanopore")) ? "--nanopore " : ""

    """

    mkdir -p /venv/opt/krona/taxonomy
    touch /venv/opt/krona/taxonomy/taxdump.tar.gz
    chmod 777 /venv/opt/krona/taxonomy/taxdump.tar.gz
    updateTaxonomy.sh

    updateAccessions.sh 
    
    microbial_profiling_configure.pl \
    $tools -bwaScoreCut $bwaScoreCut\
    $base\
    $bwa\
    $metaphlan\
    $kraken\
    $centrifuge\
    $pangia\
    $diamond\
    $gottcha_speDB_v\
    $gottcha_speDB_b\
    $gottcha_strDB_v\
    $gottcha_strDB_b\
    $gottcha_genDB_v\
    $gottcha_genDB_b\
    $gottcha2_genDB_v\
    $gottcha2_speDB_v\
    $gottcha2_speDB_b\
    $np >microbial_profiling.settings.ini 2>error.log
    """

}

//TODO: add workflow logic for retrieving unmapped reads
workflow READSTAXONOMYASSIGNMENT {
    take:
    settings
    paired
    unpaired
    avgLen

    main:
    readsTaxonomyConfig(settings, avgLen)
    readsTaxonomy(settings, paired, unpaired, readsTaxonomyConfig.out.config, readsTaxonomyConfig.out.errorlog)

}
