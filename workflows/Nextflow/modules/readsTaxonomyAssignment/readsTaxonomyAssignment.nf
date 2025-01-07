#!/usr/bin/env nextflow

//main RTA process
process readsTaxonomy {
    label 'rta'
    publishDir(
        path: "${settings["outDir"]}/ReadsBasedAnalysis/Taxonomy",
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
    def debugging = (settings["debugFlag"] != null && settings["debugFlag"]) ? "--debug " : ""
    def numCPU = settings["cpus"] != null ? settings["cpus"] : 8
    """
    cat $paired $unpaired > allReads.fastq
    microbial_profiling.pl $debugging -o . \
    -s $taxonomyConfig \
    -c $numCPU \
    $debugging \
    allReads.fastq 2>>$errorlog

    svg2pdf.sh  report/*/*/*.svg 2>>$errorlog
    """
}

//creates RTA config file based on input settings
process readsTaxonomyConfig {
    label 'rta'

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
    def tools = settings["enabledTools"] != null ? "-tools \'${settings["enabledTools"]}\' " : ""
    def template = settings["template"] != null ? "-template ${settings["template"]} " : ""
    def splitTrimMinQ = settings["splitTrimMinQ"] != null ? "-splitrim-minq ${settings["splitTrimMinQ"]} " : ""

    def bwa = settings["custom_bwa_db"] != null ? "-bwa-db ${settings["custom_bwa_db"]} " : ""
    def metaphlan = settings["custom_metaphlan_db"] != null ? "-metaphlan-db ${settings["custom_metaphlan_db"]} " : ""
    def kraken = settings["custom_kraken_db"] != null ? "-kraken-db ${settings["custom_kraken_db"]} " : ""
    def centrifuge = settings["custom_centrifuge_db"] != null ? "-centrifuge-db ${settings["custom_centrifuge_db"]} " : ""
    def pangia = settings["custom_pangia_db"] != null ? "-pangia-db ${settings["custom_pangia_db"]} " : ""
    def diamond = settings["custom_diamond_db"] != null ? "-diamond-db ${settings["custom_diamond_db"]} " : ""

    def gottcha_speDB_v = settings["custom_gottcha_speDB_v"] != null ? "-gottcha-v-speDB ${settings["custom_gottcha_speDB_v"]} " : ""
    def gottcha_speDB_b = settings["custom_gottcha_speDB_b"] != null ? "-gottcha-b-speDB ${settings["custom_gottcha_speDB_b"]} " : ""
    def gottcha_strDB_v = settings["custom_gottcha_strDB_v"] != null ? "-gottcha-v-strDB ${settings["custom_gottcha_strDB_v"]} " : ""
    def gottcha_strDB_b = settings["custom_gottcha_strDB_b"] != null ? "-gottcha-b-strDB ${settings["custom_gottcha_strDB_b"]} " : ""
    def gottcha_genDB_v = settings["custom_gottcha_genDB_v"] != null ? "-gottcha-v-genDB ${settings["custom_gottcha_genDB_v"]} " : ""
    def gottcha_genDB_b = settings["custom_gottcha_genDB_b"] != null ? "-gottcha-b-genDB ${settings["custom_gottcha_genDB_b"]} " : ""

    def gottcha2_genDB_v = settings["custom_gottcha2_genDB_v"] != null ? "-gottcha2-v-genDB ${settings["custom_gottcha2_genDB_v"]} " : ""
    def gottcha2_speDB_v = settings["custom_gottcha2_speDB_v"] != null ? "-gottcha2-v-speDB ${settings["custom_gottcha2_speDB_v"]} " : ""
    def gottcha2_speDB_b = settings["custom_gottcha2_speDB_b"] != null ? "-gottcha2-b-speDB ${settings["custom_gottcha2_speDB_b"]} " : ""

    def np = (settings["fastqSource"] != null && settings["fastqSource"].equalsIgnoreCase("nanopore")) ? "--nanopore " : ""

    """
    mkdir -p /venv/opt/krona/taxonomy
    touch /venv/opt/krona/taxonomy/taxdump.tar.gz
    chmod 777 /venv/opt/krona/taxonomy/taxdump.tar.gz
    updateTaxonomy.sh

    updateAccessions.sh 
    
    microbial_profiling_configure.pl $template \
    $tools -bwaScoreCut $bwaScoreCut\
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
