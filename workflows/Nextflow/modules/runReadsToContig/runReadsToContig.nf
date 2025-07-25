#!/usr/bin/env nextflow

process validationAlignment {
    label 'r2c'
    label 'small'
    publishDir(
        path: "${settings["assemblyOutDir"]}/readsMappingToContig",
        mode: 'copy'
    )
    input:
    val settings
    val platform
    path paired
    path unpaired
    path contigs

    output:
    path "*.sort.bam", emit: sortedBam
    path "*.bai"
    path "*.alnstats.txt", emit: alnStats
    path "*_coverage.table", emit: cov_table
    path "*_plots.pdf", emit: contigPlots
    path "magnitudes.txt", emit: magnitudes
    path "Final_contigs.fasta", emit: contig_file, optional:true //not present if using already-assembled contigs
    path "mapping.log", emit: logFile

    script:
    def outPrefix = "readsToContigs"
    def paired = paired.name != "NO_FILE" ? "-p \'${paired[0]} ${paired[1]}\' " : ""
    def unpaired = unpaired.name != "NO_FILE2" ? "-u $unpaired " : ""
    def cutoff = settings["assembledContigs"] != "${projectDir}/nf_assets/NO_FILE3" ? "-c 0 " : "-c 0.1 "
    def max_clip = settings["r2gMaxClip"] != null ? "-max_clip ${settings["r2gMaxClip"]} " : ""


    def ont_flag = (platform != null && platform.contains("NANOPORE")) ? "-x ont2d " : ""
    def pb_flag = (platform != null && platform.contains("PACBIO")) ? "-x pacbio " : ""
    
    def aligner_options = ""
    if(settings["r2c_aligner"] =~ "bowtie") {
        def bowtie_options = settings["r2c_aligner_options"].replaceAll("-p\\s*\\d+","")
        if(!(bowtie_options =~ /-k/)) {
            bowtie_options += " -k 10 "
        }
        aligner_options = "-aligner bowtie -bowtie_options \'$bowtie_options\'"
    }
    else if(settings["r2c_aligner"] =~ "bwa") {
        def bwa_options = settings["r2c_aligner_options"].replaceAll("-t\\s*\\d+","")
        if (ont_flag != "") {
            unpaired = unpaired.replaceAll("-u ","-long ")
            bwa_options += ont_flag
        }
        if (pb_flag != "") {
            unpaired = unpaired.replaceAll("-u ","-long ")
            bwa_options += pb_flag
        }
        aligner_options = "-aligner bwa -bwa_options \'$bwa_options\'"
    }
    else if (settings["r2c_aligner"] =~ "minimap") { 
        def minimap_options = settings["r2c_aligner_options"].replaceAll("-t\\s*\\d+","")
        if(ont_flag != "" || pb_flag != "") {
            unpaired = unpaired.replaceAll("-u ","-long ")
        }
        if(pb_flag != "") {
            minimap_options += " -x map-pb "
        }
        aligner_options = "-aligner minimap2 -minimap2_options \'$minimap_options\'"
    }


    """
    runReadsToContig.pl \
    $cutoff\
    -cpu ${task.cpus}\
    $paired\
    $unpaired\
    -d . -pre $outPrefix\
    -ref $contigs \
    $max_clip\
    $aligner_options &> mapping.log


    awk \'{print \$1\"\\t\"\$4}\' ${outPrefix}_coverage.table > magnitudes.txt
    """

}

process makeJSONcoverageTable {
    label 'r2c'
    label 'tiny'
    publishDir(
        path: "${settings["assemblyOutDir"]}/readsMappingToContig",
        mode: 'copy',
        pattern: "*_coverage.table.json"
    )
    publishDir(
        path: "${settings["assemblyOutDir"]}",
        mode: 'copy',
        pattern: "*stats.{pdf,txt}"
    )
    input:
    val settings
    path cov_table
    path contigFile

    output:
    path "contigs_stats.txt"
    path "contigs_stats.pdf", emit: contigStatsReport
    path "*_coverage.table.json"

    script:
    def rowLimit = settings["rowLimit"] != null ? "${settings["rowLimit"]} " : "3000"
    
    """
    tab2Json_for_dataTable.pl -project_dir . -mode contig -limit $rowLimit  \
    readsToContigs_coverage.table > readsToContigs_coverage.table.json

    contig_stats.pl -p $contigFile > contigs_stats.txt
    """
}

process extractUnmapped {
    label 'r2c'
    label 'small'
    publishDir(
        path:"${settings["assemblyOutDir"]}/readsMappingToContig/",
        mode: 'copy',
        overwrite: true
    )
    input:
    val settings
    path bamFile
    path logFile


    output:
    path "mapping.log"
    path "Unmapped*.fastq"

    script:
    """
    echo "Extract unmapped reads" >> $logFile
    bam_to_fastq.pl -unmapped -prefix Unmapped ${bamFile} >> $logFile
    """

}

workflow READSTOCONTIGS {
    take:
    settings
    platform
    paired
    unpaired
    contigs

    main:
    validationAlignment(settings, platform, paired, unpaired, contigs)
    alnStats = validationAlignment.out.alnStats
    makeJSONcoverageTable(settings, validationAlignment.out.cov_table, validationAlignment.out.contig_file)
    if(settings["extractUnmapped"]) {
        extractUnmapped(settings, validationAlignment.out.sortedBam, validationAlignment.out.logFile)
    }

    covTable = validationAlignment.out.cov_table
    magnitudes = validationAlignment.out.magnitudes
    contigPlots = validationAlignment.out.contigPlots
    contigStatsReport = makeJSONcoverageTable.out.contigStatsReport

    emit:
    covTable
    magnitudes
    contigPlots
    contigStatsReport
    alnStats


}
