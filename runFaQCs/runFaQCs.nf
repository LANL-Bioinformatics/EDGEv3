#!/usr/bin/env nextflow

//plotting for trimmed reads from ONT
process nanoplot {
    publishDir(
        path: "$params.outDir/QcReads",
        mode: 'copy'
    )
    input:
    path unpaired

    output:
    path "*" //lots of output plots

    script:
    """
    NanoPlot --fastq $unpaired --N50 --loglength -t $params.numCPU -f pdf --outdir . 2>/dev/null
    """

}


//Porechop for removing adapters from ONT or PacBio reads
process porechop {
    publishDir(
        path: "$params.outDir/QcReads",
        mode: 'copy'
    )


    input:
    path trimmed
    path log
    output:
    path "*.porechop.fastq", emit: porechopped
    
    script:
    """
    porechop -i $trimmed -o ./QC.unpaired.porechop.fastq -t $params.numCPU > $log
    """
}

//double-checks that any provided adapter file is in FASTA format
process adapterFileCheck {
    input:
    path adapterFile

    output:
    stdout

    script:
    """
    isFasta.pl $adapterFile
    """
}

//creates fastqCount.txt as preparation for downstream processes
process lenFile {
    input:
    path paired
    path unpaired

    output:
    path "fastqCount.txt"

    script:
    def paired_list = paired.name != "NO_FILE" ? "-p ${paired}" : ""
    def unpaired_list = unpaired.name != "NO_FILE2" ? "-u ${unpaired}" : ""

    """
    getAvgLen.pl\
    $paired_list\
    $unpaired_list\
    -d .
    """
}

//gets average read length from fastqCount.txt
process avgLen {
    input:

    path countFastq

    output:
    stdout

    shell:
    '''
    #!/usr/bin/env perl
    my $fastq_count_file = "./!{countFastq}";
    my $total_count = 0;
    my $total_len = 0;
    open (my $fh, "<", $fastq_count_file) or die "Cannot open $fastq_count_file\n";
    while(<$fh>){
        chomp;
        my ($name,$count,$len,$avg) = split /\t/,$_;
        $total_count += $count;
        $total_len += $len;
    }
    close $fh;
    my $avg_len = ($total_count > 0)? $total_len/$total_count : 0;
    print "$avg_len";
    '''
}

//main QC process. puts parameters together and runs FaQCs.
//EDGE currently uses a custom script (illumina_fastq_QC.pl) to handle QC for long reads,
//but it was unable to create report files when I attempted using it. For now, all input reads go through FaQCs.
process qc {
    publishDir(
        path: "$params.outDir/QcReads",
        mode: 'copy'
    )

    input:
    path paired
    path unpaired
    path adapter
    val avgLen

    output:
    path "QC.{1,2}.trimmed.fastq", optional:true
    path "QC.unpaired.trimmed.fastq", optional:true, emit: unpairedTrimmed
    path "QC_qc_report.pdf", optional: true
    path "QC.stats.txt", optional: true
    path "QC.log", emit: log

    script:
    //adjust minLength
    def min = params.minLength
    if(params.minLength < 1) {
        min = Math.abs(params.minLength * avgLen.toInteger())
    }

    def qcSoftware = "FaQCs"
    // if(params.ontFlag || params.pacbioFlag) {
    //     qcSoftware = "illumina_fastq_QC.pl"
    // }
    def pairedArg = paired.name != "NO_FILE" ? "-1 ${paired[0]} -2 ${paired[1]}" : ""
    // if(pairedArg != "" && (params.ontFlag || params.pacbioFlag)) {
    //     pairedArg = "-p $paired"
    // }
    def unpairedArg = unpaired.name != "NO_FILE2" ? "-u $unpaired" : ""

    def adapter = adapter.name != "NO_FILE3" ? "--adapter --artifactFile $adapter" : ""
    def polyA = params.polyA ? "--polyA" : ""
    def trim = ""
    // if(params.ontFlag || params.pacbioFlag) {
    //     trim = "--trim_only"
    // }
    def ascii = params.phredOffset != null ? "--ascii $params.phredOffset" : ""

    """
    $qcSoftware $pairedArg $unpairedArg \
    -q $params.qualityCutoff --min_L $min --avg_q $params.avgQuality \
    -n $params.numN --lc $params.lowComplexity --5end $params.cut5end --3end $params.cut3end \
    --split_size $params.splitSize -d . -t $params.numCPU \
    $polyA \
    $trim \
    $ascii \
    1>QC.log 2>&1
    """
}

workflow {
    //setup for optional files
    "mkdir nf_assets".execute().text
    "touch nf_assets/NO_FILE".execute().text
    "touch nf_assets/NO_FILE2".execute().text
    "touch nf_assets/NO_FILE3".execute().text 

    //input setup
    paired_ch = channel.fromPath(params.pairFile, checkIfExists:true).collect()
    unpaired_ch = channel.fromPath(params.unpairFile, checkIfExists:true).collect()

    //adapter setup
    adapter_ch = channel.empty()
    if(params.adapter != null) {
        adapter_ch = channel.fromPath(params.adapter, checkIfExists:true)
        isFasta(adapter_ch)
        if(isFasta.out != "Yes") {
            adapter_ch = channel.empty()
        }
    }

    //average read length calculation and main QC process
    avg_len_ch = avgLen(lenFile(paired_ch, unpaired_ch))
    qc(paired_ch, unpaired_ch, adapter_ch.ifEmpty("${workflow.projectDir}/nf_assets/NOFILE3"), avg_len_ch)
    
    //long read trimming and plotting
    if(params.ontFlag) {
        nanoplot_ch = channel.empty()
        if(params.porechop) {
            porechop(qc.out.unpairedTrimmed, qc.out.log)
            nanoplot(porechop.out.porechopped)
        }
        else {
            nanoplot(qc.out.unpairedTrimmed)
        }
    }
}