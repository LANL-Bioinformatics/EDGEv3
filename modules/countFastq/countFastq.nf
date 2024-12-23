//creates fastqCount.txt as preparation for downstream processes
process countFastq {
    label "countFastq"

    input:
    val settings
    path fastq

    output:
    path "fastqCount.txt", emit: counts
    path "all.*.fastq", emit: allFiles

    script:

    file_list = ""
    if(settings["pairedFile"]) {
        file_list = "-p $fastq"
    }
    else {
        file_list = "-u $fastq"
    }

    """
    getAvgLen.pl\
    $file_list\
    -d .
    """
}

//gets average read length from fastqCount.txt
process avgLen {
    label "countFastq"

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

//calculates average read length and concatenates input files
workflow COUNTFASTQ {
    take:
    settings
    inputFastq

    main:

    countFastq(settings, inputFastq)
    avgReadLen = avgLen(countFastq.out.counts)
    fastqFiles = countFastq.out.allFiles

    emit:
    avgReadLen
    fastqFiles
}