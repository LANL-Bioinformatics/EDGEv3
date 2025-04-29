process countFastq {
    label "countFastq"
    label "tiny"

    input:
    path paired
    path unpaired

    output:
    path "fastqCount.txt", emit: counts
    path "all.{1,2}.fastq", emit: paired, optional:true
    path "all.se.fastq", emit: unpaired, optional:true

    script:

    if(paired.size() > 1 && paired[0] =~ /NO_FILE/) {
        paired = paired.tail().join(" ")
    }
    else {
        paired = paired.join(" ")
    }
    if(unpaired.size() > 1 && unpaired[0] =~ /NO_FILE/) {
        unpaired = unpaired.tail().join(" ")
    }
    else {
        unpaired = unpaired.join(" ")
    }
    

    paired_list = paired.startsWith("NO_FILE") ? "" : "-p ${paired}"
    unpaired_list = unpaired.startsWith("NO_FILE2") ? "" : "-u ${unpaired}"

    """
    getAvgLen.pl\
    $paired_list\
    $unpaired_list\
    -d .
    """
}

//gets average read length from fastqCount.txt
process avgLen {
    label "countFastq"
    label "tiny"

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
    paired
    unpaired

    main:

    countFastq(paired, unpaired)
    avgReadLen = avgLen(countFastq.out.counts)
    paired = countFastq.out.paired
    unpaired = countFastq.out.unpaired
    counts = countFastq.out.counts

    emit:
    avgReadLen
    counts
    paired
    unpaired
}