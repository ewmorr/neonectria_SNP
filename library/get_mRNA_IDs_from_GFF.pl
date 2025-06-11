#!/usr/bin/perl
#Eric Morrison
#11032021
#Usage: get_mRNA_IDs_from_GFF.pl [GFF] > [GFF pos and IDs filtered for mRNA ]

use strict;
use warnings;


sub clean_gff{ #pull only those lines with mRNA as feature, extract ID, write positions and ID
    my $in = $_[0];
    open(IN, "$in") || die "can't open gff\n";
    chomp(my @in = <IN>);
    
    foreach my $line (@in){
    #    if($line =~ /^#/){next;}
            
        my @line = split("\t", $line);
        if(scalar(@line) == 1){next;} #to deal with uninitialized value warnings from sequences
        if($line[2] ne "mRNA"){next;}
        $line[8] =~ /ID=(.*?)\;/;
        print $line[0], "\t", $1, "\t", $line[3], "\t", $line[4], "\n";
    }
}
#MAIN
{
    my $gff = $ARGV[0];
    
    print "contig\tgeneID\tstart\tstop\n";
    clean_gff($gff);
}
