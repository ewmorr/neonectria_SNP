#!/usr/bin/perl
# Eric Morrison
# 9/24/2024
# Usage: cat_fasta.pl [fasta] [header text]
#This script concatenates all of the sequences in a fasta file to a single sequence. The fasta header is derived from the header text (one word).

use strict;
use warnings;

sub process_fasta{
    my $fas = $_[0];
    open(FAS, "$fas") || die "Can't open fasta\n";
    chomp(my @fas = <FAS>);
    $fas = join(":::::", @fas);
    @fas = split(">", $fas);
    shift @fas;
    my @fas_join;
    foreach my $line (@fas){
        my @line = split(":::::", $line);
        shift(@line); #rm header
        $line = join("", @line);
        push(@fas_join, $line);
    }
    my $seq = join("", @fas_join);
    return($seq);
}

#MAIN
{
    my($fas, $head) = @ARGV;
    my $seq = process_fasta($fas);
    my @seq = ( $seq =~ m/.{1,80}/g );
    print ">$head\n";
    foreach my $line (@seq){
        print "$line\n";
    }
}
