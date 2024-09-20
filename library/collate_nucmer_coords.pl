#!/usr/bin/perl
#Eric Morrison
#9/18/2024
#Usage: collate_nucmer_coords.pl [ref_genome.fasta] [coords files for pairwise collation]
# This script requires a reference (genome) sequence and at least two files of pairwise alignment coordinates from nucmer whole genome alignment. It returns the tab delimited set of coordinates that is inclusive of all pairwise alignments with respect to the reference (tig start stop length). A fasta file (core.ref_genome.fasta) with one entry per alignment named by ">reference name_start-stop" is also output. The same genome should be used as the ref in all nucemer alignments.

use strict;
use warnings;

sub hash_coords{
    my $coord = $_[0];
    open(COORD, "$coord") || die "Can't open coordinate file $coord\n";
    chomp(my @coord = <COORD>);
    
    my %coord;
    foreach my $line (@coord){
        if($line !~ /^\d+/){next;}
        my @line = split("\t", $line);
        #hash indexed by contig then start site with end site as value
        $coord{$line[11]}{$line[0]} = $line[1];
    }
    return(\%coord);
}

sub loop_coords{
    my ($ref1, $ref2) = @_;
    my %coord1 = %$ref1;
    my %coord2 = %$ref2;
    
    my %collated;
    
    foreach my $tig1 (sort{$a cmp $b} keys %coord1){
        foreach my $start1 (sort{$a <=> $b} keys %{ $coord1{$tig1} } ){
            
            foreach my $tig2 (sort{$a cmp $b} keys %coord2){
                if($tig1 ne $tig2){next;}
                foreach my $start2 (sort{$a <=> $b} keys %{ $coord2{$tig2} }){
                    if($start2 > $coord1{$tig1}{$start1} || $coord2{$tig2}{$start2} < $start1){next;} #check for complete non overlap
                    #the above was a huge speed up. nice.
                    my @seq1 = ($start1..$coord1{$tig1}{$start1});
                    my @seq2 = ($start2..$coord2{$tig2}{$start2});
                    
                    my($start,$stop) = compare_coords(\@seq1, \@seq2);
                    if($start ne 0){
                        $collated{$tig1}{$start} = $stop;
                        delete($coord2{$tig2}{$start2});
                        last;
                    }
                }
            }
        }
    }
    return(\%collated);
}
                    
sub compare_coords{
    my($seq1, $seq2) = @_;
    my @seq1 = @$seq1;
    my @seq2 = @$seq2;
    
    my $start = 0;
    my $stop = 0;
    foreach my $val1 (@seq1){
        foreach my $val2 (@seq2){
            #first test if the seq2 set is inclusive of seq1 val, if not move to next seq1 val
            if($val2 > $val1){last;}
            if($val1 == $val2){
                $start = $val1;
                last;
            }
        }
        if($start != 0){
            last;
        }
    }
    #loop in reverse for stop site search
    foreach my $val1 (reverse(@seq1)){
        foreach my $val2 (reverse(@seq2)){
            #first test if the seq2 set is inclusive of seq1 val, if not move to next seq1 val
            if($val2 < $val1){last;}
            if($val1 == $val2){
                $stop = $val1;
                last;
            }
        }
        if($stop != 0){
            last;
        }
    }
    return($start, $stop)
}

sub process_seq_print{
    my($refseq, $collated_ref) = @_;
    
    open(SEQ, "$refseq") || die "Can't open $refseq\n";
    chomp(my @refseq = <SEQ>);
    my $seqs = join("::::::::::", @refseq);
    my @seqs = split(">", $seqs);
    shift(@seqs);
    my %seq;
    foreach my $seq (@seqs){
        my @seq = split("::::::::::", $seq);
        my $head = shift(@seq);
        $seq = join("", @seq);
        $seq{$head} = $seq;
    }
    close(SEQ);
    
    open(REF, ">core.$refseq") || die "Can't open fasta ouput.\n";
    open(COL, ">collated.coords") || die "Can't open collated coordinates output.\n";
    
    my %collated = %$collated_ref;
    foreach my $tig (sort{$a cmp $b} keys %collated){
        foreach my $start (sort{$a <=> $b} keys %{ $collated{$tig}} ){
            my $offset = $start-1; #do -1 for both zero indexing and length calc (which should be *inclusive* of the start)
            my $length = $collated{$tig}{$start} - $offset;
            my $subseq = substr($seq{$tig}, $offset, $length);
            print REF ">$tig"."_$start"."-$collated{$tig}{$start}\n$subseq\n";
            print COL "$tig\t$start\t$collated{$tig}{$start}\t$length\n";
        }
    }
}

#MAIN
{
    my $refseq = shift(@ARGV);
    my @coords = @ARGV;
    
    my @coords_refs;
    for(my $i = 0; $i<@coords; $i++){
        $coords_refs[$i] = hash_coords($coords[$i]);
    }
    my $collated_ref;
    for(my $i = 1; $i<@coords_refs; $i++){
        if($i == 1){
            $collated_ref = loop_coords($coords_refs[$i-1], $coords_refs[$i]);
        } else {
            $collated_ref = loop_coords($collated_ref, $coords_refs[$i]);
        }
    }
    process_seq_print($refseq, $collated_ref);
}
