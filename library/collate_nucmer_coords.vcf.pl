#!/usr/bin/perl
#Eric Morrison
#9/18/2024
#Usage: collate_nucmer_coords.pl [coords files for pairwise collation]
# This script parses a set of nucmer coords and outputs the shared coordinates into a tab delimited file with one entry per shared position (contig \t pos). This file can then be used to parse a VCF file, e.g., bcftools view -R regions.tsv my.vcf

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
    
    foreach my $tig (sort{$a cmp $b} keys %coord1){
        if(defined($coord2{$tig}) == 0){next;} # can already move to next tig if current tig doesn't exist in second coordinate hash
        foreach my $start1 (sort{$a <=> $b} keys %{ $coord1{$tig} } ){
            
            foreach my $start2 (sort{$a <=> $b} keys %{ $coord2{$tig} }){
                if($start2 > $coord1{$tig}{$start1} || $coord2{$tig}{$start2} < $start1){next;} #check for complete non overlap
                #the above was a huge speed up. nice.
                my @seq1 = ($start1..$coord1{$tig}{$start1});
                my @seq2 = ($start2..$coord2{$tig}{$start2});
                    
                my($start,$stop) = compare_coords(\@seq1, \@seq2);
                if($start ne 0){
                    $collated{$tig}{$start} = $stop;
                    delete($coord2{$tig}{$start2});
                    last;
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
    my $collated_ref = $_[0];
    my %collated = %$collated_ref;
    open(COL, ">vcf.collated.coords") || die "Can't open collated coordinates output.\n";
    
    foreach my $tig (sort{$a cmp $b} keys %collated){
        foreach my $start (sort{$a <=> $b} keys %{ $collated{$tig}} ){
            for my $i ($start .. $collated{$tig}{$start}){
                print COL $tig, "\t", $i, "\n";
            }
        }
    }
}

#MAIN
{
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
    process_seq_print($collated_ref);
}
