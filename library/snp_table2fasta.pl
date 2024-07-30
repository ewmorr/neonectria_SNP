#!/usr/bin/perl
#Eric Morrison
#042023
#Usage: snp_table2fasta.pl [snp table as produced by gatk VariantsToTable] [number of cols of VCF info] > [fasta]
# The first row is assumed to be header info with subsequent rows being SNPs

use strict;
use warnings;

sub hash_gts{
    my ($snps, $gtIndex) = @_;
    open(IN, "$snps") || die "Can't open snp table\n";
    chomp(my @in = <IN>);
    
    my $gt_labels = shift @in;
    $gt_labels =~ s/\.GT//g;
    
    my @gt_labels = split("\t", $gt_labels);
    
    my %snps;
    
    foreach my $snp (@in){
        my @snp = split("\t", $snp);

        for(my $i = $gtIndex; $i < @gt_labels; $i++){
            if($snp[$i] !~ /^A$|^G$|^T$|^C$/){
                $snp[$i] = "N";
            }
            if($snp[$i] =~ /\./){
                $snp[$i] = "N";
            }
            push @{ $snps{$gt_labels[$i]} }, $snp[$i];
        }
    }
    return(\%snps);
}

#MAIN
{
    my $snps = $ARGV[0];
    my $gtIndex = $ARGV[1];
    my $gt_ref = hash_gts($snps, $gtIndex);
 
    my %gts = %$gt_ref;
    foreach my $gt (sort {$a cmp $b} keys %gts){
        print ">$gt\n";
        my $seq = join("", @{ $gts{$gt} });
        print join("", @{ $gts{$gt} }), "\n";
    }
}
