#!/usr/bin/perl
#Eric Morrison
#05172023
#Usage: snp_table2fasta.full_ref_seq.sep_contigs.pl [snp table as produced by gatk VariantsToTable] [ref_seq.fasta] [number of cols of VCF info]
#Any NAs and non-standard nucelotides in the SNP table are converted to N
#Any indels should be removed from the SNP table `grep -v "INDEL" snp.table > snp.noIndel.table`
# The first row is assumed to be header info with subsequent rows being SNPs
# Each contig is written to a separate fasta with the fastas named by contig name

use strict;
use warnings;

sub hash_refSeq{
    my $refSeq = $_[0];
    open(REF, "$refSeq") || die "Can't open reference fasta\n";
    
    my %refSeq;
    my $id;
    while(my $line = <REF>){
        chomp($line);
        if( $line =~ /^>/){
            $id = $line;
            $id =~ s/^>//;
            $refSeq{$id} = "";
        }else{
            $refSeq{$id} .= $line;
        }
    }
    return(\%refSeq);
}

sub hash_gts{
    my ($snps, $gtIndex) = @_;
    open(IN, "$snps") || die "Can't open snp table\n";
    chomp(my @in = <IN>);
    
    #pull the sample IDs and remove the trailing ".GT"
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
            #SNPs indexed by contig name, sample, position
            $snps{$snp[0]}{$gt_labels[$i]}{$snp[1]} = $snp[$i];
        }
    }
    return(\%snps);
}

#MAIN
{
    my $snps = $ARGV[0];
    my $refSeq = $ARGV[1];
    my $gtIndex = $ARGV[2];
    
    my $refSeq_ref = hash_refSeq($refSeq);
    my $gt_ref = hash_gts($snps, $gtIndex); #we still hash the gts so we can modify the seqs one at a time and print on the fly instead of holding in memory
 
    my %gts = %$gt_ref;
    
    #printing loop
    #open new file for each contig
    #modify sequence for each sample and then print
    foreach my $contig (sort {$a cmp $b} keys %gts){
        open(TIG, ">$contig.fasta") || die "Can't open ouput file\n";
        
        foreach my $gt (sort {$a cmp $b} keys %{ $gts{$contig} }){
            my @sequence = split("", $$refSeq_ref{$contig}); #reset to refSeq for each new genotype
            
            foreach my $pos (sort {$a <=> $b} keys %{ $gts{$contig}{$gt} }){
                $sequence[$pos - 1] = $gts{$contig}{$gt}{$pos}
            }
            print TIG ">", $gt, "\n", join("", @sequence), "\n";
        }
        close(TIG);
    }
}
