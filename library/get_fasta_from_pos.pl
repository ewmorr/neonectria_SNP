#!/usr/bin/perl
# Eric Morrison
# 9/23/2024
# Usage: get_fasta_from_pos.pl [fasta] [positions] > out.fasta
# This script slices a fasta file based on a tab-delimited list of fasta sequence headers (e.g., contigs) and coordinates, such as might be derived from a VCF file. We use this script to parse a reference sequence to those positions that are present in a VCF after filtering for quality.

use strict;
use warnings;

sub hash_pos{
    my $pos = $_[0];
    open(POS, "$pos") || die "Can't open pos\n";
    chomp(my @pos = <POS>);
    my %ps;
    my @tigOrder;
    foreach my $ps (@pos){
        if($ps =~ /CHROM/){next;}
        my @ps = split("\t", $ps);
        if(defined($tigOrder[0]) == 0){
            $tigOrder[0] = $ps[0];
        } elsif($tigOrder[-1] ne $ps[0]){
            push(@tigOrder, $ps[0]);
        }
        
        if(defined($ps{$ps[0]}) == 0){
            $ps{$ps[0]} = [$ps[1]-1]; #-1 to index correctly on array
        } else {
            push( @{ $ps{$ps[0]} }, $ps[1]-1);
        }
    }
    return(\%ps, \@tigOrder);
}

sub process_fasta{
    my $fas = $_[0];
    open(FAS, "$fas") || die "Can't open fasta\n";
    chomp(my @fas = <FAS>);
    $fas = join(":::::", @fas);
    @fas = split(">", $fas);
    shift @fas;
    my %fas;
    foreach my $line (@fas){
        my @line = split(":::::", $line);
        my $head = shift(@line);
        $fas{$head} = join("", @line);
    }
    return(\%fas);
}

sub fas_slice{
    my($fas, $ps_ref, $order_ref) = @_;
    my $fas_ref = process_fasta($fas);
    my %fas = %$fas_ref;
    my %ps = %$ps_ref;
    my @tigOrder = @$order_ref;
    foreach my $tig (@tigOrder){
        my @seq = split("", $fas{$tig});
        my @slice = @seq[@{ $ps{$tig} }];
        print ">$tig\n", join("", @slice), "\n";
    }
}

#MAIN
{
    my($fas, $pos) = @ARGV;
    my($ps_ref, $order_ref) = hash_pos($pos);
    fas_slice($fas, $ps_ref, $order_ref);
}
