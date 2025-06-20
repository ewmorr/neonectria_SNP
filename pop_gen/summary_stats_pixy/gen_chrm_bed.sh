#!/bin/bash

gunzip -c ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/FINAL_invariant.nuclear.vcf.gz | grep "##contig=<ID=" | sed 's/##contig=<ID=//' | sed 's/length=//' | sed 's/>$//' | perl -pe 's/,/\t/' > ~/repo/neonectria_SNP/data/Nf/pixy/contig_lens.tsv

gunzip -c ~/repo/neonectria_SNP/data/Nd/final_tables/rm_dups/FINAL_invariant.nuclear.vcf.gz | grep "##contig=<ID=" | sed 's/##contig=<ID=//' | sed 's/length=//' | sed 's/>$//' | perl -pe 's/,/\t/' > ~/repo/neonectria_SNP/data/Nd/pixy/contig_lens.tsv

gunzip -c ~/repo/neonectria_SNP/data/Nc/final_tables/FINAL_invariant.nuclear.vcf.gz | grep "##contig=<ID=" | sed 's/##contig=<ID=//' | sed 's/length=//' | sed 's/>$//' | perl -pe 's/,/\t/' > ~/repo/neonectria_SNP/data/Nc/pixy/contig_lens.tsv
