
#run on premise
cd ~/Nf_Nd_Nc_buscos
grep Complete ~/neonectria_minion/MAT1_polish_2/run_Nf_buscos/full_table_Nf_buscos.tsv | cut -f 1 | sort > Nf_complete_buscos.txt
grep Complete run_Nc_buscos/full_table_Nc_buscos.tsv | cut -f 1 | sort > Nc_complete_buscos.txt
grep Complete run_Nd_buscos/full_table_Nd_buscos.tsv | cut -f 1 | sort > Nd_complete_buscos.txt

wc -l Nf_complete_buscos.txt
#3665
wc -l Nc_complete_buscos.txt
#3651
wc -l Nd_complete_buscos.txt
#3655

comm -12 Nf_complete_buscos.txt Nd_complete_buscos.txt | wc -l
#3627
comm -12 Nf_complete_buscos.txt Nc_complete_buscos.txt | wc -l
# 3629
comm -12 Nd_complete_buscos.txt Nc_complete_buscos.txt | wc -l
#3612

comm -12 Nf_complete_buscos.txt Nd_complete_buscos.txt | comm -12 - Nc_complete_buscos.txt | wc -l
#3597

comm -12 Nf_complete_buscos.txt Nd_complete_buscos.txt | comm -12 - Nc_complete_buscos.txt > shared_busco_IDs.txt

busco_dir=~/neonectria_minion/MAT1_polish_2/run_Nf_buscos/single_copy_busco_sequences/
mkdir shared_busco_seqs

# a few of the BUSCO sequences have more than one sequence in the file, and they are not the correct length
# based on the position output in the fasta headers (which are all the same)
# E.g.,
perl ~/perl_scripts/get_seq_lens.pl $busco_dir/EOG09331316.fna
#check Nc and Nd
perl ~/perl_scripts/get_seq_lens.pl run_Nc_buscos/single_copy_busco_sequences/EOG09331316.fna
# this length is also incorrect based on the seq header, length is actually 1005 bp but positions suggest 1162. Not clear how busco/augustus is
# handling gaps and such. But at least there is only one seq # len based on positions is 1162, dif is 157
perl ~/perl_scripts/get_seq_lens.pl run_Nd_buscos/single_copy_busco_sequences/EOG09331316.fna
#same deal here #142 dif

#it looks like Augustus is outputting the CDS so there may different seqs dependent on any predicted introns. We should exclude any sequences with multiple CDS. Should we check this for all the spp? Probably yes
# the BWA alignment should be fine predicting gaps, and we can then exclude them in the final tree calc, so just use the single CDS Nf seqs
# BUT we check for single CD below anyways and use those

while read i
do(
    echo $i >> seqs_per_gene.txt
    numSeqsNf=$(grep ">" $busco_dir/$i.fna | wc -l)
    numSeqsNc=$(grep ">" run_Nc_buscos/single_copy_busco_sequences/$i.fna | wc -l)
    numSeqsNd=$(grep ">" run_Nd_buscos/single_copy_busco_sequences/$i.fna | wc -l)
    if [ $numSeqsNf -eq 1 ] && [ $numSeqsNc -eq 1 ] && [ $numSeqsNd -eq 1 ]
    then(
            #sed 's/pilon_.fasta:/Nf:/' $busco_dir/$i.fna >> shared_buscos.fasta
            cp $busco_dir/$i.fna shared_busco_seqs
    )fi
    #sed 's/pilon_.fasta:/Nf:/' $busco_dir/$i.fna >> shared_buscos.fasta
    #cp $busco_dir/$i.fna shared_busco_seqs
)done < shared_busco_IDs.txt

ls shared_busco_seqs/ | wc -l
#3592

cat shared_busco_seqs/* > shared_buscos.fasta
sed 's/pilon_.fasta://g' shared_buscos.fasta > shared_busco.fasta.tmp & mv shared_busco.fasta.tmp shared_buscos.fasta
grep ">" shared_buscos.fasta | wc -l
#3592
less shared_buscos.fasta

cd
mkdir Nf_Nd_Nc_buscos_SPANDx
cp Nf_Nd_Nc_buscos/shared_buscos.fasta Nf_Nd_Nc_buscos_SPANDx/ref.fasta
cd Nf_Nd_Nc_buscos_SPANDx

cp ../Nf_SPANDx_all_seqs/*fastq.gz . &
cp ../Nd_SPANDx_all_seqs/*fastq.gz . &
cp ../Nc_SPANDx_all_seqs/*fastq.gz . &

#rm all the samples that were filtered due to missing ness in the full SNP dataset
rm NG10*
rm NG55*
rm NG62*
rm NG28*
rm NG48*
rm NG139*
rm NG161*
rm NG130*
rm NG138*
