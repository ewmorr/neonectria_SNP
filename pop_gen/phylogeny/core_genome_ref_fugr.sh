conda activate bcftools

cd repo/neonectria_SNP/data/Fugr1_ref
#repeat masked Fusarium graminearum genome https://mycocosm.jgi.doe.gov/Fusgr1/Fusgr1.home.html
#using strategy described by https://phame.readthedocs.io/en/latest/introduction/runphame.html to identify and extract core genome using mummer whole-genome alignmnet
#nucmer --maxmatch refgenome.fasta genome.fasta

#the repeat masking is currently done by lowercase letters (as by convention). Convert to N for nucmer
sed -e 's/[actg]/N/g' Fusgr1_AssemblyScaffolds_Repeatmasked.fasta > Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta

sed -i -e 's/SuperNonNiN/Supercontig/g' Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta #lol

nucmer --maxmatch --prefix=Nf Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta ../Nf/ref.fasta
nucmer --maxmatch --prefix=Nd Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta ../Nd/ref.fasta
nucmer --maxmatch --prefix=Nc Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta ../Nc/ref.fasta

delta-filter -1 Nf.delta > Nf.one2oneFilter
delta-filter -1 Nd.delta > Nd.one2oneFilter
delta-filter -1 Nc.delta > Nc.one2oneFilter

show-coords -clTr Nf.one2oneFilter > Nf.coords
show-coords -clTr Nd.one2oneFilter > Nd.coords
show-coords -clTr Nc.one2oneFilter > Nc.coords

perl ../../library/collate_nucmer_coords.pl Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta Nf.coords Nd.coords Nc.coords
#334314 total positions in core alignment. This ios consistent with the phame paper which found 134062 positions and 40675 SNPs shared in an alignment of E coli Shigella and Salmonella (compared to 2.16M positions ansd 267K SNPs when just comparing E coli and Shigella
#could consider rerunning this with just the Neonectria

cd
cd repo/neonectria_SNP/data/Nf
# note that the repeatmasking proposed in phame doesn't actually do anything for Nf
#need to first run repeatmasker
