#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output Nf_buscos.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd $HOME/neonectria_minion/MAT1_polish_2
#The lineage is copied from augustus config (or downloaded from website) so that augustus can write to this dir. Write-access to config dir was required when running with --long option with older version but may not be needed with newer version of BUSCO. Have not tested, this solution works either way.

srun busco -i pilon_.fasta -l $HOME/augustus_config/lineage/sordariomyceta_odb9 -m genome -o Nf_buscos -c 1 -f

#note, the single copy busco seq are in dir single_copy_busco_sequences and the list of which BUSCO is single copy, dup'd, fragmented, with contigs and positions is in "full_table_...buscos.tsv". Can either use the positions file to pull out from the three species' SNP tables are just rerun SNP calling against the set of SCOs. The latter seems less error prone towards having a proper alignment. We just need to make sure we exclude dups (and missing, but these would just be missing data). Although, I guess the quality filtering would need to be redone if we recall.... that's probably fine bc it's just one file pass and we have the algorithm already set
