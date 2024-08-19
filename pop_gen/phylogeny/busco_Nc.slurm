#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output Nc_buscos.out
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa
conda activate busco-5.4.7

cd $HOME/Nf_Nd_Nc_buscos
#The lineage is copied from augustus config (or downloaded from website) so that augustus can write to this dir. Write-access to config dir was required when running with --long option with older version but may not be needed with newer version of BUSCO. Have not tested, this solution works either way.

srun busco -i $HOME/Nc_SPANDx_all_seqs/ref.fasta -o Nc \
--augustus -l hypocreales_odb10 -m genome -c 24 -f
