#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output Nc_buscos.out
#SBATCH --exclude=node117,node118

module purge
#module load anaconda/colsa
#conda activate busco-5.4.7
module load linuxbrew/colsa


#export AUGUSTUS_CONFIG_PATH="/mnt/home/garnas/ewj4/augustus_config/config"

cd $HOME/Nf_Nd_Nc_buscos
#The lineage is copied from augustus config (or downloaded from website) so that augustus can write to this dir. Write-access to config dir was required when running with --long option with older version but may not be needed with newer version of BUSCO. Have not tested, this solution works either way.


#busco -i $HOME/Nc_SPANDx_all_seqs/ref.fasta -o Nc --download_path /mnt/gpfs01/home/garnas/ewj4/busco_downloads --augustus -l hypocreales_odb10 -m genome -c 24 -f --offline
busco -i $HOME/Nc_SPANDx_all_seqs/ref.fasta -l $HOME/augustus_config/lineage/sordariomyceta_odb9 -m genome -o Nf_buscos -c 1 -f
