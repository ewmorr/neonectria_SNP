#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output Nf_buscos.out
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa
conda activate busco-5.4.7

export AUGUSTUS_CONFIG_PATH="/mnt/home/garnas/ewj4/augustus_config/config"

cd $HOME/Nf_Nd_Nc_buscos
#The lineage is copied from augustus config (or downloaded from website) so that augustus can write to this dir. Write-access to config dir was required when running with --long option with older version but may not be needed with newer version of BUSCO. Have not tested, this solution works either way.


busco -i $HOME/Nf_SPANDx_all_seqs/ref.fasta -o Nf --download_path /mnt/gpfs01/home/garnas/ewj4/busco_downloads --augustus -l hypocreales_odb10 -m genome -c 24 -f --offline
