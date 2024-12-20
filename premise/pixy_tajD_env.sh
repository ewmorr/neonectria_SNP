
module purge
module load anaconda/colsa
conda create -n pixy_tajD

git clone https://github.com/npb596/pixy repo/pixy_tajD_git_clone
cd repo/pixy_tajD_git_clone/
git checkout -b fix_bugs
git branch

conda activate pixy_tajD
mamba install pip
pip install -e ~/repo/pixy_tajD_git_clone

mamba install -c bioconda htslib
mamba install -c bioconda samtools=1.9 --force-reinstall -y

#Had to indent line 42 of core X0
