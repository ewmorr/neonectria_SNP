
### stand-alone install (use when root priviledges are unavailable)
```
conda create --name pixy_tajD --clone template #on premise clone the template
conda activate pixy_tajD
conda install python=3.8
conda install -c bioconda htslib
conda install -c bioconda samtools=1.9 --force-reinstall -y
conda install pip
```
clone pixy from https://github.com/npb596/pixy/tree/master (the TajD development branch)
```
cd ~/repo
mkdir pixy_tajD
cd pixy_tajD
git clone https://github.com/npb596/pixy.git .
git branch
#pip install --user -e .
python -m pip install --no-build-isolation --no-deps --user -e .
python -m pip install scikit-allel
```
DEPRECATION: Legacy editable install of pixy==1.2.6b1 from file:///mnt/gpfs01/home/garnas/ewj4/repo/pixy_tajD (setup.py develop) is deprecated. pip 25.0 will enforce this behaviour change. A possible replacement is to add a pyproject.toml or enable --use-pep517, and use setuptools >= 64. If the resulting installation is not behaving as expected, try using --config-settings editable_mode=compat. Please consult the setuptools documentation for more information. Discussion can be found at https://github.com/pypa/pip/issues/11457
Running setup.py develop for pixy
Successfully installed pixy

## Note that the pixy executable from this install is in ~/.local/bin and is available to $PATH until activating the vanilla conda env



conda create -n pixy_tajD
conda activate pixy_tajD
cd ~/repo
mkdir pixy_tajD
cd pixy_tajD
git clone https://github.com/npb596/pixy.git .

conda install --file conda.recipe/
