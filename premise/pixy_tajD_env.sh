
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
conda install scikit-allel

#Had to indent line 42 of core.py X0

#testing zero divide error inpython
# note rows (outer diminsion) equal variants and columns (second dimension) equal samples
# the las t dimension indicates ploidy which we have coerced to diploid
# we need to test what happens when a variant is all NA
# we are also biallelic and strictly homozygous
>>>
import pixy.core
import allel
import numpy as np
import warnings

from scipy import special
from itertools import combinations
from collections import Counter

g = allel.GenotypeArray( [
        [[0, 0], [0, 0]],
        [[0, 0], [1, 1]],
        [[1, 1], [1, 1]],
        [[0, 0], [-1, -1]],
        [[1, 1], [-1, -1]],
        [[-1, -1], [-1, -1]]
    ]
)

#count alelles
g.count_alleles()
# <AlleleCountsArray shape=(6, 2) dtype=int32>
4 0
2 2
0 4
2 0
0 2
0 0
>>>

#when the site is all NA  count_alleles returns a zero zero value. The current code only tests for
# incidence of 0 in the second position (ie., the variant) of the variant_sites calli
# we can fix by testing for the two indices summing to zero

# DO WE NEED TO EXCLUDE THE 0/0 ALLELE COUNTS FROM OTHER PARTS OF THE CALCS?!?!?!?!?!

# Moved the test for zerozero to the definition of the alelle counts arraoy to exclude sites
# with all missing values from the jump
# this dealt with a bunhc of the div by Zero errors buit there is still one

# new/remaining error is casued by
# watterson_theta += s/a1

#testing on inteactive python env with g array above
allele_counts = g.count_alleles(max_allele = 1)
allele_counts = allele_counts[allele_counts[:,0]+allele_counts[:,1] != 0]


variant_counts = allele_counts[allele_counts[:,1] != 0]

S = Counter(variant_counts[:,0] + variant_counts[:,1])
N = Counter(allele_counts[:,0] + allele_counts[:,1])

N_array = np.array(tuple(N.items()))

watterson_theta = 0
for n, s in S.items():
    a1 = np.sum(1 / np.arange(1, n))
    watterson_theta += s/a1


for n, s in S.items():
    print(n)
    print(s)
