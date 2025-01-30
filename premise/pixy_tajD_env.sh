
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
conda activate pixy_tajD
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

g = allel.GenotypeArray( [
        [[0], [0]],
        [[0], [1]],
        [[1], [1]],
        [[0], [-1]],
        [[1], [-1]],
        [[-1], [-1]]
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
2 0
1 1
0 2
1 0
0 1
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
#allele_counts = allele_counts[allele_counts[:,0]+allele_counts[:,1] > 1]


variant_counts = allele_counts[allele_counts[:,1] != 0]
variant_counts = variant_counts[variant_counts[:,0] + variant_counts[:,1] > 1]

S = Counter(variant_counts[:,0] + variant_counts[:,1])
N = Counter(allele_counts[:,0] + allele_counts[:,1])

N_array = np.array(tuple(N.items()))
#the following works to remove single genotype sites from the weighted site matrix (which should not be used to weight the watterson's theta calc)
N_array = N_array[N_array[:,0] != 1]
weighted_sites = np.sum(np.multiply(N_array[:,1], (N_array[:,0]/max(N))))

watterson_theta = 0
for n, s in S.items():
    a1 = np.sum(1 / np.arange(1, n))
    watterson_theta += s/a1


for n, s in S.items():
    print(n)
    print(s)

##################################
# the error looks like it is coming from np.arrange
# where n == 1. This results in a 0 float from the sum and hence illegal division
# change the declaration of variant counts to only consider sites where allele_counts[:,0] + allele_counts[:,1] > 1

# we also need to incorporate these new checks into the core.py script so that the correct number of total sites and seg sites will be reported (for incorporation into genome wide clacs)
line 514 and L515 (allele_counts and variant_sites in waterson theta def)
line 550 (allele_counts tajima's d)
