import pixy.core
import allel
import numpy as np
import warnings

from scipy import special
from itertools import combinations
from collections import Counter

g = allel.GenotypeArray( [
        [[0], [0]],
        [[1], [1]],
        [[0], [1]],
        [[0], [0]],
        [[1], [1]],
        [[0], [-1]],
        [[1], [-1]],
        [[-1], [-1]]
    ]
)

allele_counts = g.count_alleles(max_allele = 1)
variant_counts = allele_counts[allele_counts[:,1] != 0]
S = Counter(variant_counts[:,0] + variant_counts[:,1])
N = Counter(allele_counts[:,0] + allele_counts[:,1])


watterson_theta = 0
for n, s in S.items():
    a1 = np.sum(1 / np.arange(1, n))
    watterson_theta += s/a1


#check for singletons
variant_counts = variant_counts[variant_counts[:,0] + variant_counts[:,1] > 1]
S = Counter(variant_counts[:,0] + variant_counts[:,1])

watterson_theta = 0
for n, s in S.items():
    a1 = np.sum(1 / np.arange(1, n))
    watterson_theta += s/a1


N = Counter(allele_counts[:,0] + allele_counts[:,1])
N_array = np.array(tuple(N.items()))
#check for singleons
N_array = N_array[N_array[:,0] > 1]
weighted_sites = np.sum(np.multiply(N_array[:,1], (N_array[:,0]/max(N))))
