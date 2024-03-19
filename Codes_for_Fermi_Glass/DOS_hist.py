from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


L = 60  # int(input('System size: '))
W = 0  # float(input('Disorder strength: '))
tp = 0.4


def generate_disorder(L, W):
    disorder = W * ((np.random.uniform(size=L * L)).reshape((L, L)) - 0.5)
    return disorder


def generate_hamiltonian(L, W):
    H = np.zeros((L * L, L * L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        ip1 = (i + 1) % L
        im1 = (i - 1) % L
        for j in range(L):
            H[i * L + j, i * L + j] = disorder[i, j]
            jp1 = (j + 1) % L
            jm1 = (j - 1) % L
            H[ip1 * L + j, i * L + j] = 1.0
            H[i * L + j, ip1 * L + j] = 1.0
            H[i * L + jp1, i * L + j] = 1.0
            H[i * L + j, i * L + jp1] = 1.0

            H[i * L + j, ip1 * L + jp1] = tp
            H[ip1 * L + jp1, i * L + j] = tp
            H[i * L + j, ip1 * L + jm1] = tp
            H[ip1 * L + jm1, i * L + j] = tp
    return H

E = [[]for i in range(L*L)]
Energy = []
for con in range(15):   
    H = generate_hamiltonian(L,W)
    (energy_levels,eigenstates)=linalg.eigh(H)
    for i in range(L*L):
        E[i].append(energy_levels[i])

for k in range(L*L):
    avg = sum(E[k])/len(E[k])
    Energy.append(avg)
    

n, bins, patches = plt.hist(Energy, bins = 50, histtype='step',density=True)
plt.show()
bins = bins.tolist()
bins.remove(bins[50])
bins = np.array(bins)
plt.plot(bins, n)
plt.title('L='+str(L)+', '+'W='+str(W))
plt.show()