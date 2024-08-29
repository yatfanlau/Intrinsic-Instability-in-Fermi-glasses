from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

L = 60
tp = -0.6
t = -1

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
            H[ip1 * L + j, i * L + j] = t
            H[i * L + j, ip1 * L + j] = t
            H[i * L + jp1, i * L + j] = t
            H[i * L + j, i * L + jp1] = t
            H[i * L + j, ip1 * L + jp1] = tp
            H[ip1 * L + jp1, i * L + j] = tp
            H[i * L + j, ip1 * L + jm1] = tp
            H[ip1 * L + jm1, i * L + j] = tp
    return H

Energy = []
for w in range(1):
    W = 9
    E = [[] for i in range(L * L)]

    for con in range(15):
        H = generate_hamiltonian(L, W)
        (energy_levels, eigenstates) = linalg.eigh(H)
        for i in range(L * L):
            E[i].append(energy_levels[i])

    for k in range(L * L):
        avg = sum(E[k]) / len(E[k])
        Energy.append(avg)
with open(f'eigW{W}.txt', 'w') as f:
    f.write("Eigenvalues\n")  # Column name as the header
    np.savetxt(f, Energy)


# style = ["solid", "dashed", "dashdot", (0, (3, 1, 1, 1, 1, 1))]
# Disorder = [0, 4, 8, 12]

# for w in range(4):
#     W = Disorder[w]
#     E = [[] for i in range(L * L)]
#     Energy = []
#     for con in range(5):
#         H = generate_hamiltonian(L, W)
#         (energy_levels, eigenstates) = linalg.eigh(H)
#         for i in range(L * L):
#             E[i].append(energy_levels[i])

#     for k in range(L * L):
#         avg = sum(E[k]) / len(E[k])
#         Energy.append(avg)

#     n, Bins, _ = plt.hist(Energy, bins=10, density=True)
#     plt.xlim(-3, 4)
#     plt.plot(Bins[:-1], n, linestyle=style[w], label="W=" + str(W))

# plt.xlabel("Energy")
# plt.ylabel("DOS")
# plt.legend(loc='upper left')
# plt.show()