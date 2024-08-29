import numpy as np
from scipy.linalg import eigh
from random import randint

# Constants for the system size and parameters
L = 60       # System size
W = 14       # Disorder strength
tp = 0.4     # Next-nearest neighbor hopping term

def generate_disorder(L, W):
    """ Generate an LxL matrix of uniform random disorder. """
    return W * (np.random.uniform(size=(L * L)).reshape((L, L)) - 0.5)

def generate_hamiltonian(L, W):
    """ Generate the Hamiltonian matrix for a 2D lattice with periodic boundary conditions. """
    H = np.zeros((L * L, L * L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        ip1 = (i + 1) % L
        im1 = (i - 1) % L
        for j in range(L):
            idx = i * L + j
            H[idx, idx] = disorder[i, j]
            jp1 = (j + 1) % L
            jm1 = (j - 1) % L
            # Nearest neighbor hopping
            H[idx, ip1 * L + j] = H[ip1 * L + j, idx] = 1.0
            H[idx, i * L + jp1] = H[i * L + jp1, idx] = 1.0
            # Next-nearest neighbor hopping
            H[idx, ip1 * L + jp1] = H[ip1 * L + jp1, idx] = tp
            H[idx, ip1 * L + jm1] = H[ip1 * L + jm1, idx] = tp
    return H

def center_of_mass(index):
    """ Calculate the center of mass of an eigenstate. """
    s, r = 0, 0
    for i in range(L * L):
        x, y = i % L, i // L
        s += x * eigenstates2[i, index]
        r += y * eigenstates2[i, index]
    return s, r

def find_max(index):
    """ Find the position of the maximum amplitude in an eigenstate. """
    eigenstate = eigenstates2[:, int(index)]
    position = np.argmax(eigenstate)
    x, y = position % L, position // L
    return x, y

# Generate the Hamiltonian and solve for eigenstates
H = generate_hamiltonian(L, W)
energy_levels, eigenstates = eigh(H)
eigenstates2 = eigenstates ** 2

# Example use of center_of_mass and find_max for random eigenstates
for i in range(20):
    a = randint(3000, 3400)
    com = center_of_mass(a)
    max_pos = find_max(a)
    print(f"Center of Mass for state {a}: {com}")
    print(f"Max position for state {a}: {max_pos}")
