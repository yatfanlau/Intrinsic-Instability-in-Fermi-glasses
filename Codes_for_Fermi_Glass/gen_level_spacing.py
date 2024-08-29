import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson
import proplot

# Constants
L = 60
W = 9
Loc = 7.65151489
tp = 0.6
R = 0.5  # Radius of the grain

def generate_disorder(L, W):
    """Generate a disorder matrix for the Hamiltonian."""
    disorder = W * (np.random.uniform(size=(L, L)) - 0.5)
    return disorder

def generate_hamiltonian(L, W):
    """Generate the Hamiltonian matrix for a 2D lattice with periodic boundary conditions."""
    H = np.zeros((L*L, L*L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        for j in range(L):
            # Index mapping
            index = i * L + j
            H[index, index] = disorder[i, j]
            neighbors = [(i, (j+1) % L), ((i+1) % L, j), (i, (j-1) % L), ((i-1) % L, j), 
                         ((i+1) % L, (j+1) % L), ((i+1) % L, (j-1) % L)]
            for ni, nj in neighbors:
                neighbor_index = ni * L + nj
                H[index, neighbor_index] = tp if ni != i and nj != j else 1.0
    return H

def find_max(eigenstates, index, L):
    """Find the position of maximum probability in the eigenstate."""
    eigenstate = eigenstates[:, index]**2
    max_index = np.argmax(eigenstate)
    x, y = max_index % L, max_index // L + 1
    return x + 1, y if x == L - 1 else y

def distance(x1, y1, x2, y2):
    """Calculate normalized distance."""
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2) / Loc

# Main execution
energy_spacing = []
for con in range(30):
    H = generate_hamiltonian(L, W)
    energy_levels, eigenstates = linalg.eigh(H)
    energy = []
    for s in range(L * L):
        x, y = find_max(eigenstates, s, L)
        if distance(30, 30, x, y) < R:
            energy.append(energy_levels[s])
    energy.sort()
    energy_spacing.extend(np.diff(energy))

# Saving the energy level spacings
np.save(f"l={2 * R}", energy_spacing)
