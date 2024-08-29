import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# Constants
L = 60
W = 8
tp = 0.4
Loc = 7.516045912

def generate_disorder(L, W):
    """Generate a disorder matrix for the Hamiltonian."""
    return W * (np.random.rand(L, L) - 0.5)

def generate_hamiltonian(L, W):
    """Generate the Hamiltonian matrix for a 2D lattice with periodic boundary conditions and next-nearest neighbor hopping."""
    H = np.zeros((L * L, L * L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        for j in range(L):
            index = i * L + j
            H[index, index] = disorder[i, j]
            for di, dj in [(0, 1), (1, 0), (0, -1), (-1, 0), (1, 1), (1, -1)]:
                ni, nj = (i + di) % L, (j + dj) % L
                neighbor_index = ni * L + nj
                H[index, neighbor_index] = tp if abs(di + dj) % 2 == 0 else 1.0
    return H

def find_max(eigenstate, L):
    """Find the coordinates of the maximum value in a 2D eigenstate."""
    max_index = np.argmax(eigenstate)
    x, y = max_index % L, max_index // L
    return x + 1, y + 1

def distance(x1, y1, x2, y2):
    """Compute normalized Euclidean distance."""
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2) / Loc

def calculate_energy_spacing(L, energy_levels, eigenstates, Ld):
    """Calculate energy spacings for states within a radius Ld."""
    energy = []
    for s in range(L * L):
        x, y = find_max(eigenstates[:, s], L)
        if distance(30, 30, x, y) < Ld:
            energy.append(energy_levels[s])
    energy.sort()
    return np.diff(energy)

# Initialize parameters
Y = [[] for _ in range(5)]
X = np.array([3, 6, 9, 12, 15, 18, 21, 24, 27, 30]) / Loc

# Main simulation loop
for conf in range(5):
    H = generate_hamiltonian(L, W)
    energy_levels, eigenstates = linalg.eigh(H)
    for Ld in X:
        energy_spacing = calculate_energy_spacing(L, energy_levels, eigenstates, Ld)
        Y[conf].append(np.mean(energy_spacing) if energy_spacing.size > 0 else 0)

# Plotting
markers = ["^", "s", "o", "*", "D"]
for conf in range(5):
    plt.scatter(X, Y[conf], s=30, marker=markers[conf], label=f'Config {conf + 1}')

plt.xlabel("Radius of the chosen region (normalized)")
plt.ylabel("Average energy level spacing")
plt.legend()
plt.show()
