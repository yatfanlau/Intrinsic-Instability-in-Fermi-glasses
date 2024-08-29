import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Constants
L = 60  # System size
W = 8   # Disorder strength
Loc = 7.516045912  # Localization length
tp = 0.4  # Hopping term coefficient
Ld = 1.5  # Grain size (in terms of localization length)

def generate_disorder(L, W):
    """Generate a disorder matrix."""
    return W * (np.random.uniform(size=(L, L)) - 0.5)

def generate_hamiltonian(L, W):
    """Generate the Hamiltonian with given disorder and hopping terms."""
    H = np.zeros((L * L, L * L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        for j in range(L):
            index = i * L + j
            H[index, index] = disorder[i, j]
            for di, dj, coupling in [((1, 0), 1.0), ((0, 1), 1.0), ((1, 1), tp), ((1, -1), tp)]:
                ni, nj = (i + di) % L, (j + dj) % L
                neighbor_index = ni * L + nj
                H[index, neighbor_index] = coupling
                H[neighbor_index, index] = coupling
    return H

def find_max(eigenstates, index, L):
    """Find the position of maximum probability."""
    eigenstate = eigenstates[:, index]**2
    max_index = np.argmax(eigenstate)
    x, y = max_index % L, max_index // L + 1
    return x, y

def distance(x1, y1, x2, y2):
    """Compute normalized distance between two points."""
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2) / Loc

def poisson_fitting_function(x, h):
    """Poisson-like fitting function."""
    return np.exp(-x / h) / h

# Simulation
energy_spacing = []
for _ in range(15):
    H = generate_hamiltonian(L, W)
    energy_levels, eigenstates = linalg.eigh(H)
    energy = []
    for s in range(L * L):
        x, y = find_max(eigenstates, s, L)
        if distance(30, 30, x, y) < Ld:
            energy.append(energy_levels[s])
    energy.sort()
    energy_spacing.extend(np.diff(energy))

# Histogram and fitting
BIN = int(np.sqrt(len(energy_spacing))) + 1
n, bins, patches = plt.hist(energy_spacing, bins=BIN, density=True, color="b", alpha=0.7)

# Fit using the middle of bins as x values
bin_centers = 0.5 * (bins[:-1] + bins[1:])
parameters, covariance = curve_fit(poisson_fitting_function, bin_centers, n)
fit_h = parameters[0]

x_fit = np.linspace(min(energy_spacing), max(energy_spacing), 200)
y_fit = poisson_fitting_function(x_fit, fit_h)

plt.plot(x_fit, y_fit, 'r-', label=f"Poisson fit: h={fit_h:.3f}")
plt.title(f'Grain size={Ld*2}, Δ={fit_h:.3f}')
plt.xlabel("Energy level spacing / Δ")
plt.ylabel("Normalized frequency")
plt.legend()
plt.show()
