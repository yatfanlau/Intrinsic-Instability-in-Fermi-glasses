import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# System constants
L = 60
W = 8
tp = 0.4
no_of_disorder = 15
Loc = 7.516045912

def generate_disorder(L, W):
    """Generate a matrix of uniform random disorder."""
    return W * (np.random.uniform(size=(L, L)) - 0.5)

def generate_hamiltonian(L, W):
    """Generate the Hamiltonian matrix for a 2D lattice with periodic boundary conditions."""
    H = np.zeros((L * L, L * L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        for j in range(L):
            idx = i * L + j
            H[idx, idx] = disorder[i, j]
            # Apply periodic boundary conditions
            neighbors = [(i, (j + 1) % L), ((i + 1) % L, j), (i, (j - 1) % L), ((i - 1) % L, j)]
            for ni, nj in neighbors:
                nidx = ni * L + nj
                H[idx, nidx] = 1.0
                H[nidx, idx] = 1.0
                # Diagonal hopping
                if (ni, nj) != ((i - 1) % L, j) and (ni, nj) != (i, (j - 1) % L):
                    H[idx, nidx] = tp
                    H[nidx, idx] = tp
    return H

def find_max(index, eigenstates, L):
    """Find the position of the maximum amplitude in an eigenstate."""
    eigenstate = eigenstates[:, index]
    position = np.argmax(eigenstate)
    x, y = position % L, position // L
    return x + 1, y + 1

def distance(x1, y1, x2, y2):
    """Calculate normalized Euclidean distance."""
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2) / Loc

# Simulation of disorder realizations and energy level spacings
Y = [[] for _ in range(10)]
X_array = np.array([3, 6, 9, 12, 15, 18, 21, 24, 27, 30]) / Loc

for _ in range(no_of_disorder):
    H = generate_hamiltonian(L, W)
    energy_levels, eigenstates = linalg.eigh(H)
    eigenstates2 = eigenstates ** 2

    energy_spacing = [[] for _ in range(10)]
    for m, Ld in enumerate(X_array):
        for s in range(L * L):
            x, y = find_max(s, eigenstates2, L)
            if distance(30, 30, x, y) < Ld:
                energy_spacing[m].append(energy_levels[s])

        if len(energy_spacing[m]) > 1:
            energy_spacing[m].sort()
            Y[m].extend(np.diff(energy_spacing[m]))

# Fit and plot data
def fit_func(x, A, B):
    """Model function for curve fitting."""
    return A / (B + x**2)

avg_spacing = [np.mean(spacings) for spacings in Y if spacings]
X_array *= 2
parameters, _ = curve_fit(fit_func, X_array, avg_spacing)
fit_A, fit_B = parameters

x_fit = np.linspace(0.8, 8.5, 100)
y_fit = fit_func(x_fit, fit_A, fit_B)

plt.plot(x_fit, y_fit, label=r"$\langle \Delta E \rangle = \frac{A}{B + l^2}$", color="blue", linestyle="solid")
plt.scatter(X_array, avg_spacing, marker="^", color="red")
plt.xlabel("Size of the grain (l)", fontsize=12)
plt.ylabel("Average energy level spacing ($\langle \Delta E \rangle$)", fontsize=12)
plt.legend()
plt.show()
