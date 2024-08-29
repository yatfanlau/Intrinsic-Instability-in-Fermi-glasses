import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from random import randint

L = 60
W = 9
Loc = 6.02516045912
tp = 0.6

def generate_disorder(L, W):
    return W * (np.random.uniform(size=(L, L)) - 0.5)

def generate_hamiltonian(L, W):
    H = np.zeros((L*L, L*L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        for j in range(L):
            H[i*L+j, i*L+j] = disorder[i, j]
            for di, dj in [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (-1, 1), (1, -1), (-1, -1)]:
                ni, nj = (i + di) % L, (j + dj) % L
                H[i*L+j, ni*L+nj] = tp if abs(di + dj) == 2 else 1.0
    return H

def find_max_position(eigenstate):
    index = np.argmax(eigenstate)
    x, y = index % L, index // L
    return x + 1, y + 1  # Convert to 1-based index

def distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2) / Loc

def overlap(eigenstates, p, q):
    return np.dot(eigenstates[:, p], eigenstates[:, q])

# Constants for the computation
delta = 200
E1, E2 = 1800 - delta, 1800 + delta
Range = 0.1
d_list = np.array([7, 12, 17, 22, 27]) / Loc

# Simulation
final_Vnm = []
distance_list = []
data_Vnm = []
disorder_Vnm = [[] for _ in range(5)]

for con in range(5):
    H = generate_hamiltonian(L, W)
    energies, eigenstates = eigh(H)
    eigenstates2 = eigenstates**2
    for _ in range(2000):
        i, j = randint(E1, E2), randint(E1, E2)
        if i != j:
            x1, y1 = find_max_position(eigenstates2[:, i])
            x2, y2 = find_max_position(eigenstates2[:, j])
            dist = distance(x1, y1, x2, y2)
            if dist < 27 / Loc:
                distance_list.append(dist)
                V = overlap(eigenstates, i, j)
                for k, d_center in enumerate(d_list):
                    if abs(dist - d_center) < Range:
                        disorder_Vnm[k].append(V)
                        data_Vnm.append(V)
                        break

for values in disorder_Vnm:
    final_Vnm.append(np.mean(values))

# Fitting
def linear_func(x, m, c):
    return m * x + c

Vnm_ = np.log(final_Vnm)
params, _ = curve_fit(linear_func, d_list, Vnm_)
fit_m, fit_c = params

# Plotting
x_ = np.linspace(0, 27 / Loc, 100)
y_ = fit_m * x_ + fit_c
plt.plot(x_, y_, linestyle="solid", color="b", label="Fitting result: y=mx+c")
plt.scatter(d_list, Vnm_, s=10, color="red", label="Simulation result")
plt.xlabel("$Distance$ (in unit of localization length)")
plt.ylabel("$\ln(U_{kp})$")
plt.title(f"m={fit_m:.2f}, c={fit_c:.2f}")
plt.legend()
plt.show()
