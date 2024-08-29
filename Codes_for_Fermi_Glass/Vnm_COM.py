import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from random import randint

# Constants
L = 60
W = 10
tp = 0.4
DELTA = 100
E_CENTER = 3400
N = 30
NUM_CONFIGURATIONS = 10

def generate_disorder(L, W):
    """ Generate a disorder matrix for the Hamiltonian. """
    return W * (np.random.uniform(size=(L, L)) - 0.5)

def generate_hamiltonian(L, W):
    """ Generate the Hamiltonian matrix for the system. """
    H = np.zeros((L*L, L*L))
    disorder = generate_disorder(L, W)
    for i in range(L):
        for j in range(L):
            index = i * L + j
            right = i * L + (j + 1) % L
            down = ((i + 1) % L) * L + j
            diag_right_down = ((i + 1) % L) * L + (j + 1) % L
            diag_right_up = ((i + 1) % L) * L + (j - 1) % L

            H[index, index] = disorder[i, j]
            H[index, right] = H[right, index] = 1.0
            H[index, down] = H[down, index] = 1.0
            H[index, diag_right_down] = H[diag_right_down, index] = tp
            H[index, diag_right_up] = H[diag_right_up, index] = tp

    return H

def center_of_mass(state, L):
    """ Calculate the center of mass for a given state. """
    grid = np.arange(L*L)
    x = grid % L
    y = grid // L
    cm_x = np.sum(x * state) / np.sum(state)
    cm_y = np.sum(y * state) / np.sum(state)
    return cm_x, cm_y

def distance(x1, y1, x2, y2):
    """ Compute Euclidean distance between two points. """
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

def overlap(state1, state2):
    """ Compute the overlap integral between two states. """
    return np.dot(state1, state2)

# Main data collection
d_list = []
v_list = []

for _ in range(NUM_CONFIGURATIONS):
    H = generate_hamiltonian(L, W)
    energies, eigenstates = eigh(H)
    probabilities = eigenstates**2

    selected_indices = np.random.choice(np.where((energies > E_CENTER - DELTA) & (energies < E_CENTER + DELTA))[0], 2*N, replace=False)

    for k in selected_indices:
        for m in selected_indices:
            if k != m:
                xk, yk = center_of_mass(probabilities[:, k], L)
                xm, ym = center_of_mass(probabilities[:, m], L)
                dist = distance(xk, yk, xm, ym) / 6
                if dist < 5:
                    vnm = overlap(probabilities[:, k], probabilities[:, m])
                    d_list.append(dist)
                    v_list.append(vnm)

# Regression and plotting
d_array = np.array(d_list)
v_array = np.array(v_list)
ln_v_array = np.log(v_array)
A, B = np.polyfit(d_array, ln_v_array, 1, w=np.sqrt(v_array))

x = np.linspace(0, 5, 100)
y = np.exp(B) * np.exp(A * x)
plt.plot(x, y, label=f'U_kp = exp({B:.2f})exp({A:.2f}d)')
plt.scatter(d_list, v_list, s=3, color='r')
plt.xlabel("Distance between two states (In units of localization length)")
plt.ylabel("U_kp")
plt.title(f'Total pairs = {len(d_list)}, {L}x{L}, W={W}, Configurations = {NUM_CONFIGURATIONS}')
plt.legend()
plt.show()
