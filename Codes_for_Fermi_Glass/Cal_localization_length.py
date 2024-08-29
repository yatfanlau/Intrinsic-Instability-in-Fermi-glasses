import numpy as np
from scipy.linalg import eigh

# Constants for the system size, disorder strength, and hopping term
L = 60
W = 9
tp = 0.6

def generate_disorder(L, W):
    """Generate a LxL matrix of uniform random disorder values."""
    return W * (np.random.uniform(size=(L, L)) - 0.5)
  
def generate_hamiltonian(L, W):
    """Generate the Hamiltonian for a 2D lattice with disorder and periodic boundary conditions."""
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
            # Setting up nearest neighbor hopping terms
            H[idx, ip1 * L + j] = 1.0
            H[idx, i * L + jp1] = 1.0
            H[i * L + jp1, idx] = 1.0
            H[ip1 * L + j, idx] = 1.0
            # Setting up next-nearest neighbor hopping terms
            H[idx, ip1 * L + jp1] = tp
            H[ip1 * L + jp1, idx] = tp
            H[idx, ip1 * L + jm1] = tp
            H[ip1 * L + jm1, idx] = tp
    return H

def cal_loc_length(eigenstates):
    """Calculate the average localization length from the eigenstates."""
    Loc_list = []
    for k in range(1600, 2000):
        IPR = np.dot(eigenstates[:, k], eigenstates[:, k])
        Loc = np.sqrt(1 / IPR)
        Loc_list.append(Loc)
    avg_Loc = sum(Loc_list) / len(Loc_list)
    return avg_Loc

# List to store localization lengths for different configurations
Loc_con_list = []
for _ in range(15):
    H = generate_hamiltonian(L, W)
    energy_levels, eigenstates = eigh(H)
    eigenstates2 = eigenstates**2
    Loc_con = cal_loc_length(eigenstates2)
    Loc_con_list.append(Loc_con)

# Calculate the mean of the localization lengths
average_loc_length = sum(Loc_con_list) / len(Loc_con_list)
print(f"Average Localization Length: {average_loc_length}")
