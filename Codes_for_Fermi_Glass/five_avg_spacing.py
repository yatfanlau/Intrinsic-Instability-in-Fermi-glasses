from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from random import randint
from scipy.optimize import curve_fit

# This code computes the average level spacing(<ΔE>) within a region of size L and see how <ΔE> changes with L

L = 60     #int(input('System size: ')) 
W = 8      #float(input('Disorder strength: '))
tp = 0.4

Loc = 7.516045912

def generate_disorder(L,W):
    disorder=W*((np.random.uniform(size=L*L)).reshape((L,L))-0.5)
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
  



def findmax(index):
    eigenstate = []
    for i in range(L*L):
        eigenstate.append(eigenstates2[i,int(index)])
    ep = np.argsort(eigenstate)
    position = ep[L*L-1]+1
    if position%L==0:
        x,y = L,(position//L)
    else:
        x,y = position%L,(position//L)+1
    return x,y

def d(x1,y1,x2,y2):
    dist = np.sqrt((x1-x2)**2+(y1-y2)**2)/Loc
    return dist

def Cal_Y(con):
    energy = [[]for i in range(10)]
    energy_spacing = [[]for j in range(10)]
    for m in range(10):
        Ld = X[m]
        for s in range(L*L):
            distance = d(30,30,findmax(s)[0],findmax(s)[1])
            if distance<Ld:
                energy[m].append(energy_levels[s])
        energy[m].sort()
        for k in range(len(energy[m])-1):
            level_spacing = abs(energy[m][k]-energy[m][k+1])
            energy_spacing[m].append(level_spacing)

    for p in range(10):
        avg = sum(energy_spacing[p])/len(energy_spacing[p])
        Y[con].append(avg)
    return Y[con]
        
# def func(x, A, B):
    # y = A/(B+x**2)
    # return y


# def fitting(con):
    # parameters, covariance = curve_fit(func, X, Y[con])
    # fit_A = parameters[0]
    # fit_B = parameters[1]
    # x_ = np.linspace(0.5,5,100)
    # y_ = fit_A/(fit_B+x_**2)
    # return x_,y_


Y = [[] for con in range(5)]
X_array = np.array([3,6,9,12,15,18,21,24,27,30]) # 10 numbers
X_array = X_array/Loc
X = X_array.tolist()

plot = ["^","s","o","*","D"]

for conf in range(5):
    H = generate_hamiltonian(L,W)
    (energy_levels,eigenstates) = linalg.eigh(H)
    eigenstates2 = eigenstates**2
    Y[conf] = Cal_Y(conf)
    #x_,y_ = fitting(conf)
    #plt.plot(x_,y_,label='y=A/(B+x^2)')
    plt.scatter(X, Y[conf],s=4,plot[conf])

plt.xlabel("Radius of the chosen region")
plt.ylabel("Average energy level spacing")
#plt.title(str(L)+'x'+str(L)+', '+'W='+str(W))
plt.legend()
plt.show()

