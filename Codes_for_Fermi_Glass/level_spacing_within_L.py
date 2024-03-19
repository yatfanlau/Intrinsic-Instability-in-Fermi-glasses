from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from random import randint
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson
import proplot


L = 60 #int(input('System size: ')) 
W = 8 #float(input('Disorder strength: '))
Loc = 7.516045912
tp = 0.4

def generate_disorder(L,W):
  disorder=W*((np.random.uniform(size=L*L)).reshape((L,L))-0.5)
  return disorder
  
def generate_hamiltonian(L,W):
  H=np.zeros((L*L,L*L))
  disorder=generate_disorder(L,W)
  for i in range(L):
    ip1=(i+1)%L
    im1=(i-1)%L
    for j in range(L):
      H[i*L+j,i*L+j]=disorder[i,j]
      jp1=(j+1)%L
      jm1=(j-1)%L
      H[ip1*L+j  ,i  *L+j  ]=1.0
      H[i*  L+j  ,ip1*L+j  ]=1.0
      H[i  *L+jp1,i  *L+j  ]=1.0
      H[i*  L+j  ,i  *L+jp1]=1.0
      
      H[i*L+j,ip1*L+jp1] = tp
      H[ip1*L+jp1,i*L+j] = tp
      H[i*L+j,ip1*L+jm1] = tp
      H[ip1*L+jm1,i*L+j] = tp
      
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
Ld = 1.5
energy_spacing = []
for con in range(15):
    H=generate_hamiltonian(L,W)
    (energy_levels,eigenstates)=linalg.eigh(H)
    eigenstates2 = eigenstates**2
    energy = []
    for s in range(L*L):
        distance = d(30,30,findmax(s)[0],findmax(s)[1])
        if distance<Ld:
            energy.append(energy_levels[s])
    energy.sort()
    for k in range(len(energy)-1):
        level_spacing = abs(energy[k]-energy[k+1])
        energy_spacing.append(level_spacing)

def func(x,h):
    y = (np.exp(-x/h))/h
    return y

BIN = int(np.sqrt(len(energy_spacing)))+1
n, binss, patches = plt.hist(energy_spacing, bins = BIN,density = True, color="b")

binss = binss.tolist()
binss.remove(binss[BIN])
binss = np.array(binss)
parameters, covariance = curve_fit(func, binss, n)
fit_h = parameters[0]
x_ = np.linspace(min(binss/fit_h),max(binss/fit_h),100)
y_ = (np.exp(-x_/fit_h))/fit_h
plt.plot(x_,y_,label="Poisson fitting "+r"$P(s)=\frac{1}{\Delta}e^{-s/\Delta}$",color="r",linestyle="solid")


# lam = np.mean(energy_spacing)
# x_plot = np.arange(0, max(energy_spacing) + 1)
# plt.plot(binss, poisson.pmf(binss, lam), label='Poisson fitting')


# plt.plot(bins, n)

plt.title('Grain size='+str(Ld*2)+r", $\Delta$="+str(round(fit_h,3)))
plt.xlabel("Energy level spacing"+r"/$\Delta$")
plt.ylabel("Normalized frequency")
plt.legend()
plt.show()



        
    