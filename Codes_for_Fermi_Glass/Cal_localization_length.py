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

#Calculate the average localization length for a L*L lattice

L = 60 #int(input('System size: ')) 
W = 9 #float(input('Disorder strength: '))
tp = 0.6

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

def cal_loc_length():
    Loc_list = []
    for k in range(1600,2000):
        IPR = np.dot(eigenstates2[:,k],eigenstates2[:,k])
        Loc = np.sqrt(1/IPR)
        Loc_list.append(Loc)
    avg_Loc = sum(Loc_list)/len(Loc_list)
    return avg_Loc
 

Loc_con_list = [] 
for con in range(15):
    H=generate_hamiltonian(L,W)
    (energy_levels,eigenstates)=linalg.eigh(H)
    eigenstates2 = eigenstates**2
    Loc_con = cal_loc_length()
    Loc_con_list.append(Loc_con)
LL = sum(Loc_con_list)/len(Loc_con_list)
print(LL)
    