import math
import random
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from random import randint

 

L = 60 #int(input('System size: ')) 
W = 14 #float(input('Disorder strength: '))

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

def COM(index):
    s,r = 0,0
    for i in range(L*L):
        ip1 = i + 1
        if ip1 % L == 0:
            x,y = L, (i//L)
        else:
            x,y = i%L,(i//L)+1
        s = s + x*eigenstates2[i,index]
        r = r + y*eigenstates2[i,index]
    return s,r

def findmax(index):
    eigenstate = []
    for i in range(L*L):
        eigenstate.append(eigenstates2[i,int(index)])
    ep = np.argsort(eigenstate)
    position = ep[L*L-1]+1
    x,y = position%L,(position//L)+1
    return x,y

H = generate_hamiltonian(L, W)
(energy_levels, eigenstates) = linalg.eigh(H)
eigenstates2 = eigenstates ** 2

for i in range(20):
    a = randint(3000,3400)
    print(COM(a))
    print(findmax(a))