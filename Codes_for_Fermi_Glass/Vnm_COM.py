from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from random import randint

# Within an energy range, randomly find N pairs to calculate Vnm vs distance, 
# where distance beyond half of the period is not considered.
# Use "centre of mass" definition to find the 

L = 60 #int(input('System size: ')) 
W = 10 #float(input('Disorder strength: '))

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


def d(x1,y1,x2,y2):
    dist = np.sqrt((x1-x2)**2+(y1-y2)**2)
    return dist
def V(p,q):
    V = np.dot(eigenstates2[:,p],eigenstates2[:,q])
    return V
    

dlist = []
Vlist = []
delta = 100
E1 = 3400-delta
E2 = 3400+delta

N = 30
for con in range(10):
    H=generate_hamiltonian(L,W)
    (energy_levels,eigenstates) = linalg.eigh(H)
    eigenstates2 = eigenstates**2
    for j in range(N):
        k = randint(E1,E2)
        m = randint(E1,E2)
        if k != m:
            distance = d(COM(k)[0],COM(k)[1],COM(m)[0],COM(m)[1])/6
            if distance < 5:
                Vnm = V(k,m)
                dlist.append(distance)
                Vlist.append(Vnm)
    

d_array = np.array(dlist)
V_array = np.array(Vlist)
ln_V_array = np.log(V_array)

A,B = np.polyfit(d_array, ln_V_array, 1, w=np.sqrt(V_array))



        
    
x = np.linspace(0,5,100)
y = np.exp(B)* np.exp(A*x)
plt.plot(x,y,label='U_kp = exp('+str(round(B,2))+')'+'exp('+str(round(A,2))+'d'+')')          
plt.scatter(dlist, Vlist,s=3,color='r')
plt.xlabel("Distance between two states(In units of localization length)")
plt.ylabel("U_kp")
plt.title('No. of total pairs = '+str(len(dlist))+', '+str(L)+'x'+str(L)+", W="+str(W)+', '+'no. of disorder confuguration = 10')
plt.legend()
plt.show()