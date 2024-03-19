from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from random import randint
import time
from numba import jit

L = 60  # int(input('System size: '))
W = 5  # float(input('Disorder strength: '))
Number = 10

delta = 180
no_of_disorder = 8
N = 1500
tp = 0.4


def generate_disorder(L, W):
    disorder = W * ((np.random.uniform(size=L * L)).reshape((L, L)) - 0.5)
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

#@jit(nopython=True)
def energy_change_same(m,y1):
    sum1,sum2=0,0
    for l in range(int(N/2+delta+2)):
        if n[l][y1]!=0:
            sum1 = sum1 + U*np.dot(eigenstates2[:,m],eigenstates2[:,l])*n[l][y1]  #for new n
    for l in range(int(N/2+delta+2)):
        if n[l][(y1+1)%2]!=0:
            sum2 = sum2 + U*np.dot(eigenstates2[:,m],eigenstates2[:,l])*n[l][(y1+1)%2]  #for old n
    return sum1-(sum2-U*np.dot(eigenstates2[:,m],eigenstates2[:,m]))


#@jit(nopython=True)
def energy_change(k,p,x1,y1,x0,y0):
    sum1,sum2,sum_k,sum_p = 0,0,0,0
    if x1==0:  #k is old
        if y0==y1:
            for l in range(int(N/2+delta+2)):
                if n[l][(y1+1)%2] != 0:
                    sum_k = sum_k + U*np.dot(eigenstates2[:,k],eigenstates2[:,l])*n[l][(y1+1)%2]  
            for l in range(int(N/2+delta+2)):
                if n[l][(y0+1)%2] != 0:
                    sum_p = sum_p + U*np.dot(eigenstates2[:,p],eigenstates2[:,l])*n[l][(y0+1)%2] 
            total_change = ((energy_levels[p]-energy_levels[k]) + sum_p - sum_k)  
        elif y0!=y1:
            for l in range(int(N/2+delta+2)):
                if n[l][(y1+1)%2] != 0:
                    sum_k = sum_k + U*np.dot(eigenstates2[:,k],eigenstates2[:,l])*n[l][(y1+1)%2]  
            for l in range(int(N/2+delta+2)):
                if n[l][(y0+1)%2] != 0:
                    sum_p = sum_p + U*np.dot(eigenstates2[:,p],eigenstates2[:,l])*n[l][(y0+1)%2] 
            total_change = ((energy_levels[p]-energy_levels[k]) + sum_p - (sum_k-U*np.dot(eigenstates2[:,k],eigenstates2[:,p])))
            
    elif x1==1:  #p is old
        if y0==y1:
            for l in range(int(N/2+delta+2)):
                if n[l][(y0+1)%2] != 0:
                    sum_k = sum_k + U*np.dot(eigenstates2[:,k],eigenstates2[:,l])*n[l][(y0+1)%2]  
            for l in range(int(N/2+delta+2)):
                if n[l][(y1+1)%2] != 0:
                    sum_p = sum_p + U*np.dot(eigenstates2[:,p],eigenstates2[:,l])*n[l][(y1+1)%2] 
            total_change = ((energy_levels[k]-energy_levels[p]) + sum_k - sum_p)  
        elif y0!=y1:
            for l in range(int(N/2+delta+2)):
                if n[l][(y0+1)%2] != 0:
                    sum_k = sum_k + U*np.dot(eigenstates2[:,k],eigenstates2[:,l])*n[l][(y0+1)%2]  
            for l in range(int(N/2+delta+2)):
                if n[l][(y1+1)%2] != 0:
                    sum_p = sum_p + U*np.dot(eigenstates2[:,p],eigenstates2[:,l])*n[l][(y1+1)%2] 
            total_change = ((energy_levels[k]-energy_levels[p]) + sum_k - (sum_p-U*np.dot(eigenstates2[:,k],eigenstates2[:,p])))       
    return total_change


U_list = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8,1.0,1.2,1.4] #13
spin_polar = [[] for i in range(13)]

done = 0
for con in range(no_of_disorder):
    H = generate_hamiltonian(L, W)
    (energy_levels, eigenstates) = linalg.eigh(H)
    eigenstates2 = eigenstates ** 2         
    for u in range(13):
        U=U_list[u]
        count = 0   
        n = np.zeros((L*L,2))
        for state in range(int(N / 2)):
            for spin in range(2): 
                n[state][spin] = 1
        for sample in range(Number):
            k,p = np.random.randint(int(N / 2 - delta),int(N / 2 + delta)),np.random.randint(int(N / 2 - delta),int(N / 2 + delta))
            #print(k,p)
            if k != p and n[k][0]+n[k][1]+n[p][0]+n[p][1] != 4 and n[k][0]+n[k][1]+n[p][0]+n[p][1] != 0:
                A = np.array([[n[k][0],n[k][1]],[n[p][0],n[p][1]]])
                q1,q0 = np.where(A==1),np.where(A==0)  #Positions of 1 and 0
                r1,r0 = np.random.randint(0,len(q1[0])),np.random.randint(0,len(q0[0])) # Randomly choose a 1 and 0
                x1,y1 = q1[0][r1],q1[1][r1]      #position of the randomly chosen 1
                x0,y0 = q0[0][r0],q0[1][r0]      #position of the randomly chosen 0
                A[x0][y0] = 1
                A[x1][y1] = 0
                n[k][0],n[k][1],n[p][0],n[p][1] = A[0][0],A[0][1],A[1][0],A[1][1]     
                if x0 == 0 and x1 == 0:
                    Delta_E = energy_change_same(k,y1)
                elif x0 == 1 and x1 == 1:
                    Delta_E = energy_change_same(p,y1)
                else:
                    Delta_E = energy_change(k,p,x1,y1,x0,y0)
                if Delta_E >= 0:
                    A[x0][y0] = 0
                    A[x1][y1] = 1
                    n[k][0],n[k][1],n[p][0],n[p][1] = A[0][0],A[0][1],A[1][0],A[1][1]
                elif Delta_E < 0:
                    count = count + 1
        print(count)
        print(n[int(N/2-delta):int(N/2+delta)])
        num = 0
        for l in range(int(N/2-delta),int(N/2+delta)):
            if n[l][0]+n[l][1] == 1:
                num = num + 1
        spin_polar[u].append(num)  
        done = done + 1
        print('finish',done)

Avg_spin_polar = []
for j in range(13):
    Avg = sum(spin_polar[j])/len(spin_polar[j])
    Avg_spin_polar.append(Avg)

plt.scatter(U_list,Avg_spin_polar,s=5)
plt.xlabel("Manitude of U")
plt.ylabel("Number of spin polarized states")
plt.title(str(60)+'x'+str(60)+', '+'W='+str(5))
plt.show()    


 
    
    

