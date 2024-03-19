from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from random import randint
import proplot
from scipy.optimize import curve_fit

L = 60 #int(input('System size: ')) 
W = 9 #float(input('Disorder strength: '))
Loc = 6.02516045912
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
def V(p,q):
    V = np.dot(eigenstates2[:,p],eigenstates2[:,q])
    return V

def idk(Ld):
    Vnm_list = []
    for n in range(5000):
        i = randint(E1,E2)
        j = randint(E1,E2)
        if i!=j:
            distance = d(findmax(i)[0],findmax(i)[1],findmax(j)[0],findmax(j)[1])
            for s in range(5):
                if (d_list[s]-Range)<distance<(d_list[s]+Range):
                    Vnm = V(i,j)
                    Vnm_list.append(Vnm)
    avg_Vnm = sum(Vnm_list)/len(Vnm_list)
    return avg_Vnm

d_list = np.array([7,12,17,22,27])
d_list = d_list/Loc

Range = 0.1
delta = 200
E1 = 1800-delta  #1800,3000
E2 = 1800+delta

disorder_Vnm = [[]for i in range(5)]
final_Vnm = []
count=0


##### Raw data plot
Data_Vnm = []
Distance_list = []
#### 
 
for con in range(5):
    H=generate_hamiltonian(L,W)
    (energy_levels,eigenstates) = linalg.eigh(H)
    eigenstates2 = eigenstates**2
    for n in range(2000):
        i = randint(E1,E2)
        j = randint(E1,E2)
        if i!=j:
            distance = d(findmax(i)[0],findmax(i)[1],findmax(j)[0],findmax(j)[1])
            if distance < 27/Loc:
                Distance_list.append(distance)
                if (d_list[0]-Range) < distance < (d_list[0]+Range):
                    disorder_Vnm[0].append(V(i,j))
                    Data_Vnm.append(V(i,j))
                elif (d_list[1]-Range) < distance < (d_list[1]+Range):
                    disorder_Vnm[1].append(V(i,j))
                    Data_Vnm.append(V(i,j))
                elif (d_list[2]-Range) < distance < (d_list[2]+Range):
                    disorder_Vnm[2].append(V(i,j))
                    Data_Vnm.append(V(i,j))
                elif (d_list[3]-Range) < distance < (d_list[3]+Range):
                    disorder_Vnm[3].append(V(i,j))
                    Data_Vnm.append(V(i,j))
                elif (d_list[4]-Range) < distance < (d_list[4]+Range):
                    disorder_Vnm[4].append(V(i,j))
                    Data_Vnm.append(V(i,j))
                else:
                    Data_Vnm.append(V(i,j))
    count = count + 1
    print(count)
for k in range(5):
    temp_Vnm = sum(disorder_Vnm[k])/len(disorder_Vnm[k])
    final_Vnm.append(temp_Vnm)


def func(x, m, c):
    y = m*x+c
    return y


    
Vnm_ = np.array(final_Vnm)
Vnm_ = np.log(Vnm_)

parameters, covariance = curve_fit(func, d_list,Vnm_)
fit_m = parameters[0]
fit_c = parameters[1]
print(fit_m,fit_c)
x_ = np.linspace(0,27/Loc,100)
y_ = fit_m*x_+fit_c
plt.plot(x_,y_,linestyle=("solid"),color="b",label="Fitting result: y=mx+c")
plt.scatter(d_list,Vnm_,s=10,color="red",label="Simulation result")
plt.xlabel(r"$Distance$(in unit of localization length)",fontsize=10)
plt.ylabel(r"$ln(U_{kp})$",fontsize=10)
plt.title("m="+str(round(fit_m,2))+", c="+str(round(fit_c,2)))
plt.legend()
plt.show()
       
########## Raw data plot
d_array = np.array(Distance_list)
V_array = np.array(Data_Vnm)
ln_V_array = np.log(V_array)

A,B = np.polyfit(d_array, ln_V_array, 1, w=np.sqrt(V_array))

  
x = np.linspace(0,27/Loc,100)
y = np.exp(B)* np.exp(A*x)
print(A,B)
plt.plot(x,y,linestyle="solid",color="r",label=r'$U_{kp}$='+str(round(np.exp(B),2))+'exp('+str(round(A,2))+'d'+')')          
plt.scatter(Distance_list, Data_Vnm, s=3,color='b')
plt.xlabel(r"$d$",fontsize=10)
plt.ylabel(r"$U_{kp}$",fontsize=10)
#plt.title('No. of total pairs = '+str(len(dlist))+', '+str(L)+'x'+str(L)+", W="+str(W)+', '+'no. of disorder confuguration = 10')
plt.legend()
plt.show() 
###############




    