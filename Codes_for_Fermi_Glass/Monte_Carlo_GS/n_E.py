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
import matplotlib.ticker as ticker
import scienceplots

# Set the style as 'science' without LaTeX
plt.style.use('science')
# Now adjust usetex option to False
plt.rcParams['text.usetex'] = False

L = 60
no_of_disorder = 15
no_of_U = 3                     #need filling


E_dict = {}
for con in range(no_of_disorder):
    E_dict["E"+str(con)] = np.load("E_"+str(con)+".npy")  #need filling

E_arr = np.zeros(L*L)
for con in range(no_of_disorder):
    E_arr = E_arr + E_dict["E"+str(con)]
E_arr = E_arr/no_of_disorder
print(E_arr[0])

U_list = [0.2,0.8,1.4]         #need filling

U_dict = {}

for u in range(no_of_U):
    U = U_list[u]
    for con in range(no_of_disorder):
        U_dict[str(U)+"_" + str( con ) ] = np.load(str(U)+"_"+str(con)+".npy")


Plot_y = [[]for i in range(no_of_U)]

for u in range(no_of_U):
    U_ = U_list[u]
    occ_arr = np.zeros((L*L,2))
    for con in range(no_of_disorder):
        occ_arr = occ_arr + U_dict[str(U_) + "_" + str(con)]
    occ_arr = occ_arr/no_of_disorder
    for t in range(L*L):
        Plot_y[u].append(occ_arr[t][0]+occ_arr[t][1])

def func(x,mu,beta):
    y = 2/(np.exp(beta*(x-mu))+1)
    return y
 
parameters0, covariance0 = curve_fit(func, E_arr, Plot_y[0])
fit_mu_0 = parameters0[0]
fit_beta_0 = parameters0[1]
print("U=0.2:",fit_mu_0,fit_beta_0)  
x_0 = np.linspace(-2,-0.4,5000)
y_0 = 2/(np.exp(fit_beta_0*(x_0-fit_mu_0))+1)
plt.plot(x_0,y_0,label='Circle: U=0.2')

parameters1, covariance1 = curve_fit(func, E_arr, Plot_y[1])
fit_mu_1 = parameters1[0]
fit_beta_1 = parameters1[1]
print("U=0.8:",fit_mu_1,fit_beta_1)  
x_1 = np.linspace(-2,-0.4,5000)
y_1 = 2/(np.exp(fit_beta_1*(x_1-fit_mu_1))+1)
plt.plot(x_1,y_1,linestyle=':',label='Star: U=0.8')

parameters2, covariance2 = curve_fit(func, E_arr, Plot_y[2])
fit_mu_2 = parameters2[0]
fit_beta_2 = parameters2[1]
print("U=1.4:",fit_mu_2,fit_beta_2)  
x_2 = np.linspace(-2,-0.4,5000)
y_2 = 2/(np.exp(fit_beta_2*(x_2-fit_mu_2))+1)
plt.plot(x_2,y_2,linestyle='dashed',label='Square: U=1.4')
plt.legend()

marker_list = ['o','*','+']
for j in range(no_of_U):
    Ulab = U_list[j]
    plt.scatter(E_arr,Plot_y[j],s=2,label="U="+str(Ulab),marker=marker_list[j])
plt.xlabel("Energy (t)",fontsize=11)
plt.ylabel(r"$\langle n_{a\uparrow}+n_{a\downarrow} \rangle_{dis}$",fontsize=11,labelpad=0.5)
plt.xlim(-2, -0.6)
plt.savefig("ne_plot.png", dpi=300)
plt.show()


# for e in range(no_of_disorder):
    # E_arr = np.load("E"+str(e)+".npy")
    