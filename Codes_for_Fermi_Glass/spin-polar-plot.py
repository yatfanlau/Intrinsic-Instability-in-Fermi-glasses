# N=1500, no_of_disorder=10, W=8, Number of Monte Carlo = 80000
# delta = 165 for U>=0.8
# delta = 150 for U<0.8

from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import proplot

def func(x, m, c):
    y = m*x+c
    return y

U_list = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8,1.0,1.2,1.4] #13

new_U = np.array(U_list)
new_U = 0.0925*new_U

no_of_spin_polar = [1.4, 2.6, 4.4, 7.0, 8.6, 9.8, 12.0, 13.4, 22.0, 32.0, 44.8, 61.4, 77.4]

fit_X = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]
new_fit_X = np.array(fit_X)
new_fit_X = 0.0925*new_fit_X
fit_Y = [1.4, 2.6, 4.4, 7.0, 8.6, 9.8, 12.0, 13.4]

#[3.8, 6.0, 11.0, 14.4, 17.6, 22.8, 27.0, 32.8]

parameters, covariance = curve_fit(func, new_fit_X, fit_Y)
fit_m = parameters[0]
fit_c = parameters[1]

print(fit_m,fit_c)
fig, ax = plt.subplots()
x_ = np.linspace(0,0.7*0.0925,100)
y_ = fit_m*x_+fit_c
plt.plot(x_,y_,linestyle=("solid"),color="b",label="Fitting result")
plt.scatter(new_U,no_of_spin_polar,s=15,color="r",label="Simulation result")
plt.xlabel(r'$UN(0)$',fontsize=12)
plt.ylabel("Number of spin polarized states",fontsize=12)
plt.legend()
#plt.grid(True)
#plt.axhline(0,ls='--',color='k')
#plt.axvline(0,ls='--',color='k')
fig.savefig('spin_new.png', dpi=300)
plt.show()
