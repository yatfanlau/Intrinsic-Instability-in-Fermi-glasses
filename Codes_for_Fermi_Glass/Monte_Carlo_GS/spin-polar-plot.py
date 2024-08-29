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

no_of_spin_polar = [0.6666666666666666, 2.533333333333333, 4.0, 6.133333333333334, 7.6, 8.933333333333334, 10.933333333333334, 12.8, 23.333333333333332, 35.2, 53.06666666666667, 68.8, 87.73333333333333]

fit_X = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]
fit_Y = [0.6666666666666666, 2.533333333333333, 4.0, 6.133333333333334, 7.6, 8.933333333333334, 10.933333333333334, 12.8]

#[3.8, 6.0, 11.0, 14.4, 17.6, 22.8, 27.0, 32.8]

parameters, covariance = curve_fit(func, fit_X, fit_Y)
fit_m = parameters[0]
fit_c = parameters[1]

print(fit_m,fit_c)
x_ = np.linspace(0,0.7,100)
y_ = fit_m*x_+fit_c
plt.plot(x_,y_,linestyle=("solid"),color="b",label="Fitting result")
plt.scatter(U_list,no_of_spin_polar,s=15,color="r",label="Simulation result")
plt.xlabel("Magnitude of U",fontsize=12)
plt.ylabel("Number of spin polarized states",fontsize=12)
plt.legend()
#plt.grid(True)
#plt.axhline(0,ls='--',color='k')
#plt.axvline(0,ls='--',color='k')
plt.show()
