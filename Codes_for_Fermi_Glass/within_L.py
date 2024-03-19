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
import proplot

def func(x,h):
    y = (np.exp(-x/h))/h
    return y
Ld = 3
a = np.load("S=3.npy")
energy_spacing = a.tolist()

BIN = int(np.sqrt(len(energy_spacing)))+1
n, binss, batches = plt.hist(energy_spacing, bins = BIN,density = True,color="b") 

binss = binss.tolist()
binss.remove(binss[BIN])
binss = np.array(binss)
parameters, covariance = curve_fit(func, binss, n)
fit_h = parameters[0]
h = round(fit_h,2)

x_ = np.linspace(0,max(energy_spacing),100)
y_ = np.exp(-x_/fit_h)/fit_h
#plt.plot(x_,y_,label="Poisson fitting "+r"$P(s)=\frac{1}{\Delta}e^{-s/\Delta}$",color="r",linestyle="solid")
xm = 6*h 
ym = np.max(n)
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    xit = xm/5
    if N/(xit).is_integer():
        return str(np.round(N/h,2))
        
def format_func_y(value, tick_number):
    # find number of multiples of pi/2
    N = value
    yit = ym/5
    if N/(yit).is_integer():
        return str(np.round(N*h,2))

#plt.title('Grain size='+str(Ld)+r", $\Delta$="+str(round(fit_h,2)))
#plt.xlabel("Energy level spacing"+r"/$\Delta$")
#plt.ylabel("Normalized frequency")


ax = plt.axes()
ax.plot(x_,y_,label="Poisson fitting "+r"$P(s)=\frac{1}{\Delta}e^{-s/\Delta}$",color="r",linestyle="solid")
ax.set_title('Grain size='+str(Ld)+r", $\Delta$="+str(round(fit_h,2)))
ax.set_xlabel("Energy level spacing"+r"/$\Delta$")
ax.set_ylabel("Normalized frequency")
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func) )
ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func_y) )
plt.xlim(0,xm)
#ax.xaxis.set_major_locator(plt.MultipleLocator())
#plt.xlim([0,(8*h)])
#plt.xticks(np.linspace(0,2.5,5)/h)
#plt.xlim(0,5*h)
#xm = 1.5/h
#plt.xticks(np.linspace(0,xm,5)*h)
#plt.xticks(np.linspace(0,xm,5))



# #Changing scales for x-axis and y-axis
# locs, labels = plt.xticks()
# labels = [((item)*(1/h)) for item in locs]
# plt.xticks(locs, labels)
# locsy, labelsy = plt.yticks()
# labelsy = [((itemy)*(h)) for itemy in locs]
# plt.yticks(locsy, labelsy)

  
plt.legend()
plt.show()