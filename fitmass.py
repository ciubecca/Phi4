import numpy as np
from time import time
from numpy import exp
import os
from integrator import *
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy

x = np.linspace(1,8,20)
y = np.loadtxt("g2.txt")

print(y)

imin = 5
m, b = np.polyfit(x[imin:], y[imin:], deg=1)

norm = 1/(12*(4*pi)**2)

print("norm = {}, m/norm = {}, m = {}, b= {}, b/norm ={},\
        e^(b/norm) = {}".format(norm,m/norm,m,b,
    b/norm, exp(b/norm)))

plt.scatter(x, y)

xx = np.linspace(min(x), max(x), 100)
plt.plot(xx, b + m*xx)

plt.ylim(min(y), max(y))

plt.xlabel(r"$\log \Lambda$")
plt.ylabel(r"$g_2$")
plt.savefig("g2.pdf")
