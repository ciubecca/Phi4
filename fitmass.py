import numpy as np
from time import time
from numpy import exp, e
import os
from integrator import *
import matplotlib.pyplot as plt
from matplotlib import rc
from paramplots import *
import scipy

# Need to fit m/Lambda corrections
x = np.linspace(0.01, 0.1, 20)
reslist = []

y = np.loadtxt("g2.txt")

m, b = np.polyfit(x, y, deg=1)

norm = 1/(12*(4*pi)**2)
# b = e**((res+coef*log(lam))/coef)

print("m={}, b={}, m/norm={}, e**(b/norm)={}".format(m,b,m/norm,e**(b/norm)))

plt.scatter(x, y)

xx = np.linspace(min(x), max(x), 100)
plt.plot(xx, b + m*xx)

plt.ylim(min(y), max(y))

plt.xlabel(r"$1/\Lambda$")
plt.title(r"$a = \frac{1}{12 (4 \pi)^2}$")
plt.ylabel(r"$g_2 - a \log (\Lambda)$")
plt.savefig("g2fitLambda.pdf")
