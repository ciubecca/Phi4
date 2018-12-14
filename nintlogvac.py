import numpy as np
from time import time
import os
from numpy import log, e
from integrator import *


nitn = 20
neval = 50000


# print("Computing O(VV) vacuum diagram...")
print("Computing O(g^2) logarithmic correction to vacuum...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

loglamlist = np.linspace(1,8,20)
reslist = []

integ = Phi0_1(nitn, neval)

for lam in e**(loglamlist):
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam) - integ.counterterm(lam)

    print("Result: {}".format(res))

    reslist.append(res.mean)


end = time()
print("Time passed: {}".format(end-start))

np.savetxt("logvac.txt", np.array(reslist))
