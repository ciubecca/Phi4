import numpy as np
from time import time
import os
from numpy import log, e
from integrator import *

nitn = 16
neval = 50000


# print("Computing O(VV) vacuum diagram...")
print("Computing O(g^2) vacuum energy with m=0 ...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

lamlist = np.linspace(10,200,30)
reslist = []

integ = Phi0_1_m0(nitn, neval)

for lam in lamlist:
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam)
    print("Result: {}".format(res))
    reslist.append(res.mean)

end = time()
print("Time passed: {}".format(end-start))

np.savetxt("Vac.txt", np.array(reslist))
