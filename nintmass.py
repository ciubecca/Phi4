import numpy as np
from time import time
import os
from numpy import log, e
from integrator import *

nitn = 16
neval = 10000


# print("Computing O(VV) vacuum diagram...")
print("Computing O(g^2) correction to g2...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

loglamlist = np.linspace(1,8,20)
reslist = []

integ = Phi1_1(nitn, neval)

for lam in e**(loglamlist):
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam)
    print("Result: {}".format(res))
    reslist.append(res.mean)
    print("Result - ct: {}".format(res-integ.counterterm(lam)))

end = time()
print("Time passed: {}".format(end-start))

np.savetxt("g2.txt", np.array(reslist))
