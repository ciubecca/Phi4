import numpy as np
from time import time
import os
from integrator import *

nitn = 20
neval = 50000

# nitn = 10
# neval = 5000

print("Computing O(VV) vacuum diagram...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

lamList = [10, 20, 43, 63]

integ = Phi0_1(nitn, neval)

for lam in lamList:
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam)

    print("Result: {}".format(res))


end = time()
print("Time passed: {}".format(end-start))
