import numpy as np
from time import time
import os
from numpy import log, e, pi
from integrator import *

nitn = 20
neval = 100000

print("Computing O(g^2) as function of 1/Lambda...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

# Need to fit m/Lambda corrections
invlamlist = np.linspace(0.01, 0.1, 20)
reslist = []

integ = Phi1_1(nitn, neval)

coef = 1/(12*(4*pi)**2)

for lam in 1/invlamlist:
    res = integ.do(lam)
    y = res+coef*log(lam)
    print(res, y)
    print("Result: {}, subtracted: {}".format(res, y))
    # b = e**((res+coef*log(lam))/coef)
    # print("b = {}".format(b))
    reslist.append(y.mean)

end = time()
print("Time passed: {}".format(end-start))

np.savetxt("g2.txt", np.array(reslist))
