import numpy as np
from time import time
import os
from numpy import log, e, pi
from integrator import *

nitn = 40
neval = 1000000


print("Computing O(g^2) order 1 correction to g2...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

# Need to suppress m/Lambda corrections
lam = 100000
reslist = []

integ = Phi1_1(nitn, neval)

coef = 1/(12*(4*pi)**2)

res = integ.do(lam)
print("Result: {}".format(res))
b = e**((res+coef*log(lam))/coef)
print("b = {}".format(b))

end = time()
print("Time passed: {}".format(end-start))
