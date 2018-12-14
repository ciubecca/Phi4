import numpy as np
from time import time
import os
from numpy import log, e
from integrator import *
import matplotlib.pyplot as plt
from matplotlib import rc


nitn = 20
neval = 50000


plt.style.use('ggplot')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# print("Computing O(VV) vacuum diagram...")
print("Computing O(g^2) overlap ...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

loglamlist = np.linspace(1,8,20)
reslist = []

integ = Phi0_2(nitn, neval)

for lam in e**(loglamlist):
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam)

    print("Result: {}".format(res))

    reslist.append(res.mean)


end = time()
print("Time passed: {}".format(end-start))

np.savetxt("Overlap.txt", np.array(reslist))
