import numpy as np
from time import time
import os
from integrator import *
import matplotlib.pyplot as plt
from matplotlib import rc


nitn = 20
neval = 50000

plt.style.use('ggplot')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


print("Computing O(VV) vacuum diagram...")

print("nitn={}, neval={}".format(nitn, neval))

start = time()

lamList = np.linspace(10,200,20)
reslist = []

integ = Phi0_1(nitn, neval)

for lam in lamList:
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam)
    res2 = res - integ.counterterm(lam)

    print("Result: {}".format(res))
    print("Result with counterterm: {}".format(res2))

    reslist.append(res2.mean)


end = time()
print("Time passed: {}".format(end-start))

plt.plot(lamList, reslist)
plt.xlabel(r"$\Lambda$")
plt.ylabel("$\mathcal{E}_0$")
plt.savefig("V2.pdf")
